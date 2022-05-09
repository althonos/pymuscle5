import configparser
import glob
import os
import platform
import re
import subprocess
import sys

import setuptools
import setuptools.extension
from distutils import log
from distutils.errors import CompileError
from distutils.command.clean import clean as _clean
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Constants -----------------------------------------------------------------

SETUP_FOLDER = os.path.realpath(os.path.join(__file__, os.pardir))

SYSTEM  = platform.system()
MACHINE = platform.machine()
if re.match("^mips", MACHINE):
    TARGET_CPU = "mips"
elif re.match("^arm", MACHINE):
    TARGET_CPU = "arm"
elif re.match("^aarch64", MACHINE):
    TARGET_CPU = "aarch64"
elif re.match("(x86_64)|(AMD64|amd64)|(^i.86$)", MACHINE):
    TARGET_CPU = "x86"
elif re.match("^(powerpc|ppc)", MACHINE):
    TARGET_CPU = "ppc"

# --- Utils ------------------------------------------------------------------

def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

# --- Commands ------------------------------------------------------------------

class Extension(setuptools.extension.Extension):

    def __init__(self, *args, **kwargs):
        self._needs_stub = False
        self.platform_sources = kwargs.pop("platform_sources", {})
        super().__init__(*args, **kwargs)


# --- Commands ------------------------------------------------------------------

class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    # --- Compatibility with `setuptools.Command`

    user_options = _build_ext.user_options + [
        ("disable-avx2", None, "Force compiling the extension without AVX2 instructions"),
        ("disable-sse2", None, "Force compiling the extension without SSE2 instructions"),
        ("disable-neon", None, "Force compiling the extension without NEON instructions"),
    ]

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.disable_avx2 = False
        self.disable_sse2 = False
        self.disable_neon = False

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # record SIMD-specific options
        self._simd_supported = dict(AVX2=False, SSE2=False, NEON=False)
        self._simd_defines = dict(AVX2=[], SSE2=[], NEON=[])
        self._simd_flags = dict(AVX2=[], SSE2=[], NEON=[])
        self._simd_disabled = {
            "AVX2": self.disable_avx2,
            "SSE2": self.disable_sse2,
            "NEON": self.disable_neon,
        }
        # transfer arguments to the build_clib method
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.compiler = self.compiler

    # --- Autotools-like helpers ---

    def _check_simd_generic(self, name, flags, header, vector, set, op, extract):
        _eprint('checking whether compiler can build', name, 'code', end="... ")

        base = "have_{}".format(name)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write("""
                #include <{}>
                int main() {{
                    {}      a = {}(1);
                            a = {}(a);
                    short   x = {}(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """.format(header, vector, set, op, extract))

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_preargs=flags)
            self.compiler.link_executable(objects, base, output_dir=self.build_temp)
            subprocess.run([binfile], check=True)
        except CompileError:
            _eprint("no")
            return False
        except subprocess.CalledProcessError:
            _eprint("yes, but cannot run code")
            return True  # assume we are cross-compiling, and still build
        else:
            if not flags:
                _eprint("yes")
            else:
                _eprint("yes, with {}".format(" ".join(flags)))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _avx2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:AVX2"]
        return ["-mavx", "-mavx2"]

    def _check_avx2(self):
        return self._check_simd_generic(
            "AVX2",
            self._avx2_flags(),
            header="immintrin.h",
            vector="__m256i",
            set="_mm256_set1_epi16",
            op="_mm256_abs_epi32",
            extract="_mm256_extract_epi16",
        )

    def _sse2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:SSE2"]
        return ["-msse", "-msse2"]

    def _check_sse2(self):
        return self._check_simd_generic(
            "SSE2",
            self._sse2_flags(),
            header="emmintrin.h",
            vector="__m128i",
            set="_mm_set1_epi16",
            op="_mm_move_epi64",
            extract="_mm_extract_epi16",
        )

    def _neon_flags(self):
        return ["-mfpu=neon"] if TARGET_CPU == "arm" else []

    def _check_neon(self):
        return self._check_simd_generic(
            "NEON",
            self._neon_flags(),
            header="arm_neon.h",
            vector="int16x8_t",
            set="vdupq_n_s16",
            op="vabsq_s16",
            extract="vgetq_lane_s16"
        )

    # --- Build code ---

    def build_simd_code(self, ext):
        # build platform-specific code
        for simd, sources in ext.platform_sources.items():
            if self._simd_supported[simd] and not self._simd_disabled[simd]:
                objects = [
                    os.path.join(self.build_temp, s.replace(".c", self.compiler.obj_extension))
                    for s in sources
                ]
                for source, object in zip(sources, objects):
                    self.make_file(
                        [source],
                        object,
                        self.compiler.compile,
                        (
                            [source],
                            self.build_temp,
                            ext.define_macros + self._simd_defines[simd],
                            ext.include_dirs,
                            self.debug,
                            ext.extra_compile_args + self._simd_flags[simd],
                            None,
                            ext.depends
                        )
                    )
                ext.extra_objects.extend(objects)

    def build_extension(self, ext):
        # show the compiler being used
        _eprint("building", ext.name, "with", self.compiler.compiler_type, "compiler")
        # add debug symbols if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Z7")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))
        # add C++11 flags
        if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
            ext.extra_compile_args.append("-std=c++11")
        elif self.compiler.compiler_type == "msvc":
            ext.extra_compile_args.append("/std:c11")
        # update link and include directories
        for name in ext.libraries:
            lib = self._clib_cmd.get_library(name)
            libfile = self.compiler.library_filename(
                lib.name, output_dir=self._clib_cmd.build_clib
            )
            ext.depends.append(libfile)
            ext.extra_objects.append(libfile)
        # build platform-specific code
        self.build_simd_code(ext)
        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)

    def build_extensions(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "compiler_directives": {
                "cdivision": True,
                "nonecheck": False,
            },
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
                "TARGET_CPU": TARGET_CPU,
                "AVX2_BUILD_SUPPORT": False,
                "NEON_BUILD_SUPPORT": False,
                "SSE2_BUILD_SUPPORT": False,
            }
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["cdivision_warnings"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False

        # compile the C library
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # check if we can build platform-specific code
        # if TARGET_CPU == "x86":
        #     if not self._simd_disabled["AVX2"] and self._check_avx2():
        #         cython_args["compile_time_env"]["AVX2_BUILD_SUPPORT"] = True
        #         self._simd_supported["AVX2"] = True
        #         self._simd_flags["AVX2"].extend(self._avx2_flags())
        #         self._simd_defines["AVX2"].append(("__AVX2__", 1))
        #     if not self._simd_disabled["SSE2"] and self._check_sse2():
        #         cython_args["compile_time_env"]["SSE2_BUILD_SUPPORT"] = True
        #         self._simd_supported["SSE2"] = True
        #         self._simd_flags["SSE2"].extend(self._sse2_flags())
        #         self._simd_defines["SSE2"].append(("__SSE2__", 1))
        # elif TARGET_CPU == "arm" or TARGET_CPU == "aarch64":
        #     if not self._simd_disabled["NEON"] and self._check_neon():
        #         cython_args["compile_time_env"]["NEON_BUILD_SUPPORT"] = True
        #         self._simd_supported["NEON"] = True
        #         self._simd_flags["NEON"].extend(self._neon_flags())
        #         self._simd_defines["NEON"].append(("__ARM_NEON__", 1))

        # add the include dirs
        for ext in self.extensions:
            ext.include_dirs.append(self._clib_cmd.build_clib)

        # cythonize the extensions (retaining platform-specific sources)
        platform_sources = [ext.platform_sources for ext in self.extensions]
        self.extensions = cythonize(self.extensions, **cython_args)
        for ext, plat_src in zip(self.extensions, platform_sources):
            ext.platform_sources = plat_src

        # build the extensions as normal
        _build_ext.build_extensions(self)


class build_clib(_build_clib):
    """A custom `build_clib` that splits the `training.c` file from Prodigal.
    """

    # --- Autotools-like helpers ---

    def _check_function(self, funcname, header, args="()"):
        _eprint('checking whether function', repr(funcname), 'is available', end="... ")

        base = "have_{}".format(funcname)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        with open(testfile, "w") as f:
            f.write("""
                #include <{}>
                int main() {{
                    {}{};
                    return 0;
                }}
            """.format(header, funcname, args))
        try:
            objects = self.compiler.compile([testfile], debug=self.debug)
            self.compiler.link_executable(objects, base, output_dir=self.build_temp)
        except CompileError:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _publicize(self, input, output):
        with open(input, "rb") as src:
            with open(output, "wb") as dst:
                for line in src:
                    dst.write(line.replace(b"private:", b"public:"))

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    def get_library(self, name):
        return next(lib for lib in self.libraries if lib.name == name)

    # --- Build code ---

    def build_libraries(self, libraries):
        # check for functions required for libcpu_features on OSX
        # if SYSTEM == "Darwin":
        #     if self._check_function("sysctlbyname", "sys/sysctl.h", args="(NULL, NULL, 0, NULL, 0)"):
        #         self.compiler.define_macro("HAVE_SYSCTLBYNAME", 1)

        # build each library only if the sources are outdated
        self.mkpath(self.build_clib)
        for library in libraries:
            libname = self.compiler.library_filename(library.name, output_dir=self.build_clib)
            self.make_file(library.sources, libname, self.build_library, (library,))

    def build_library(self, library):
        # show the compiler being used
        _eprint("building", library.name, "with", self.compiler.compiler_type, "compiler")

        # add debug flags if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Z7")

        # add C++11 flags
        if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
            library.extra_compile_args.append("-std=c++11")
        elif self.compiler.compiler_type == "msvc":
            library.extra_compile_args.append("/std:c11")

        # expose all private members and copy headers to build directory
        for header in library.depends:
            output = os.path.join(self.build_clib, os.path.basename(header))
            self.make_file(
                [header],
                output,
                self._publicize,
                (header, output)
            )

        # copy sources to build directory
        sources = [
            os.path.join(self.build_temp, os.path.basename(source))
            for source in library.sources
        ]
        for source, source_copy in zip(library.sources, sources):
            self.make_file(
                [source],
                source_copy,
                self.copy_file,
                (source, source_copy)
            )

        # store compile args
        compile_args = (
            library.define_macros,
            library.include_dirs,
            self.debug,
            library.extra_compile_args,
            None,
            library.depends,
        )
        # manually prepare sources and get the names of object files
        objects = [
            s.replace(".cpp", self.compiler.obj_extension)
            for s in sources
        ]
        # only compile outdated files
        for source, object in zip(sources, objects):
            self.make_file(
                [source],
                object,
                self.compiler.compile,
                ([source], None, *compile_args),
            )

        # link into a static library
        libfile = self.compiler.library_filename(
            library.name,
            output_dir=self.build_clib,
        )
        self.make_file(
            objects,
            libfile,
            self.compiler.create_static_lib,
            (objects, library.name, self.build_clib, None, self.debug)
        )


class clean(_clean):
    """A `clean` that removes intermediate files created by Cython.
    """

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pymuscle")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c", "*.cpp"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                log.info("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)

# --- Setup ---------------------------------------------------------------------

setuptools.setup(
    libraries=[
        Library(
            "muscle",
            language="c++",
            sources=[
                x
                for x in glob.glob(os.path.join(SETUP_FOLDER, "vendor", "muscle", "src", "*.cpp"))
                if os.path.basename(x) not in {"main.cpp", "make_a2m.cpp"}
            ],
            depends=[
                x
                for x in glob.glob(os.path.join(SETUP_FOLDER, "vendor", "muscle", "src", "*.h"))
            ],
            include_dirs=[
                os.path.join(SETUP_FOLDER, "pymuscle"),
                os.path.join(SETUP_FOLDER, "include"),
                os.path.join(SETUP_FOLDER, "pymuscle", "patch"),
            ],
        ),
    ],
    ext_modules=[
        Extension(
            "pymuscle._pymuscle",
            language="c++",
            sources=[
                os.path.join(SETUP_FOLDER, "pymuscle", "_pymuscle.pyx"),
            ],
            platform_sources={
                # "AVX2": ["pymuscle/impl/avx.c"],
                # "NEON": ["pymuscle/impl/neon.c"],
                # "SSE2": ["pymuscle/impl/sse.c"],
            },
            include_dirs=[
                os.path.join(SETUP_FOLDER, "pymuscle"),
                os.path.join(SETUP_FOLDER, "include"),
                os.path.join(SETUP_FOLDER, "pymuscle", "patch"),
            ],
            libraries=[
                "muscle",
            ],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "build_clib": build_clib,
        "clean": clean
    }
)
