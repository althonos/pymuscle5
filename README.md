# pyMUSCLE5 [![Stars](https://img.shields.io/github/stars/althonos/pymuscle.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pymuscle/stargazers)

*Cython bindings and Python interface to [MUSCLE v5](https://www.drive5.com/muscle/), a highly efficient and accurate multiple sequence alignment software.*

[![Actions](https://img.shields.io/github/workflow/status/althonos/pymuscle5/Test/main?logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pymuscle5/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pymuscle5?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pymuscle5/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pymuscle5.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pymuscle5)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pymuscle5?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pymuscle5)
[![AUR](https://img.shields.io/aur/version/python-pymuscle?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pymuscle5)
[![Wheel](https://img.shields.io/pypi/wheel/pymuscle5.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pymuscle5/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pymuscle5.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pymuscle5/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pymuscle5.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pymuscle5/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pymuscle5/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pymuscle5.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pymuscle5/issues)
[![Docs](https://img.shields.io/readthedocs/pymuscle5/latest?style=flat-square&maxAge=600)](https://pymuscle5.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pymuscle5/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpymuscle5)](https://pepy.tech/project/pymuscle5)

## 🗺️ Overview

MUSCLE is widely-used software for making multiple alignments of biological
sequences. Version 5 of MUSCLE achieves highest scores on several benchmark
tests and scales to thousands of sequences on a commodity desktop computer.

pyMUSCLE5 is a Python module that provides bindings to MUSCLE v5 using
[Cython](https://cython.org/). It directly interacts with the MUSCLE
internals, which has the following advantages:

- **single dependency**: If your software or your analysis pipeline is
  distributed as a Python package, you can add `pymuscle5` as a dependency to
  your project, and stop worrying about the MUSCLE binaries being properly
  setup on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you fully control, so you don't have to invoke the MUSCLE CLI using a
  sub-process and temporary files. Sequences can be passed directly as
  strings or bytes, which avoids the overhead of formatting your input to
  FASTA for MUSCLE.
- **no OpenMP-dependency**: The original MUSCLE code uses [OpenMP](https://www.openmp.org/)
  to parallelize embarassingly-parallel tasks. In pyMUSCLE5 the dependency on
  OpenMP has been removed in favor of the Python `threading` module for better
  portability.

## 💡 Example

Let's load some sequences sequence from a FASTA file, use an `Aligner` to
align proteins together, and print the alignment in two-line FASTA format.

### 🔬 [Biopython](https://github.com/biopython/biopython)

```python
import os

import Bio.SeqIO
import pymuscle

path = os.path.join("pymuscle", "tests", "data", "swissprot-halorhodopsin.faa")
records = list(Bio.SeqIO.parse(path, "fasta"))

sequences = [
    pymuscle.Sequence(record.id.encode(), bytes(record.seq))
    for record in records
]

aligner = pymuscle.Aligner()
msa = aligner.align(sequences)

for seq in msa.sequences:
    print(f">{seq.name.decode()}")
    print(seq.sequence.decode())
```

### 🧪 [Scikit-bio](https://github.com/biocore/scikit-bio)

```python
import os

import skbio.io
import pymuscle

path = os.path.join("pymuscle", "tests", "data", "swissprot-halorhodopsin.faa")
records = list(skbio.io.read(path, "fasta"))

sequences = [
    pymuscle.Sequence(record.metadata["id"].encode(), record.values.view('B'))
    for record in records
]

aligner = pymuscle.Aligner()
msa = aligner.align(sequences)

for seq in msa.sequences:
    print(f">{seq.name.decode()}")
    print(seq.sequence.decode())
```

*We need to use the [`view`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.view.html)
method to get the sequence viewable by Cython as an array of `unsigned char`.*


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/pymuscle/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pymuscle/blob/main/CONTRIBUTING.md)
for more details.

## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pymuscle/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
The MUSCLE code was written by [Robert Edgar](https://github.com/rcedgar) and is distributed under the
terms of the GPLv3 as well. See `vendor/muscle/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original MUSCLE authors](https://github.com/rcedgar). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
