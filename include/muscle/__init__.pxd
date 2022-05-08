from libcpp cimport bool

cdef extern from "myutils.h" nogil:
    ctypedef unsigned char byte

cdef extern from "types.h" nogil:
    cdef enum LINKAGE:
        LINKAGE_Min
        LINKAGE_Max
        LINKAGE_Avg
        LINKAGE_Biased

# cdef extern from "myutils.h" nogil:
#     pass

cdef extern from * nogil:
    extern bool         opt_quiet
    extern unsigned int opt_threads
