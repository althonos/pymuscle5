
cdef extern from "myutils.h" nogil:
    ctypedef unsigned char byte

cdef extern from "types.h" nogil:
    cdef enum LINKAGE:
        LINKAGE_Min
        LINKAGE_Max
        LINKAGE_Avg
        LINKAGE_Biased
