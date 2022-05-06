from libcpp cimport bool

from muscle cimport byte


cdef extern from "alpha.h":

    cdef enum ALPHA:
        ALPHA_Undefined
        ALPHA_Nucleo
        ALPHA_Amino


    cdef bool StrHasAmino(const char* Str)
    cdef bool StrHasGap(const char* Str)

    cdef void ClearInvalidLetterWarning()
    cdef void InvalidLetterWarning()
    cdef void ReportInvalidLetters()

    cdef extern unsigned g_CharToLetter[];
    cdef extern unsigned g_CharToLetterEx[];




    cdef void SetAlpha(ALPHA Alpha)
