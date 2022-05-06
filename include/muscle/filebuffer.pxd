from libcpp cimport bool





cdef extern from "filebuffer.h" nogil:

    cdef cppclass FileBuffer:
        FileBuffer(const char* filename)

        bool fail() const
