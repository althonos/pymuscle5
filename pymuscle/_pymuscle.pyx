# distutils: language = c++
# cython: language_level=3, linetrace=True

from libc.stdio cimport printf, fprintf, FILE, fopen, fclose
from libc.string cimport memcpy
from libcpp.vector cimport vector
from libcpp.string cimport string

from cpython.buffer cimport PyBUF_READ, PyBUF_FORMAT
from cpython.memoryview cimport PyMemoryView_FromMemory

from muscle.alpha cimport ALPHA, SetAlpha
from muscle.hmmparams cimport HMMParams
from muscle.sequence cimport Sequence as _Sequence
from muscle.msa cimport MSA as _MSA
from muscle.multisequence cimport MultiSequence as _MultiSequence
from muscle.treeperm cimport TREEPERM, TREEPERMToStr
from muscle.mpcflat cimport MPCFlat

import os
import threading
import warnings


cdef class Sequence:

    cdef _Sequence* _seq

    def __cinit__(self):
        self._seq = NULL

    def __dealloc__(self):
        _Sequence.DeleteSequence(self._seq)
        self._seq = NULL

    def __init__(
        self,
        bytes name not None,
        const unsigned char[::1] sequence not None,
    ):

        if self._seq != NULL:
            _Sequence.DeleteSequence(self._seq)
        self._seq = _Sequence.NewSequence()

        self._seq.m_GSI = 0
        self._seq.m_SMI = 0
        self._seq.m_Label = <string> name

        self._seq.m_CharVec = vector[char](len(sequence) + 1)
        cdef char* dst = self._seq.m_CharVec.data()
        dst[0] = b'@'
        memcpy(&dst[1], &sequence[0], len(sequence))

    def __copy__(self):
        return self.copy()

    def __len__(self):
        assert self._seq != NULL
        return self._seq.GetLength()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        assert self._seq is not NULL

        if flags & PyBUF_FORMAT:
            buffer.format = b"B"
        else:
            buffer.format = NULL
        buffer.buf = &self._seq.m_CharVec.data()[1]
        buffer.internal = NULL
        buffer.itemsize = sizeof(char)
        buffer.len = (len(self._seq.m_CharVec) - 1) * sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        # buffer.shape = self._shape
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def name(self):
        assert self._seq is not NULL
        return <bytes> self._seq.m_Label

    cpdef Sequence copy(self):
        assert self._seq != NULL
        cdef Sequence copy = Sequence.__new__(Sequence)
        copy._seq = self._seq


cdef class MultiSequence:

    cdef _MultiSequence* _mseq

    @classmethod
    def from_file(cls, filename):
        cdef bytes _filename = os.fsencode(filename)
        cdef MultiSequence obj = cls.__new__(cls)
        obj._mseq = new _MultiSequence(<string> _filename)
        return obj

    def __cinit__(self):
        self._mseq = NULL

    def __dealloc__(self):
        del self._mseq

    def __init__(self):
        raise NotImplementedError("MultiSequence.__init__")

    def __len__(self):
        assert self._mseq != NULL
        return self._mseq.m_Seqs.size()


cdef class Alignment:

    cdef _MSA _msa

    def __cinit__(self):
        self._msa = _MSA()


cdef class Aligner:

    cdef MPCFlat _mpcflat

    def __cinit__(self):
        self._mpcflat.Clear()

    def __init__(
        self,
        *,
        object consistency_iterations = None,
        object refine_iterations = None,
    ):
        self._mpcflat.Clear()

        if consistency_iterations is not None:
            self._mpcflat.m_ConsistencyIterCount = int(consistency_iterations)
        if refine_iterations is not None:
            self._mpcflat.m_RefineIterCount = int(refine_iterations)

        self._mpcflat.m_TreePerm = TREEPERM.TP_None

    cpdef object align(
        self,
        MultiSequence sequences,
        unsigned int seed = 0,
    ):
        self._mpcflat.Run(sequences._mseq)

        cdef Alignment msa = Alignment.__new__(Alignment)
        self._mpcflat.m_MSA.ToMSA(msa._msa)

        return msa
