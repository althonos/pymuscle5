# distutils: language = c++
# cython: language_level=3, linetrace=True

from libc.stdlib cimport malloc, free
from libc.stdio cimport printf, fprintf, FILE, fopen, fclose
from libc.string cimport memcpy, strlen
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
        buffer.internal = <void*> malloc(sizeof(Py_ssize_t))
        buffer.itemsize = sizeof(char)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = True
        buffer.shape = <Py_ssize_t*> buffer.internal
        buffer.len = (len(self._seq.m_CharVec) - 1) * sizeof(char)
        buffer.shape[0] = buffer.len
        buffer.suboffsets = NULL
        buffer.strides = NULL

    def __releasebuffer__(self, Py_buffer* buffer):
        free(buffer.internal)

    def __repr__(self):
        cdef type  ty  = type(self)
        cdef bytes seq = memoryview(self).tobytes()
        return f"{ty.__module__}.{ty.__name__}({self.name!r}, {seq!r})"

    @property
    def name(self):
        assert self._seq is not NULL
        return <bytes> self._seq.m_Label

    @property
    def sequence(self):
        assert self._seq is not NULL
        return <bytes> self._seq.m_CharVec

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


cdef class _AlignmentSequences:
    cdef Alignment _alignment

    def __cinit__(self, Alignment alignment):
        self._alignment = alignment

    def __len__(self):
        return self._alignment.m_uSeqCount

    def __getitem__(self, ssize_t i):

        cdef Sequence seq
        cdef ssize_t  index     = i
        cdef unsigned seq_count = self._alignment._msa.m_uSeqCount
        cdef unsigned col_count = self._alignment._msa.m_uColCount

        if index < 0:
            index += seq_count
        if index < 0 or index >= seq_count:
            raise IndexError(i)

        seq = Sequence.__new__(Sequence)
        seq._seq = _Sequence.NewSequence()
        seq._seq.m_Label = string(self._alignment._msa.m_szNames[index])

        seq._seq.m_CharVec.reserve(col_count + 1)
        seq._seq.m_CharVec[0] = b'@'
        memcpy(&seq._seq.m_CharVec[1], self._alignment._msa.m_szSeqs[index], col_count*sizeof(char))

        seq._seq.m_CharVec.assign(
            self._alignment._msa.m_szSeqs[index],
            self._alignment._msa.m_szSeqs[index] + col_count
        )

        return seq


cdef class Alignment:

    cdef _MSA _msa

    def __cinit__(self):
        self._msa = _MSA()

    def __copy__(self):
        return self.copy()

    @property
    def names(self):
        return tuple(
            <bytes> self._msa.m_szNames[i]
            for i in range(self._msa.m_uSeqCount)
        )

    @property
    def sequences(self):
        return _AlignmentSequences.__new__(_AlignmentSequences, self)

    cpdef Alignment copy(self):
        cdef Alignment copy = Alignment.__new__(Alignment)
        copy._msa.Copy(self._msa)
        return copy


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

    cdef _alloc_pair_count(self, unsigned int pair_count):
        if pair_count < self._mpcflat.m_ptrSparsePosts.size():
            return
        self._mpcflat.m_SparsePosts1.resize(pair_count)
        self._mpcflat.m_SparsePosts2.resize(pair_count)

    cdef _init_seqs(self, MultiSequence sequences):
        cdef unsigned int i
        cdef _Sequence*   seq
        cdef string       label
        cdef unsigned int seq_count = sequences._mseq.GetSeqCount()

        self._mpcflat.m_InputSeqs = sequences._mseq
        self._mpcflat.m_Labels.clear()
        self._mpcflat.m_LabelToIndex.clear()

        for i in range(seq_count):
            seq = <_Sequence*> sequences._mseq.GetSequence(i)
            seq.m_SMI = i
            label = seq.GetLabel()
            self._mpcflat.m_Labels.push_back(label)
            if self._mpcflat.m_LabelToIndex.find(label) != self._mpcflat.m_LabelToIndex.end():
                raise KeyError("Duplicate label: {!r}".format(label.decode('utf-8', 'replace')))
            self._mpcflat.m_LabelToIndex[label] = i

    cpdef object align(
        self,
        MultiSequence sequences,
        unsigned int seed = 0,
    ):
        cdef Alignment    msa        = Alignment.__new__(Alignment)
        cdef unsigned int seq_count  = sequences._mseq.GetSeqCount()
        cdef unsigned int pair_count = seq_count*(seq_count - 1) / 2

        if seq_count == 0:
            msa._msa.Clear()
        elif seq_count == 1:
            msa._msa.FromSequence(sequences._mseq.m_Seqs[0][0])
        else:
            self._mpcflat.Clear()
            self._alloc_pair_count(pair_count)
            self._init_seqs(sequences._mseq)
            self._mpcflat.InitPairs()
            self._mpcflat.InitDistMx()
            self._mpcflat.CalcPosteriors()
            self._mpcflat.Consistency()
            self._mpcflat.CalcGuideTree()
            self._mpcflat.CalcJoinOrder()
            self._mpcflat.ProgressiveAlign()
            self._mpcflat.Refine()
            self._mpcflat.m_MSA.ToMSA(msa._msa)

        return msa
