from libc.stdio cimport FILE
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector

from muscle.sequence cimport Sequence
from muscle.alpha cimport ALPHA


cdef extern from "msa.h" nogil:

    cdef cppclass MSA:
        unsigned m_uSeqCount
        unsigned m_uColCount
        unsigned m_uCacheSeqLength
        unsigned m_uCacheSeqCount

        char** m_szSeqs
        char** m_szNames

        unsigned* m_IdToSeqIndex
        unsigned* m_SeqIndexToId

        MSA()
        void FromStrings(const vector[string]& Strings)
        void FromStrings2(const vector[string]& Labels, vector[string]& Seqs)
        # void FromFile(TextFile& File)
        void FromFASTAFile(const string& FileName)
        # void FromFASTAFile(TextFile& File)
        void FromFASTAFile_PreserveCase(const string& FileName)
        void FromSeq(const Sequence& s)
        void FromSequence(const Sequence& s)

        void GetLabelToSeqIndex(vector[string]& Labels, map[string, unsigned int]& LabelToIndex) const

        # void ToFile(TextFile& File) const
        # void ToFASTAFile(TextFile& File) const
        void ToFASTAFile(FILE* f) const
        void ToFASTAFile(const string& FileName) const

        void GetSeqLabel(unsigned int SeqIndex, string& Label) const

        void Copy(const MSA &msa)

        ALPHA GuessAlpha() const
        void FixAlpha()

        void LogMe() const
        void Clear()

        unsigned GetSeqCount() const
        unsigned GetColCount() const
