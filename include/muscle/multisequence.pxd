# distutils: sources = multisequence.cpp


from libc.stdio cimport FILE
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.set cimport set

from muscle cimport byte
from muscle.filebuffer cimport FileBuffer
from muscle.sequence cimport Sequence
from muscle.msa cimport MSA

cdef extern from "multisequence.h" nogil:

    cdef cppclass MultiSequence:
        vector[const Sequence*] m_Seqs
        vector[bool] m_Owners

        MultiSequence()
        MultiSequence(FileBuffer& infile)
        MultiSequence(const string& filename)

        void Clear()
        void FromStrings(const vector[string]& Labels, const vector[string]& Seqs)
        void Copy(const MultiSequence &rhs)
        void LoadMFA(const string& filename, bool stripGaps)
        void LoadMFA(FileBuffer& infile, bool stripGaps)
        void FromFASTA(const string& filename, bool stripGaps = false)

        unsigned int GetSeqIndex(const string &Label, bool FailOnError = true) const
        void ToMSA(MSA& msa) const

        void AddSequence(const Sequence* sequence, bool Owner)

        void WriteMFA(FILE* f) const
        void WriteMFA(const string& FileName) const

        const Sequence* GetSequence(int i) const
        int GetNumSequences() const
        unsigned int GetSeqCount() const

        double GetMeanSeqLength() const
        unsigned int GetMaxSeqLength() const

        unsigned int GetChar(unsigned int SeqIndex, unsigned int ZeroBasedPos) const

        MultiSequence* Project(const set[int]& indices)

        bool IsAligned() const
        unsigned int GetColCount() const

        const string& GetLabelStr(unsigned int SeqIndex) const
        const char* GetLabel(unsigned int SeqIndex) const
        const byte* GetBytePtr(unsigned int SeqIndex) const

        bool GuessIsNucleo() const
        void LogGSIs(const char* Msg = NULL) const
        void AssertGSIs() const
        void GetLengthOrder(vector[unsigned int]& SeqIndexes) const
        unsigned int GetSeqLength(unsigned int SeqIndex) const
        unsigned int GetGSI(unsigned int SeqIndex) const
        void LogMe() const
