from libcpp.vector cimport vector

from muscle cimport byte
from muscle.pairhmm cimport DEFAULT_CONSISTENCY_ITERS


cdef extern from "mysparsemx.h":

    cdef cppclass MySparseMx:
        unsigned int m_LX
        unsigned int m_LY
        unsigned int m_VecSize
        unsigned int m_MaxVecSize
        byte* m_ValueVec
        unsigned int m_MaxLX
        unsigned int* m_Offsets

        const byte* m_X
        const byte* m_Y

        MySparseMx()

        void Clear()
        float GetProb_Offset(unsigned int Offset) const
        unsigned int GetCol_Offset(unsigned int Offset) const
        void SetProb_Offset(unsigned int Offset, float P)
        void SetCol_Offset(unsigned int Offset, unsigned int Col) const

        void AllocLX(unsigned int LX)
        void AllocVec(unsigned int Size)
        void FromPost(const float* Post, unsigned int LX, unsigned int LY)
        void UpdateFromPost(const MySparseMx& OldMX, const float* Post, unsigned int SeqCount)
        void GetColToRowLoHi(vector[unsigned int]& ColToRowLo, vector[unsigned int]& ColToRowHi) const
        void ToPost(float* Post) const
        float GetProb(unsigned int i, unsigned int j) const
        unsigned int GetOffset(unsigned int i) const
        unsigned int GetSize(unsigned int i) const
        float GetMaxProbRow(unsigned int i) const

        void LogMe() const
        void LogStats(const char** Msg) const

        const unsigned int GetLX() const
        const unsigned int GetLY() const
