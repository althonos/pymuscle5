from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef extern from "hmmparams.h":

    cdef cppclass HMMParams:
        bool m_Logs
        unsigned int m_LineNr
        float m_Var

        vector[float] m_Trans
        vector[vector[float]] m_Emits
        vector[string] m_Lines
        string m_Alpha

        HMMParams()
        void Clear()
        unsigned int GetAlphaSize() const

        void FromParams(const HMMParams& Params, bool AsProbs)
        void FromStrings(const vector[string]& Lines)
        void FromFile(const string& FileName)
        void FromDefaults(bool Nucleo)

        void PerturbProbs(unsigned int Seed)
        void ToSingleAffineProbs(HMMParams& Params)

        void Normalize()
        void NormalizeStart()
        void NormalizeMatch()
        void NormalizeShortGap()
        void NormalizeLongGap()
        void NormalizeEmit()

        void AssertProbsValid() const

        void ToFile(const string& FileName) const
        void ToPairHMM() const

        void GetProbs(HMMParams& Probs) const
        void GetScore(HMMParams& Scores) const

        
