from libc.stdio cimport FILE
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from muscle cimport byte
from muscle.filebuffer cimport FileBuffer

cdef extern from "sequence.h" nogil:

    cdef cppclass Sequence:
        unsigned int m_GSI
        unsigned int m_SMI

        string       m_Label
        vector[char] m_CharVec

        bool FromFileBuffer(FileBuffer& infile, bool stripGaps = false)
        void Create(const vector[char]* m_CharVec, string m_Label, unsigned int GSI, unsigned int SMI)

        void InitData()
        void AppendChar()

        const string& GetLabel() const
        const char* GetLabelCStr() const

        char* GetCharPtr1()
        const char* GetCharPtr1() const
        const byte* GetBytePtr() const

        char GetPosition(int i) const
        char GetChar(unsigned int ZeroBasedPos) const

        void SetGSI(unsigned int GSI)
        void OverwriteGSI(unsigned int GSI)
        void OverwriteLabel(const string& Label)

        unsigned int GetSMI() const
        unsigned int GetGSI() const
        unsigned int GetLength() const

        void WriteMFA(FILE* f) const

        Sequence* Clone() const
        Sequence* AddGaps(const vector[char]* alignment, char id) const
        Sequence* AddGapsPath(const string& Path, char id) const
        Sequence* DeleteGaps() const

        void GetPosToCol_OneBased(vector[unsigned int]& PosToCol) const
        void GetPosToCol(vector[unsigned int]& PosToCol) const
        void GetColToPos(vector[unsigned int]& ColToPos) const
        void LogMe() const

        # static Sequence* _NewSequence()
        @staticmethod
        Sequence* NewSequence()
        # static void _DeleteSequence(const Sequence* Seq)
        @staticmethod
        void DeleteSequence(const Sequence* Seq)
