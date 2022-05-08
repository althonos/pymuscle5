from libcpp.string cimport string
from libcpp.vector cimport vector


cdef extern from "tree.h" nogil:

    cdef enum NEWICK_TOKEN_TYPE:
        NTT_Unknown
        NTT_Lparen
        NTT_Rparen
        NTT_Colon
        NTT_Comma
        NTT_Semicolon
        NTT_String
        NTT_SingleQuotedString
        NTT_DoubleQuotedString
        NTT_Comment

    cdef cppclass Tree:
        Tree()
        void Clear()

        void CreateRooted()
        void CreateUnrooted(double dEdgeLength)

        void FromFile(const string& FileName)
        # void FromFile(TextFile& File)
        # void FromClust(Clust& C)

        void Copy(const Tree& tree)

        void Create(
            unsigned uLeafCount,
            unsigned uRoot,
            const unsigned Left[],
            const unsigned Right[],
            const float LeftLength[],
            const float RightLength[],
            const unsigned LeafIds[],
            char* LeafNames[]
        )
        void FromVectors(
            const vector[string]& Labels,
            const vector[unsigned int]& Parents,
            const vector[float]& Lengths
        )
        void ToVectors(
            vector[string]& Labels,
            vector[float]& Lengths
        ) const
