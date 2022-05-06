from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map

from muscle cimport LINKAGE
from muscle.tree cimport Tree

cdef extern from "upgma5.h":

    cdef cppclass UPGMA5:
        unsigned int m_LeafCount
        unsigned int m_TriangleSize
        unsigned int m_InternalNodeCount
        unsigned int m_InternalNodeIndex

        float* m_Dist
        float* m_MinDist
        unsigned int* m_NearestNeighbor
        unsigned int* m_NodeIndex

        unsigned int* m_Left
        unsigned int* m_Right
        float* m_Height
        float* m_LeftLength
        float* m_RightLength

        vector[string] m_Labels
        vector[vector[float]] m_DistMx
        map[string, unsigned int] m_LabelToIndex

        void Clear()
        void Init(const vector[string]& Labels, const vector[vector[float]]& DistMx)
        void Run(LINKAGE Linkage, Tree& tree)
        void ReadDistMx(const string& FileName)
        void ScaleDistMx()
        void FixEADistMx()
        void LogMe() const
        void AddLabel(const string& Label)
        unsigned int GetLabelIndex(const string& Label) const

        unsigned int TriangleSubscript(unsigned int uIndex1, unsigned int uIndex2) const
