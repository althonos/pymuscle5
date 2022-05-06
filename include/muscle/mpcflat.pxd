from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair

from muscle.multisequence cimport MultiSequence
from muscle.mysparsemx cimport MySparseMx
from muscle.tree cimport Tree
from muscle.treeperm cimport TREEPERM
from muscle.upgma5 cimport UPGMA5


cdef extern from "mpcflat.h":

    cdef cppclass MPCFlat:
        MultiSequence* m_InputSeqs
        MultiSequence* m_MSA

        unsigned int m_ConsistencyIterCount
        unsigned int m_RefineIterCount
        TREEPERM m_TreePerm
        vector[string] m_Labels
        map[string, unsigned int] m_LabelToIndex

        UPGMA5 m_Upgma5
        Tree m_GuideTree
        vector[MultiSequence*] m_ProgMSAs

        vector[vector[float]] m_DistMx
        vector[pair[unsigned int, unsigned int]] m_Pairs
        map[pair[unsigned int, unsigned int], unsigned int] m_PairToIndex
        vector[unsigned int] m_JoinIndexes1
        vector[unsigned int] m_JoinIndexes2

        vector[MySparseMx*] m_SparsePosts1
        vector[MySparseMx*] m_SparsePosts2
        vector[MySparseMx*] m_ptrSparsePosts
        vector[MySparseMx*] m_ptrUpdatedSparsePosts

        void Clear()
        void Run(MultiSequence* InputSeqs)
        unsigned int GetSeqCount() const
        void Run_Super4(MultiSequence* InputSeqs)
