from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.pair cimport pair

from muscle cimport byte
from muscle.sequence cimport Sequence
from muscle.multisequence cimport MultiSequence
from muscle.mysparsemx cimport MySparseMx
from muscle.tree cimport Tree
from muscle.treeperm cimport TREEPERM
from muscle.upgma5 cimport UPGMA5


cdef extern from "mpcflat.h" nogil:

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

        vector[MySparseMx*]  m_SparsePosts1
        vector[MySparseMx*]  m_SparsePosts2
        vector[MySparseMx*]* m_ptrSparsePosts
        vector[MySparseMx*]* m_ptrUpdatedSparsePosts

        void Clear()
        void Run(MultiSequence* InputSeqs)
        unsigned int GetSeqCount() const
        void Run_Super4(MultiSequence* InputSeqs)

        void AllocPairCount(unsigned int SeqCount)
        void FreeProgMSAs()
        void FreeSparsePosts()
        unsigned int GetL(unsigned int SeqIndex) const
        void InitSeqs(MultiSequence *InputSeqs)
        void InitPairs();
        void InitDistMx()
        void CalcPosteriors()
        void CalcPosterior(unsigned int PairIndex)
        void Consistency()
        void ConsIter(unsigned int Iter)
        void ConsPair(unsigned int PairIndex)
        void CalcGuideTree()
        void CalcGuideTree_RandomChain()
        void CalcJoinOrder()
        void ProgressiveAlign()
        void Refine()
        void RefineIter()
        void ProgAln(unsigned int JoinIndex)
        const pair[unsigned int, unsigned int] &GetPair(unsigned int PairIndex) const
        const char* GetLabel(unsigned int SeqIndex) const
        const byte* GetBytePtr(unsigned int SeqIndex) const
        unsigned int GetPairIndex(unsigned int SMI1, unsigned int SMI2) const
        MySparseMx& GetSparsePost(unsigned int PairIndex)
        MySparseMx& GetUpdatedSparsePost(unsigned int PairIndex)
        MultiSequence* AlignAlns(const MultiSequence &MSA1, const MultiSequence &MSA2)
        void BuildPost(const MultiSequence &MSA1, const MultiSequence &MSA2, float *Post)
        unsigned int GetSeqLength(unsigned int SeqIndex) const
        const Sequence* GetSequence(unsigned int SeqIndex) const
