from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.vector cimport vector

from muscle.tree cimport Tree
from muscle.multisequence cimport MultiSequence


cdef extern from "pprog.h":

    cdef cppclass PProg:
        unsigned int m_InputMSACount
        unsigned int m_JoinCount
        unsigned int m_NodeCount
        unsigned int m_TargetPairCount
        unsigned int m_MaxCoarseSeqs

        map[string, unsigned int] m_MSALabelToIndex
        vector[string] m_MSALabels
        vector[const MultiSequence*] m_MSAs

        vector[unsigned int]  m_Pending
        vector[vector[float]] m_ScoreMx
        vector[vector[float]] m_PathMx

        unsigned int m_JoinIndex
        vector[unsigned int] m_JoinMSAIndexes1
        vector[unsigned int] m_JoinMSAIndexes2


    cdef void MakeGuideTreeFromJoinOrder(
        const vector[unsigned int] &Indexes1,
        const vector[unsigned int] &Indexes2,
        const map[string, uint] &LabelToIndex,
        Tree &GuideTree
    )
    cdef void GetGuideTreeJoinOrder(
        const Tree &GuideTree,
        const map[string, unsigned int] &LabelToIndex,
        vector[unsigned int] &Indexes1,
        vector[unsigned int] &Indexes2
    )
    cdef void ValidateJoinOrder(
        const vector[unsigned int] &Indexes1,
        const vector[unsigned int] &Indexes2
    )
