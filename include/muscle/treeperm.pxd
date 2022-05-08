from libcpp.string cimport string
from libcpp.vector cimport vector

from muscle cimport LINKAGE
from muscle.tree cimport Tree

cdef extern from "treeperm.h" nogil:

    cdef enum TREEPERM:
        TP_None = 0
        TP_ABC = 1
        TP_ACB = 2
        TP_BCA = 3
        TP_All = 4

    TREEPERM StrToTREEPERM(const string& s)
    const char* TREEPERMToStr(TREEPERM TP)


cdef extern from "muscle.h" nogil:

    void PermTree(Tree &InputTree, TREEPERM TP)
    void PermuteTree(
        const Tree &InputTree,
        Tree &TreeABC,
        Tree &TreeACB,
        Tree &TreeBCA,
        vector[string] &LabelsA,
        vector[string] &LabelsB,
        vector[string] &LabelsC
    )
