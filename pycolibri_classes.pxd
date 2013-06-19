from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "ngram.h":
    cdef cppclass EncAnyGram:
        string decode(ClassDecoder&)

cdef extern from "patternmodel.h":
    cdef cppclass IndexedPatternModel:
        IndexedPatternModel(string, bool) except +
        bool exists(EncAnyGram *)

cdef extern from "classdecoder.h":
    cdef cppclass ClassDecoder:
        ClassDecoder(string) except +
        int size()

cdef extern from "classencoder.h":
    cdef cppclass ClassEncoder:
        ClassEncoder(string) except +
        int size()

