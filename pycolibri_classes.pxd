from libcpp.string cimport string

cdef extern from "patternmodel.h":
    cdef cppclass IndexedPatternModel:
        IndexedPatternModel(string, bool) except +

cdef extern from "classdecoder.h":
    cdef cppclass ClassDecoder:
        ClassDecoder(string) except +

cdef extern from "classencoder.h":
    cdef cppclass ClassEncoder:
        ClassEncoder(string) except +


