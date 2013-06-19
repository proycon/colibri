# cython: c_string_type=unicode, c_string_encoding=utf8
from libcpp.string cimport string
from libcpp.set cimport set as cppset
from libcpp cimport bool
from unordered_map cimport unordered_map
from libc.stdint cimport *

cdef extern from "ngram.h":
    cdef cppclass EncAnyGram:
        string decode(ClassDecoder&)
        int n()

    cdef cppclass EncNGram(EncAnyGram):
        pass

cdef extern from "patternmodel.h":
    cdef cppclass CorpusReference:
        uint8_t token
        uint32_t sentence

    cdef cppclass NGramData:
        cppset[CorpusReference] refs
        int count()

    cdef cppclass IndexedPatternModel:
        IndexedPatternModel(string, bool) except +
        bool exists(EncAnyGram*)
        unordered_map[EncNGram,NGramData] ngrams
        int types()
        int tokens()
        int occurrencecount(EncAnyGram*)

cdef extern from "classdecoder.h":
    cdef cppclass ClassDecoder:
        ClassDecoder(string) except +
        int size()

cdef extern from "classencoder.h":
    cdef cppclass ClassEncoder:
        ClassEncoder(string) except +
        int size()
        EncAnyGram* input2anygram(string , bool allowunknown, bool autoaddunknown)
