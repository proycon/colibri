from libcpp.string cimport string
from libcpp.set cimport set as cppset
from libcpp.vector cimport vector
from libcpp cimport bool
from unordered_map cimport unordered_map
from libc.stdint cimport *

#cdef extern from *:
#    ctypedef char* const_EncAnyGram "const EncAnyGram"

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

    cdef cppclass AnyGramData:
        pass

    cdef cppclass NGramData(AnyGramData):
        cppset[CorpusReference] refs
        int count()

    cdef cppclass IndexedPatternModel:
        IndexedPatternModel(string, bool, bool) except +
        bool exists(EncAnyGram*)
        unordered_map[EncNGram,NGramData] ngrams
        int types()
        int tokens()
        int occurrencecount(EncAnyGram*) except +
        AnyGramData * getdata(EncAnyGram*)
        vector[const EncAnyGram*] get_reverse_index(int i)
        #unordered_map[int, vector[EncAnyGram*]] reverse_index

cdef extern from "classdecoder.h":
    cdef cppclass ClassDecoder:
        ClassDecoder(string) except +
        int size()

cdef extern from "classencoder.h":
    cdef cppclass ClassEncoder:
        ClassEncoder(string) except +
        int size()
        EncAnyGram* input2anygram(string , bool allowunknown, bool autoaddunknown)
