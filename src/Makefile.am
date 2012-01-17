AM_CPPFLAGS = -I@top_srcdir@/include

bin_PROGRAMS = classdecode classencode patternfinder graphmodel #aligner

lib_LTLIBRARIES = libcolibri.la

libcolibri_la_CXXFLAGS = -Wall -std=c++0x
libcolibri_la_SOURCES = classdecoder.cpp classencoder.cpp ngram.cpp algorithms.cpp patternmodel.cpp alignmodel.cpp common.cpp

classdecode_CXXFLAGS = -Wall -std=c++0x
classdecode_LDADD = libcolibri.la
classdecode_SOURCES = classdecode.cpp
classdecode_LDFLAGS = -static

classencode_CXXFLAGS = -Wall -std=c++0x
classencode_LDADD = libcolibri.la
classencode_SOURCES = classencode.cpp
classencode_LDFLAGS = -static


graphmodel_CXXFLAGS = -Wall -std=c++0x
graphmodel_LDADD = libcolibri.la
graphmodel_SOURCES = graphmodel.cpp
graphmodel_LDFLAGS = -static

patternfinder_CXXFLAGS = -Wall -std=c++0x
patternfinder_LDADD = libcolibri.la
patternfinder_SOURCES = patternfinder.cpp
patternfinder_LDFLAGS = -static

aligner_CXXFLAGS = -Wall -std=c++0x
aligner_LDADD = libcolibri.la
aligner_SOURCES = aligner.cpp
aligner_LDFLAGS = -static

testskipgrams_CXXFLAGS = -Wall -std=c++0x
testskipgrams_LDADD = libcolibri.la
testskipgrams_SOURCES = testskipgrams.cpp
testskipgrams_LDFLAGS = -static

#lib_LTLIBRARIES = libcppnlp.la
#libcppnlp_la_LDFLAGS = -version-info 1:0:0
#libcppnlp_la_SOURCES = classdecoder.cpp ngram.cpp
