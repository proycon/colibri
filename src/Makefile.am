AM_CPPFLAGS = -I@top_srcdir@/include

bin_PROGRAMS = classdecode modeldecode patternfinder graphmodel

lib_LTLIBRARIES = libcppnlp.la

libcppnlp_la_CXXFLAGS = -Wall -std=c++0x
libcppnlp_la_SOURCES = classdecoder.cpp ngram.cpp algorithms.cpp

classdecode_CXXFLAGS = -Wall -std=c++0x
classdecode_LDADD = libcppnlp.la
classdecode_SOURCES = classdecode.cpp
classdecode_LDFLAGS = -static

modeldecode_CXXFLAGS = -Wall -std=c++0x
modeldecode_LDADD = libcppnlp.la
modeldecode_SOURCES = modeldecode.cpp
modeldecode_LDFLAGS = -static

graphmodel_CXXFLAGS = -Wall -std=c++0x
graphmodel_LDADD = libcppnlp.la
graphmodel_SOURCES = graphmodel.cpp
graphmodel_LDFLAGS = -static

patternfinder_CXXFLAGS = -Wall -std=c++0x
patternfinder_LDADD = libcppnlp.la
patternfinder_SOURCES = patternfinder.cpp
patternfinder_LDFLAGS = -static
#LDADD = libcppnlp.la

testskipgrams_CXXFLAGS = -Wall -std=c++0x
testskipgrams_LDADD = libcppnlp.la
testskipgrams_SOURCES = testskipgrams.cpp
testskipgrams_LDFLAGS = -static

#lib_LTLIBRARIES = libcppnlp.la
#libcppnlp_la_LDFLAGS = -version-info 1:0:0
#libcppnlp_la_SOURCES = classdecoder.cpp ngram.cpp
