AM_CPPFLAGS = -I@top_srcdir@/include

bin_PROGRAMS = classdecode classencode patternfinder grapher aligner querier giza2html decoder trainclassifiers test

lib_LTLIBRARIES = libcolibri.la

libcolibri_la_CXXFLAGS = -Wall -std=c++11
libcolibri_la_SOURCES = classdecoder.cpp classencoder.cpp ngram.cpp algorithms.cpp patternmodel.cpp alignmodel.cpp gizamodel.cpp lm.cpp common.cpp classifiers.cpp

classdecode_CXXFLAGS = -Wall -std=c++11
classdecode_LDADD = libcolibri.la
classdecode_SOURCES = classdecode.cpp
classdecode_LDFLAGS = -static

classencode_CXXFLAGS = -Wall -std=c++11
classencode_LDADD = libcolibri.la
classencode_SOURCES = classencode.cpp
classencode_LDFLAGS = -static


grapher_CXXFLAGS = -Wall -std=c++11
grapher_LDADD = libcolibri.la
grapher_SOURCES = grapher.cpp
grapher_LDFLAGS = -static

patternfinder_CXXFLAGS = -Wall -std=c++11
patternfinder_LDADD = libcolibri.la
patternfinder_SOURCES = patternfinder.cpp
patternfinder_LDFLAGS = -static

aligner_CXXFLAGS = -Wall -std=c++11
aligner_LDADD = libcolibri.la
aligner_SOURCES = aligner.cpp
aligner_LDFLAGS = -static


decoder_CXXFLAGS = -Wall -std=c++11
decoder_LDADD = libcolibri.la
decoder_SOURCES = decoder.cpp
decoder_LDFLAGS = -static


trainclassifiers_CXXFLAGS = -Wall -std=c++11
trainclassifiers_LDADD = libcolibri.la
trainclassifiers_SOURCES = trainclassifiers.cpp
trainclassifiers_LDFLAGS = -static

querier_CXXFLAGS = -Wall -std=c++11
querier_LDADD = libcolibri.la
querier_SOURCES = querier.cpp
querier_LDFLAGS = -static

giza2html_CXXFLAGS = -Wall -std=c++11
giza2html_LDADD = libcolibri.la
giza2html_SOURCES = giza2html.cpp
giza2html_LDFLAGS = -static


test_CXXFLAGS = -Wall -std=c++11
test_LDADD = libcolibri.la
test_SOURCES = test.cpp
test_LDFLAGS = -static

#testskipgrams_CXXFLAGS = -Wall -std=c++11
#testskipgrams_LDADD = libcolibri.la
#testskipgrams_SOURCES = testskipgrams.cpp
#testskipgrams_LDFLAGS = -static

#lib_LTLIBRARIES = libcppnlp.la
#libcppnlp_la_LDFLAGS = -version-info 1:0:0
#libcppnlp_la_SOURCES = classdecoder.cpp ngram.cpp
