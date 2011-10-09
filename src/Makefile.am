AM_CPPFLAGS = -I@top_srcdir@/include

bin_PROGRAMS = classdecode

lib_LTLIBRARIES = libcppnlp.la

libcppnlp_la_CXXFLAGS = -Wall -std=c++0x
libcppnlp_la_SOURCES = classdecoder.cpp ngram.cpp

classdecode_CXXFLAGS = -Wall -std=c++0x
classdecode_LDADD = libcppnlp.la
classdecode_SOURCES = classdecoder.cpp classdecode.cpp
classdecode_LDFLAGS = -static
#LDADD = libcppnlp.la

#lib_LTLIBRARIES = libcppnlp.la
#libcppnlp_la_LDFLAGS = -version-info 1:0:0
#libcppnlp_la_SOURCES = classdecoder.cpp ngram.cpp
