# Automake template for xlink-search library
AM_CXXFLAGS = -Wall -W $(USEDOUBLES) $(USEDARWIN) $(SEARCH_ENABLED) \
            -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE


libqranker_barista_la_SOURCES = \
	Barista.cpp \
	Barista.h \
	DataSet.h \
	ProtScores.cpp \
	QRanker.h \
	NeuralNet.cpp \
        ProtScores.h \
        SpecFeatures.cpp \
	BipartiteGraph.cpp \
	NeuralNet.h\
	PSMScores.cpp\
	SpecFeatures.h\
	BipartiteGraph.h\
	PepScores.cpp\
	PSMScores.h\
	SQTParser.cpp\
	TabDelimParser.cpp\
	TabDelimParser.h\
	DataSet.cpp\
	PepScores.h\
	QRanker.cpp\
	SQTParser.h

libqranker_barista_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/

noinst_LTLIBRARIES = libqranker_barista.la

