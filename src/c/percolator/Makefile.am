distdir = $(PACKAGE)_$(VERSION)_src

AM_CPPFLAGS=$(DEBUG_OPTIONS) -I$(srcdir)/ArrayLibrary

libqvality_la_SOURCES = \
  BaseSpline.cpp \
  BaseSpline.h \
  Globals.cpp \
  Globals.h \
  LogisticRegression.cpp \
  LogisticRegression.h \
  Option.cpp \
  Option.h \
  PosteriorEstimator.cpp \
  PosteriorEstimator.h \
  Transform.h 


libArrayLibrary_la_SOURCES = \
  ArrayLibrary/ArrayTypes.h \
  ArrayLibrary/Array.h \
  ArrayLibrary/Arithmetic.h \
  ArrayLibrary/ArrayLibrary.h \
  ArrayLibrary/MatrixFunctions.cpp \
  ArrayLibrary/MatrixFunctions.h \
  ArrayLibrary/Numerical.cpp \
  ArrayLibrary/Numerical.h \
  ArrayLibrary/OrderedArray.h \
  ArrayLibrary/PackedMatrix.cpp \
  ArrayLibrary/PackedMatrix.h \
  ArrayLibrary/PackedSetFunctions.cpp \
  ArrayLibrary/PackedSetFunctions.h \
  ArrayLibrary/PackedVec.cpp \
  ArrayLibrary/PackedVec.h \
  ArrayLibrary/Random.cpp \
  ArrayLibrary/Random.h \
  ArrayLibrary/VectorFunctions.cpp \
  ArrayLibrary/VectorFunctions.h


libpercolator_la_SOURCES = \
  Caller.cpp \
  Caller.h \
  DataSet.cpp \
  DataSet.h \
  IntraSetRelation.cpp \
  IntraSetRelation.h \
  Normalizer.cpp \
  Normalizer.h \
  Option.cpp \
  Option.h \
  PercolatorCInterface.cpp \
  PercolatorCInterface.h \
  ResultHolder.cpp \
  ResultHolder.h \
  SanityCheck.h \
  SanityCheck.cpp \
  Scores.cpp \
  Scores.h \
  SetHandler.cpp \
  SetHandler.h \
  SqtSanityCheck.h \
  SqtSanityCheck.cpp \
  StdvNormalizer.cpp \
  StdvNormalizer.h \
  UniNormalizer.cpp \
  UniNormalizer.h \
  ssl.cpp \
  ssl.h

noinst_LTLIBRARIES = libqvality.la libArrayLibrary.la libpercolator.la

ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

EXTRA_DIST = \
  ArrayLibrary/Array.cpp \
  ArrayLibrary/Arithmetic.cpp \
  ArrayLibrary/OrderedArray.cpp

