AM_CXXFLAGS=$(DEBUG_OPTIONS) -I$(top_srcdir)/src/c -I$(top_srcdir)/src/c/ArrayLibrary -D_LARGEFILE_SOURCE

libpercolator_la_SOURCES = \
  BaseSpline.cpp \
  BaseSpline.h \
  Caller.cpp \
  Caller.h \
  DataSet.cpp \
  DataSet.h \
  Globals.cpp \
  Globals.h \
  IntraSetRelation.cpp \
  IntraSetRelation.h \
  LogisticRegression.cpp \
  LogisticRegression.h \
  Normalizer.cpp \
  Normalizer.h \
  Option.cpp \
  Option.h \
  PercolatorCInterface.cpp \
  PercolatorCInterface.h \
  PosteriorEstimator.cpp \
  PosteriorEstimator.h \
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
  ssl.h \
  Transform.h 

noinst_LTLIBRARIES = libpercolator.la

ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

