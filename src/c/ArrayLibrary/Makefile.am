libArrayLibrary_la_SOURCES = \
  ArrayTypes.h \
  Array.h \
  Arithmetic.h \
  ArrayLibrary.h \
  MatrixFunctions.cpp \
  MatrixFunctions.h \
  Numerical.cpp \
  Numerical.h \
  OrderedArray.h \
  PackedMatrix.cpp \
  PackedMatrix.h \
  PackedSetFunctions.cpp \
  PackedSetFunctions.h \
  PackedVec.cpp \
  PackedVec.h \
  Random.cpp \
  Random.h \
  VectorFunctions.cpp \
  VectorFunctions.h

noinst_LTLIBRARIES = libArrayLibrary.la

# These files are template implementations that 
# are included like header files rather then being compiled,
# so they need to be expilicity included in the distribution.
EXTRA_DIST = \
  Arithmetic.cpp \
  Array.cpp \
  OrderedArray.cpp \
  PackedArray.cpp

ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

