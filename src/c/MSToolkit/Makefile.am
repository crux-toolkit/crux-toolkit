libmstoolkit_la_SOURCES = \
  MSConvertFile.cpp \
  MSObject.cpp \
  MSObject.h \
  MSReader.cpp \
  MSReader.h \
  MSSingleScan.cpp \
  MSTDemo.cpp \
  MSToolkitTypes.h \
  RAMP/ramp_base64.cpp \
  RAMP/ramp_base64.h \
  RAMP/ramp.cpp \
  RAMP/ramp.h \
  Spectrum.cpp \
  Spectrum.h \
  zLib/adler32.c \
  zLib/compress.c \
  zLib/crc32.c \
  zLib/crc32.h \
  zLib/deflate.c \
  zLib/deflate.h \
  zLib/gzio.c \
  zLib/infback.c \
  zLib/inffast.c \
  zLib/inffast.h \
  zLib/inffixed.h \
  zLib/inflate.c \
  zLib/inflate.h \
  zLib/inftrees.c \
  zLib/inftrees.h \
  zLib/trees.c \
  zLib/trees.h \
  zLib/uncompr.c \
  zLib/zconf.h \
  zLib/zlib.h \
  zLib/zutil.c \
  zLib/zutil.h

libmstoolkit_la_CPPFLAGS = $(AM_CPPFLAGS) \
  -I$(top_srcdir)/src/c/MSToolkit/zLib \
  -I$(top_srcdir)/src/c/MSToolkit/RAMP \
  -D_LARGEFILE_SOURCE \
  -D_FILE_OFFSET_BITS=64

noinst_LTLIBRARIES = libmstoolkit.la
