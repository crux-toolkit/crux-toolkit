#Set your paths here.
ZLIB_PATH = ./src/zLib-1.2.5
MZPARSER_PATH = ./src/mzParser
EXPAT_PATH = ./src/expat-2.0.1
SQLITE_PATH = ./src/sqlite-3.7.7.1
MST_PATH = ./src/MSToolkit

libmstoolkit_la_SOURCES = \
  $(ZLIB_PATH)/adler32.c \
  $(ZLIB_PATH)/compress.c \
  $(ZLIB_PATH)/crc32.c \
  $(ZLIB_PATH)/deflate.c \
  $(ZLIB_PATH)/inffast.c \
  $(ZLIB_PATH)/inflate.c \
  $(ZLIB_PATH)/infback.c \
  $(ZLIB_PATH)/inftrees.c \
  $(ZLIB_PATH)/trees.c \
  $(ZLIB_PATH)/uncompr.c \
  $(ZLIB_PATH)/zutil.c \
  $(MZPARSER_PATH)/mzp_base64.cpp \
  $(MZPARSER_PATH)/BasicSpectrum.cpp \
  $(MZPARSER_PATH)/mzParser.cpp \
  $(MZPARSER_PATH)/RAMPface.cpp \
  $(MZPARSER_PATH)/saxhandler.cpp \
  $(MZPARSER_PATH)/saxmzmlhandler.cpp \
  $(MZPARSER_PATH)/saxmzxmlhandler.cpp \
  $(MZPARSER_PATH)/Czran.cpp \
  $(EXPAT_PATH)/xmlparse.c \
  $(EXPAT_PATH)/xmlrole.c \
  $(EXPAT_PATH)/xmltok.c \
  $(SQLITE_PATH)/sqlite3.c \
  $(MST_PATH)/Spectrum.cpp \
  $(MST_PATH)/MSReader.cpp \
  $(MST_PATH)/MSObject.cpp

libmstoolkit_la_CXXFLAGS = $(AM_CPPFLAGS) \
  -I$(top_srcdir)/src/external/MSToolkit/include \
  -D_LARGEFILE_SOURCE \
  -DGCC \
  -DHAVE_EXPAT_CONFIG_H \
  -D_FILE_OFFSET_BITS=64

libmstoolkit_la_CFLAGS = $(AM_CFLAGS) \
  -I$(top_srcdir)/src/external/MSToolkit/include \
  -D_LARGEFILE_SOURCE \
  -DGCC \
  -DHAVE_EXPAT_CONFIG_H \
  -D_FILE_OFFSET_BITS=64


noinst_LTLIBRARIES = libmstoolkit.la
