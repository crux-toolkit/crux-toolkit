#Set your paths here.
ZLIB_PATH = ./src/zLib-1.2.5
MZPARSER_PATH = ./src/mzParser
EXPAT_PATH = ./src/expat-2.0.1
SQLITE_PATH = ./src/sqlite-3.7.7.1
MST_PATH = ./src/MSToolkit
INCLUDE_PATH=./include

libmstoolkit_la_SOURCES = \
  $(INCLUDE_PATH)/ascii.h \
  $(INCLUDE_PATH)/asciitab.h \
  $(INCLUDE_PATH)/crc32.h \
  $(INCLUDE_PATH)/deflate.h \
  $(INCLUDE_PATH)/expat_config.h \
  $(INCLUDE_PATH)/expat_external.h \
  $(INCLUDE_PATH)/expat.h \
  $(INCLUDE_PATH)/gzguts.h \
  $(INCLUDE_PATH)/iasciitab.h \
  $(INCLUDE_PATH)/inffast.h \
  $(INCLUDE_PATH)/inffixed.h \
  $(INCLUDE_PATH)/inflate.h \
  $(INCLUDE_PATH)/inftrees.h \
  $(INCLUDE_PATH)/internal.h \
  $(INCLUDE_PATH)/latin1tab.h \
  $(INCLUDE_PATH)/MSObject.h \
  $(INCLUDE_PATH)/MSReader.h \
  $(INCLUDE_PATH)/MSToolkitTypes.h \
  $(INCLUDE_PATH)/mzParser.h \
  $(INCLUDE_PATH)/nametab.h \
  $(INCLUDE_PATH)/Spectrum.h \
  $(INCLUDE_PATH)/sqlite3.h \
  $(INCLUDE_PATH)/trees.h \
  $(INCLUDE_PATH)/utf8tab.h \
  $(INCLUDE_PATH)/winconfig.h \
  $(INCLUDE_PATH)/xmlrole.h \
  $(INCLUDE_PATH)/xmltok.h \
  $(INCLUDE_PATH)/xmltok_impl.h \
  $(INCLUDE_PATH)/zconf.h \
  $(INCLUDE_PATH)/zlib.h \
  $(INCLUDE_PATH)/zutil.h \
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
  $(EXPAT_PATH)/xmltok_impl.c \
  $(EXPAT_PATH)/xmltok_ns.c \
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
