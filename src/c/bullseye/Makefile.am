# Automake template for xlink-search library


AM_CXXFLAGS = -Wall \
              -W \
              $(USEDOUBLES) \
              $(USEDARWIN) \
              -D_LARGEFILE_SOURCE=1 \
              -D_FILE_OFFSET_BITS=64 \
              -D_GNU_SOURCE


libbullseye_la_SOURCES = \
	CKronik2.h \
	CKronik2.cpp \
	bullseye.cpp \
	CruxBullseyeApplication.h \
        CruxBullseyeApplication.cpp


libbullseye_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/ \
  -I$(top_srcdir)/src/external/MSToolkit/include \
  -I$(top_srcdir)/src/c/hardklor \
  -DGCC \
  -D_LARGEFILE_SOURCE \
  -DHAVE_EXPAT_CONFIG_H \
  -D_FILE_OFFSET_BITS=64 \
  -DCRUX


noinst_LTLIBRARIES = libbullseye.la

