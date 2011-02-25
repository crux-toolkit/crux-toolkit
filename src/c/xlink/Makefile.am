# Automake template for xlink-search library
AM_CXXFLAGS = -Wall -W -Werror  $(USEDOUBLES) $(USEDARWIN) $(SEARCH_ENABLED) \
            -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE


libxlink_la_SOURCES = \
  xhhc.cpp \
  xhhc.h \
  xhhc_ion_series.cpp \
  xhhc_ion_series.h \
  xhhc_make_histogram.cpp \
  xhhc_predict_peptide_ions.cpp \
  xhhc_score_peptide_spectrum.cpp \
  xhhc_scorer.cpp \
  xhhc_score_peptide_spectrum.cpp \
  xhhc_scorer.cpp \
  xhhc_scorer.h \
  xhhc_search.cpp \
  SearchForXLinks.cpp \
  SearchForXLinks.h \
  xlink_assign_ions.cpp \
  xlink_compute_qvalues.cpp \
  XLinkSite.cpp \
  XLinkSite.h \
  XLinkBondMap.cpp \
  XLinkBondMap.h




libxlink_la_CXXFLAGS = $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/ \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/MSToolkit/zLib \
  -I$(top_srcdir)/src/c/MSToolkit/RAMP

noinst_LTLIBRARIES = libxlink.la

