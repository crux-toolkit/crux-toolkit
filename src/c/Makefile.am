distdir = $(PACKAGE)_$(VERSION)

SUBDIRS = ArrayLibrary MSToolkit percolator q-ranker xlink .

# Crux can use either single or double precision floating poinnt.
if USEDOUBLES
  USEDOUBLES = -DUSE_DOUBLES
else
  USEDOUBLES=
endif

# Compiling under OS X requires a preprocessor flag.
if DARWIN
  USEDARWIN = -DDARWIN
else
  USEDARWIN = 
endif

# We disable search and exclude search related code
# in open-source distribution of Crux.
if SEARCH
  SEARCH_ENABLED = -DSEARCH_ENABLED
else
  SEARCH_ENABLED = 
endif

AM_CFLAGS = -Wall -W -Werror $(USEDOUBLES) $(USEDARWIN) $(SEARCH_ENABLED) \
            -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
AM_CXXFLAGS = -Wall -W -Werror  $(USEDOUBLES) $(USEDARWIN) $(SEARCH_ENABLED) \
            -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

LIBFILES= \
  alphabet.cpp \
  alphabet.h \
  carp.cpp \
  carp.h \
  create_index.cpp \
  create_index.h \
  crux-utils.cpp \
  crux-utils.h \
  database.cpp \
  database.h \
  DelimitedFile.cpp \
  DelimitedFile.h \
  DelimitedFileReader.cpp \
  DelimitedFileReader.h \
  MatchColumns.cpp \
  MatchColumns.h \
  MatchFileReader.cpp \
  MatchFileReader.h \
  generate_peptides_iterator.cpp \
  generate_peptides_iterator.h \
  hash.cpp \
  hash.h \
  index.cpp \
  index.h \
  ion.cpp \
  ion.h \
  ion_series.cpp \
  ion_series.h \
  linked_list.cpp \
  linked_list.h \
  mass.cpp \
  mass.h \
  match.cpp \
  match.h \
  match_collection.cpp \
  match_collection.h \
  match_search.cpp \
  modifications.cpp \
  modifications.h \
  modified_peptides_iterator.cpp \
  modified_peptides_iterator.h \
  search.h \
  objects.h \
  output-files.cpp \
  output-files.h \
  parameter.cpp \
  parameter.h \
  parse_arguments.cpp \
  parse_arguments.h \
  peak.cpp \
  peak.h \
  peptide.cpp \
  peptide_constraint.cpp \
  peptide_constraint.h \
  peptide.h \
  peptide_modifications.cpp \
  peptide_modifications.h \
  peptide_src.cpp \
  peptide_src.h \
  print-processed-spectra.cpp \
  print-processed-spectra.h \
  protein.cpp \
  protein.h \
  protein_index.cpp \
  protein_index.h \
  q-value.cpp \
  q-value.h \
  search-progress.h \
  sequest-search.cpp \
  sequest-search.h \
  scorer.cpp \
  scorer.h \
  sorter.cpp \
  sorter.h \
  spectrum.cpp \
  spectrum_collection.cpp \
  spectrum_collection.h \
  spectrum.h \
  utils.cpp \
  utils.h \
  version.h \
  xlink_search.h \
  xlink_compute_qvalues.h

noinst_LTLIBRARIES = libcrux.la
libcrux_la_CFLAGS = \
  $(AM_CFLAGS) \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/percolator/ \
  -I$(top_srcdir)/src/c/q-ranker/

libcrux_la_CXXFLAGS = \
  $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/MSToolkit/RAMP \
  -I$(top_srcdir)/src/c/percolator/ \
  -I$(top_srcdir)/src/c/q-ranker/

libcrux_la_SOURCES = $(LIBFILES)

PROGRAM_FILES = \
  crux \
  crux-get-ms2-spectrum \
  crux-generate-peptides \
  crux-predict-peptide-ions \
  xlink-predict-peptide-ions \
  xlink-score-spectrum \
  xlink-assign-ions

bin_PROGRAMS = $(PROGRAM_FILES)

crux_get_ms2_spectrum_SOURCES = get_ms2_spectrum.cpp
crux_get_ms2_spectrum_CXXFLAGS = $(AM_CFLAGS) -DMAIN -I$(top_srcdir)/src/c//MSToolkit
crux_get_ms2_spectrum_LDADD = \
  $(AM_LDFLAGS) \
  libcrux.la -L$(top_builddir)/src/c/MSToolkit -lm -lmstoolkit

crux_generate_peptides_SOURCES = generate_peptides.cpp
crux_generate_peptides_CXXFLAGS = $(AM_CFLAGS) -DMAIN
crux_generate_peptides_LDADD = $(AM_LDFLAGS) libcrux.la -lm  

crux_predict_peptide_ions_SOURCES = predict_peptide_ions.cpp
crux_predict_peptide_ions_CXXFLAGS = $(AM_CFLAGS) -DMAIN
crux_predict_peptide_ions_LDADD = \
  $(AM_LDFLAGS) \
  libcrux.la -L$(top_builddir)/src/c/MSToolkit -lm -lmstoolkit

crux_SOURCES = crux-main.cpp crux-main.h
crux_CXXFLAGS = \
  $(AM_CXXFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/percolator \
  -I$(top_srcdir)/src/c/q-ranker \
  -I$(top_srcdir)/src/c/xlink
crux_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/c/xlink \
  -lxlink \
  libcrux.la \
  -lstdc++ \
  -L$(top_builddir)/src/c/ArrayLibrary \
  -L$(top_builddir)/src/c/MSToolkit \
  -L$(top_builddir)/src/c/percolator \
  -L$(top_builddir)/src/c/q-ranker \
  -lm \
  -lpercolator \
  -lqranker \
  -lmstoolkit \
  -lArrayLibrary

xlink_predict_peptide_ions_SOURCES=./xlink/xhhc_predict_peptide_ions.cpp
xlink_predict_peptide_ions_CXXFLAGS= \
  $(AM_CFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/c/xlink

xlink_predict_peptide_ions_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/c/MSToolkit \
  -L$(top_builddir)/src/c/xlink \
  -lm \
  -lxlink \
  libcrux.la \
  -lmstoolkit

xlink_score_spectrum_SOURCES=./xlink/xhhc_score_peptide_spectrum.cpp
xlink_score_spectrum_CXXFLAGS= \
  $(AM_CFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/c/xlink

xlink_score_spectrum_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/c/MSToolkit \
  -L$(top_builddir)/src/c/xlink \
  -lm \
  -lxlink \
  libcrux.la \
  -lmstoolkit

xlink_assign_ions_SOURCES=./xlink/xlink_assign_ions.cpp
xlink_assign_ions_CXXFLAGS= \
  $(AM_CFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/c/xlink

xlink_assign_ions_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/c/MSToolkit \
  -L$(top_builddir)/src/c/xlink \
  -lm \
  -lxlink \
  libcrux.la \
  -lmstoolkit


paramdir = $(prefix)/doc
param_DATA = default.params

EXTRA_DIST = default.params

doc:	*.h *.cpp
	doxygen Doxyfile

ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

# In the distribution hook we create the version header,
# and if SEARCH_ENABLED is not defined, we remove peptide-spectrum database 
# search code  from the files in the distribution.
# This is needed to support open-source distribution.
dist-hook:
	grep '#define VERSION' ../../config.h > version.h
	if test -z "$(SEARCH_ENABLED)"; then \
    sed -e '/#ifdef SEARCH_ENABLED/,/#else/d' \
        -e '/#endif \/\/ SEARCH_ENABLED/d' \
        < match_search.cpp > $(distdir)/match_search.cpp; \
    sed -e '/#ifdef SEARCH_ENABLED/,/#else/d' \
        -e '/#endif \/\/ SEARCH_ENABLED/d' \
        < sequest-search.cpp > $(distdir)/sequest-search.cpp; fi
