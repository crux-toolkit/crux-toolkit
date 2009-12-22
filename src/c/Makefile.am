distdir = $(PACKAGE)_$(VERSION)

SUBDIRS = ArrayLibrary MSToolkit percolator q-ranker .

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
  alphabet.c \
  alphabet.h \
  array.c \
  array.h \
  carp.c \
  carp.h \
  create_index.c \
  create_index.h \
  crux-utils.c \
  crux-utils.h \
  database.c \
  database.h \
  generate_peptides_iterator.c \
  generate_peptides_iterator.h \
  hash.c \
  hash.h \
  hit.c \
  hit.h \
  index.c \
  index.h \
  ion.c \
  ion.h \
  ion_series.c \
  ion_series.h \
  linked_list.c \
  linked_list.h \
  mass.c \
  mass.h \
  match.c \
  match.h \
  match_collection.c \
  match_collection.h \
  match_search.cpp \
  modifications.c \
  modifications.h \
  modified_peptides_iterator.c \
  modified_peptides_iterator.h \
  search.h \
  objects.h \
  output-files.cpp \
  output-files.h \
  parameter.c \
  parameter.h \
  parse_arguments.c \
  parse_arguments.h \
  peak.c \
  peak.h \
  peptide.c \
  peptide_constraint.c \
  peptide_constraint.h \
  peptide.h \
  peptide_modifications.c \
  peptide_modifications.h \
  peptide_src.c \
  peptide_src.h \
  percolator.cpp \
  percolator.h \
  print-processed-spectra.cpp \
  print-processed-spectra.h \
  protein.c \
  protein.h \
  protein_index.c \
  protein_index.h \
  q-ranker.cpp \
  q-ranker.h \
  q-value.cpp \
  q-value.h \
  sequest-search.cpp \
  sequest-search.h \
  scorer.c \
  scorer.h \
  sorter.c \
  sorter.h \
  spectrum.c \
  spectrum_collection.c \
  spectrum_collection.h \
  spectrum.h \
  utils.c \
  utils.h \
  version.h

noinst_LTLIBRARIES = libcrux.la
libcrux_la_CFLAGS = \
  $(AM_CFLAGS) \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/percolator/ \
  -I$(top_srcdir)/src/c/q-ranker/

libcrux_la_CXXFLAGS = \
  $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/c/MSToolkit \
  -I$(top_srcdir)/src/c/percolator/ \
  -I$(top_srcdir)/src/c/q-ranker/

libcrux_la_SOURCES = $(LIBFILES)

PROGRAM_FILES = \
  crux \
  crux-get-ms2-spectrum \
  crux-generate-peptides \
  crux-predict-peptide-ions 

bin_PROGRAMS = $(PROGRAM_FILES)

crux_get_ms2_spectrum_SOURCES = get_ms2_spectrum.cpp
crux_get_ms2_spectrum_CXXFLAGS = $(AM_CFLAGS) -DMAIN -I$(top_srcdir)/src/c//MSToolkit
crux_get_ms2_spectrum_LDADD = \
  $(AM_LDFLAGS) \
  libcrux.la -L$(top_builddir)/src/c/MSToolkit -lm -lmstoolkit

crux_generate_peptides_SOURCES = generate_peptides.c
crux_generate_peptides_CFLAGS = $(AM_CFLAGS) -DMAIN
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
  -I$(top_srcdir)/src/c/q-ranker
crux_LDADD = \
  $(AM_LDFLAGS) \
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

paramdir = $(prefix)/doc
param_DATA = default.params

EXTRA_DIST = default.params

doc:	*.h *.c
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
        < match_search.c > $(distdir)/match_search.c; fi
