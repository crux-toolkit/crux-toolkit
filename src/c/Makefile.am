distdir = $(PACKAGE)_$(VERSION)

SUBDIRS = percolator .

if USEDOUBLES
  USEDOUBLES='-DUSE_DOUBLES'
else
  USEDOUBLES=
endif

if DARWIN
  USEDARWIN = '-DDARWIN'
else
  USEDARWIN = 
endif

AM_CFLAGS = -Wall -W -Werror  -fno-strict-aliasing $(USEDOUBLES) $(USEDARWIN) \
            -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

noinst_LTLIBRARIES = libcrux.la

NON_PERCOLATOR_LIB_FILES= \
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
  hit_collection.c \
  hit_collection.h \
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
  match_search.c \
  modifications.c \
  modifications.h \
  modified_peptides_iterator.c \
  modified_peptides_iterator.h \
  search.h \
  objects.h \
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
  protein.c \
  protein.h \
  protein_index.c \
  protein_index.h \
  q-value.c \
  q-value.h \
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

LIBFILES=$(PERCOLATOR_LIB_FILES) $(NON_PERCOLATOR_LIB_FILES)
libcrux_la_CFLAGS = $(AM_CFLAGS)

libcrux_la_SOURCES = $(LIBFILES)
libcrux_la_LIBADD = -lm

PROGRAM_FILES = \
  crux \
  crux-get-ms2-spectrum \
  crux-generate-peptides \
  crux-predict-peptide-ions 
#  crux-search-for-matches \
#  crux-create-index \
#  create_protein_index \
#  transform_ms2_spectra \
#  create_psm_files

bin_PROGRAMS = $(PROGRAM_FILES)

crux_get_ms2_spectrum_SOURCES = get_ms2_spectrum.c
crux_get_ms2_spectrum_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_get_ms2_spectrum_LDADD = libcrux.la 

crux_generate_peptides_SOURCES = generate_peptides.c
crux_generate_peptides_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_generate_peptides_LDADD = libcrux.la 

#crux_create_index_SOURCES = create_index.c
#crux_create_index_CFLAGS = $(AM_CFLAGS) -DMAIN 
#crux_create_index_LDADD = libcrux.la 

#crux_search_for_matches_SOURCES = match_search.c
#crux_search_for_matches_CFLAGS = $(AM_CFLAGS) -DMAIN -I$(PERCOLATOR_DIR)/
#crux_search_for_matches_LDADD = libcrux.la 

#crux_analyze_matches_SOURCES = match_analysis.c
#crux_analyze_matches_CFLAGS = $(AM_CFLAGS) -DMAIN -I$(PERCOLATOR_DIR)
#crux_analyze_matches_LDADD = \
#  -lstdc++ -L$(PERCOLATOR_DIR)/ -lpercolator -lqvality -lArrayLibrary -lcrux 

crux_predict_peptide_ions_SOURCES = predict_peptide_ions.c
crux_predict_peptide_ions_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_predict_peptide_ions_LDADD = libcrux.la 

#create_protein_index_SOURCES = create_protein_index.c
#create_protein_index_CFLAGS = $(AM_CFLAGS) -DMAIN 
#create_protein_index_LDADD = libcrux.la 

#transform_ms2_spectra_SOURCES = transform_ms2_spectra.c
#transform_ms2_spectra_CFLAGS = $(AM_CFLAGS) -DMAIN 
#transform_ms2_spectra_LDADD = libcrux.la 

#create_psm_files_SOURCES = create_psm_files.c
#create_psm_files_CFLAGS = $(AM_CFLAGS) -DMAIN 
#create_psm_files_LDADD = libcrux.la 

crux_SOURCES = crux-main.cpp  crux-main.h percolator.c percolator.h
crux_CFLAGS = $(AM_CFLAGS) -DMAIN -Ipercolator
crux_CPPFLAGS = $(AM_CFLAGS) -DMAIN -Ipercolator
crux_LDADD = \
  libcrux.la -lstdc++ -L./percolator -lpercolator -lqvality -lArrayLibrary

paramdir = $(prefix)/doc
param_DATA = default.params

EXTRA_DIST = default.params ../../bin/spectacle.pl

doc:	*.h *.c
	doxygen Doxyfile

ignore-files:
	svn update .svnignore
	svn propset svn:ignore -F .svnignore .

dist-hook:
	grep '#define VERSION' ../../config.h > version.h
