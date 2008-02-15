distdir = $(PACKAGE)_$(VERSION)
CFLAGS = -Wall -W -Werror  -fno-strict-aliasing \
         -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

noinst_LTLIBRARIES = libcrux.la

libcrux_la_SOURCES = \
  alphabet.c \
  alphabet.h \
  array.c \
  array.h \
  carp.c \
  carp.h \
  crux-utils.c \
  crux-utils.h \
  database.c \
  database.h \
  generate_peptides_iterator.c \
  generate_peptides_iterator.h \
  hash.c \
  hash.h \
  index.c \
  index.h \
  ion.c \
  ion.h \
  ion_series.c \
  ion_series.h \
  mass.c \
  mass.h \
  match.c \
  match_collection.c \
  match_collection.h \
  match.h \
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
  peptide_src.c \
  peptide_src.h \
  protein.c \
  protein.h \
  protein_index.c \
  protein_index.h \
  scorer.c \
  scorer.h \
  sorter.c \
  sorter.h \
  spectrum.c \
  spectrum_collection.c \
  spectrum_collection.h \
  spectrum.h \
  utils.c \
  utils.h

libcrux_la_CFLAGS = $(AM_CFLAGS) -I$(PERCOLATOR_DIR)/src
libcrux_la_LIBADD = -lm $(PERCOLATOR_DIR)/src/Release/lib/libpercolator.a

bin_PROGRAMS = \
  crux-get-ms2-spectrum \
  crux-generate-peptides \
  crux-create-index \
  crux-search-for-matches \
  crux-analyze-matches \
  crux-predict-peptide-ions \
  create_protein_index \
  transform_ms2_spectra \
  create_psm_files

crux_get_ms2_spectrum_SOURCES = get_ms2_spectrum.c
crux_get_ms2_spectrum_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_get_ms2_spectrum_LDADD = libcrux.la 

crux_generate_peptides_SOURCES = generate_peptides.c
crux_generate_peptides_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_generate_peptides_LDADD = libcrux.la 

crux_create_index_SOURCES = create_index.c
crux_create_index_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_create_index_LDADD = libcrux.la 

crux_search_for_matches_SOURCES = match_search.c
crux_search_for_matches_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_search_for_matches_LDADD = libcrux.la 

crux_analyze_matches_SOURCES = match_analysis.c
crux_analyze_matches_CFLAGS = $(AM_CFLAGS) -DMAIN -I$(PERCOLATOR_DIR)/src
crux_analyze_matches_LDADD = \
  -lstdc++ $(PERCOLATOR_DIR)/src/Release/lib/libpercolator.a libcrux.la

crux_predict_peptide_ions_SOURCES = predict_peptide_ions.c
crux_predict_peptide_ions_CFLAGS = $(AM_CFLAGS) -DMAIN 
crux_predict_peptide_ions_LDADD = libcrux.la 

create_protein_index_SOURCES = create_protein_index.c
create_protein_index_CFLAGS = $(AM_CFLAGS) -DMAIN 
create_protein_index_LDADD = libcrux.la 

transform_ms2_spectra_SOURCES = transform_ms2_spectra.c
transform_ms2_spectra_CFLAGS = $(AM_CFLAGS) -DMAIN 
transform_ms2_spectra_LDADD = libcrux.la 

create_psm_files_SOURCES = create_psm_files.c
create_psm_files_CFLAGS = $(AM_CFLAGS) -DMAIN 
create_psm_files_LDADD = libcrux.la 

doc:	*.h *.c
	doxygen Doxyfile
