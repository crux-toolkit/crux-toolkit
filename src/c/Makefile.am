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
  Alphabet.cpp \
  Alphabet.h \
  analyze_psms.cpp \
  analyze_psms.h \
  carp.cpp \
  carp.h \
  ComputeQValues.cpp \
  ComputeQValues.h \
  CreateIndex.cpp \
  CreateIndex.h \
  crux-file-utils.h \
  crux-utils.cpp \
  crux-utils.h \
  Database.cpp \
  Database.h \
  DatabaseProteinIterator.cpp \
  DatabaseProteinIterator.h \
  DatabasePeptideIterator.cpp \
  DatabasePeptideIterator.h \
  CruxApplication.cpp \
  CruxApplication.h \
  CruxApplicationList.cpp \
  CruxApplicationList.h \
  DelimitedFile.cpp \
  DelimitedFile.h \
  DelimitedFileReader.cpp \
  DelimitedFileReader.h \
  DelimitedFileWriter.cpp \
  DelimitedFileWriter.h \
  ExtractColumns.cpp \
  ExtractColumns.h \
  ExtractRows.cpp \
  ExtractRows.h \
  FilteredSpectrumChargeIterator.cpp \
  FilteredSpectrumChargeIterator.h \
  GeneratePeptidesIterator.cpp \
  GeneratePeptidesIterator.h \
  hash.cpp \
  hash.h \
  Index.cpp \
  Index.h \
  IndexFile.cpp \
  IndexFile.h \
  IndexPeptideIterator.cpp \
  IndexPeptideIterator.h \
  IndexFilteredPeptideIterator.cpp \
  IndexFilteredPeptideIterator.h \
  Ion.cpp \
  Ion.h \
  IonConstraint.cpp \
  IonConstraint.h \
  IonFilteredIterator.cpp \
  IonFilteredIterator.h \
  IonSeries.cpp \
  IonSeries.h \
  linked_list.cpp \
  linked_list.h \
  mass.cpp \
  mass.h \
  Match.cpp \
  Match.h \
  MatchColumns.cpp \
  MatchColumns.h \
  MatchFileReader.cpp \
  MatchFileReader.h \
  MatchFileWriter.cpp \
  MatchFileWriter.h \
  MatchSearch.cpp \
  MatchSearch.h \
  MatchCollection.cpp \
  MatchCollection.h \
  MatchCollectionIterator.cpp \
  MatchCollectionIterator.h \
  MatchIterator.cpp \
  MatchIterator.h \
  MGFSpectrumCollection.cpp \
  MGFSpectrumCollection.h \
  modifications.cpp \
  modifications.h \
  ModifiedPeptidesIterator.cpp \
  ModifiedPeptidesIterator.h \
  MS2SpectrumCollection.cpp \
  MS2SpectrumCollection.h \
  MSToolkitSpectrumCollection.cpp \
  MSToolkitSpectrumCollection.h \
  objects.h \
  OutputFiles.cpp \
  OutputFiles.h \
  parameter.cpp \
  parameter.h \
  parse_arguments.cpp \
  parse_arguments.h \
  Peak.cpp \
  Peak.h \
  Peptide.cpp \
  Peptide.h \
  PeptideConstraint.cpp \
  PeptideConstraint.h \
  PeptideIterator.h \
  peptide_modifications.cpp \
  peptide_modifications.h \
  PeptideSrc.cpp \
  PeptideSrc.h \
  Percolator.cpp \
  Percolator.h \
  PrintProcessedSpectra.cpp \
  PrintProcessedSpectra.h \
  PrintVersion.cpp \
  PrintVersion.h \
  Protein.cpp \
  Protein.h \
  ProteinPeptideIterator.cpp \
  ProteinPeptideIterator.h \
  ProteinIndex.cpp \
  ProteinIndex.h \
  ProteinIndexIterator.cpp \
  ProteinIndexIterator.h \
  QRanker.cpp \
  QRanker.h \
  q-value.cpp \
  q-value.h \
  SearchProgress.h \
  SequestSearch.cpp \
  SequestSearch.h \
  Scorer.cpp \
  Scorer.h \
  SortColumn.cpp \
  SortColumn.h \
  sorter.cpp \
  sorter.h \
  SpectralCounts.cpp \
  SpectralCounts.h \
  SpectrumCollection.cpp \
  SpectrumCollection.h \
  SpectrumCollectionFactory.cpp \
  SpectrumCollectionFactory.h \
  Spectrum.cpp \
  Spectrum.h \
  SpectrumZState.cpp \
  SpectrumZState.h \
  StatColumn.cpp \
  StatColumn.h \
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
  -I$(top_srcdir)/src/c/q-ranker/ \
  -I$(top_srcdir)/src/external/

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
        < MatchSearch.cpp > $(distdir)/MatchSearch.cpp; \
    sed -e '/#ifdef SEARCH_ENABLED/,/#else/d' \
        -e '/#endif \/\/ SEARCH_ENABLED/d' \
        < SequestSearch.cpp > $(distdir)/SequestSearch.cpp; fi
