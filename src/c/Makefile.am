distdir = $(PACKAGE)_$(VERSION)

SUBDIRS = ArrayLibrary percolator xlink posterior-error qranker-barista hardklor bullseye .

AM_CPPFLAGS = -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 \
 -D_GNU_SOURCE -DGCC -DHAVE_EXPAT_CONFIG_H
AM_CFLAGS = -Wall -W -Werror
AM_CXXFLAGS = -Wall -W 

# Crux can use either single or double precision floating poinnt.
if USEDOUBLES
  AM_CPPFLAGS += -DUSE_DOUBLES
endif

# Compiling under OS X requires a preprocessor flag.
if DARWIN
  AM_CPPFLAGS += -DDARWIN
endif

# Compiling under OS X pre-Lion requires another preprocessor flag.
if PRE_LION
  AM_CPPFLAGS += -DPRE_LION
endif

# We disable search and exclude search related code
# in open-source distribution of Crux.
if SEARCH
  AM_CPPFLAGS += -DSEARCH_ENABLED
  SEARCH_ENABLED = 1
endif

LIBFILES= \
  Alphabet.cpp \
  Alphabet.h \
  analyze_psms.cpp \
  analyze_psms.h \
  carp.cpp \
  carp.h \
  ComputeQValues.cpp \
  ComputeQValues.h \
  ComputeQValuesLegacy.h \
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
  GetMs2Spectrum.h\
  GetMs2Spectrum.cpp\
  hash.cpp \
  hash.h \
  Index.cpp \
  Index.h \
  IndexFile.cpp \
  IndexFile.h \
  IndexMap.cpp \
  IndexMap.h \
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
  MatchCollectionParser.cpp \
  MatchCollectionParser.h \
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
  MzIdentMLReader.cpp \
  MzIdentMLReader.h \
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
  PepXMLReader.cpp \
  PepXMLReader.h \
  PepXMLWriter.cpp \
  PepXMLWriter.h \
  Percolator.cpp \
  Percolator.h \
  GeneratePeptides.cpp \
  GeneratePeptides.h \
  PredictPeptideIons.cpp \
  PredictPeptideIons.h \
  PrintProcessedSpectra.cpp \
  PrintProcessedSpectra.h \
  PrintVersion.cpp \
  PrintVersion.h \
  Protein.cpp \
  Protein.h \
  PostProcessProtein.cpp \
  PostProcessProtein.h \
  ProteinPeptideIterator.cpp \
  ProteinPeptideIterator.h \
  ProteinIndex.cpp \
  ProteinIndex.h \
  PWIZSpectrumCollection.cpp \
  PWIZSpectrumCollection.h \
  ProteinIndexIterator.cpp \
  ProteinIndexIterator.h \
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
  xlink_compute_qvalues.h \
  WinCrux.h

noinst_LTLIBRARIES = libcrux.la
libcrux_la_CFLAGS = \
  $(AM_CFLAGS) \
  -I$(top_srcdir)/src/external/MSToolkit/include \
  -I$(top_srcdir)/src/external/proteowizard/install/include/ \
  -I$(top_srcdir)/src/c/percolator/ \
  -I$(top_srcdir)/src/c/qranker-barista/ \
  -I$(top_srcdir)/src/c/posterior-error/

libcrux_la_CXXFLAGS = \
  $(AM_CXXFLAGS) \
  -I$(top_srcdir)/src/external/MSToolkit/include \
  -I$(top_srcdir)/src/external/MSToolkit/RAMP \
  -I$(top_srcdir)/src/external/proteowizard/install/include/ \
  -I$(top_srcdir)/src/c/percolator/ \
  -I$(top_srcdir)/src/c/qranker-barista/ \
  -I$(top_srcdir)/src/c/posterior-error/ \
  -I$(top_srcdir)/src/c/hardklor/ \
  -I$(top_srcdir)/src/c/bullseye/ \
  -I$(top_srcdir)/src/external/

libcrux_la_SOURCES = $(LIBFILES)

PWIZ_LIB_FILES = \
  -lpwiz_data_msdata \
  -lpwiz_data_common \
  -lpwiz_data_identdata \
  -lpwiz_data_msdata_version \
  -lpwiz_data_misc \
  -lpwiz_utility_chemistry \
  -lpwiz_utility_misc \
  -lpwiz_utility_minimxml \
  -lpwiz_data_common \
  -lpwiz_data_identdata \
  -lpwiz_data_identdata_version \
  -lpwiz_data_proteome \
  -lboost_system-mt-s \
  -lboost_iostreams-mt-s \
  -lz-mt-s-1_2 \
  -lboost_filesystem-mt-s \
  -lboost_regex-mt-s \
  -lboost_thread-mt-s \
  -lpthread


PROGRAM_FILES = \
  crux \
  xlink-predict-peptide-ions \
  xlink-score-spectrum \
  xlink-assign-ions

bin_PROGRAMS = $(PROGRAM_FILES)

crux_SOURCES = crux-main.cpp crux-main.h
crux_CXXFLAGS = \
  $(AM_CXXFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/external/MSToolkit \
  -I$(top_srcdir)/src/external/proteowizard/install/include/ \
  -I$(top_srcdir)/src/c/percolator \
  -I$(top_srcdir)/src/c/posterior-error \
  -I$(top_srcdir)/src/c/xlink \
  -I$(top_srcdir)/src/c/qranker-barista \
  -I$(top_srcdir)/src/c/hardklor \
  -I$(top_srcdir)/src/c/bullseye

crux_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/c/xlink \
  -lxlink \
  libcrux.la \
  -lstdc++ \
  -L$(top_builddir)/src/c/ArrayLibrary \
  -L$(top_builddir)/src/external/MSToolkit \
  -L$(top_builddir)/src/c/percolator \
  -L$(top_builddir)/src/c/posterior-error \
  -L$(top_builddir)/src/c/qranker-barista \
  -L$(top_builddir)/src/c/posterior-error \
  -L$(top_builddir)/src/c/hardklor \
  -L$(top_builddir)/src/c/bullseye \
  -lpercolator \
  -lqranker_barista \
  -lhardklor \
  -lbullseye \
  -lpep \
  -lmstoolkit \
  -lm -lpthread -ldl \
  -lArrayLibrary \
  -L$(top_builddir)/src/external/proteowizard/install/lib \
  $(PWIZ_LIB_FILES)

xlink_predict_peptide_ions_SOURCES=./xlink/xhhc_predict_peptide_ions.cpp
xlink_predict_peptide_ions_CXXFLAGS= \
  $(AM_CFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/external/proteowizard/install/include/ \
  -I$(top_srcdir)/src/c/xlink

xlink_predict_peptide_ions_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/external/MSToolkit \
  -L$(top_builddir)/src/c/xlink \
  -lm \
  -lxlink \
  libcrux.la \
  -L$(top_builddir)/src/external/proteowizard/install/lib $(PWIZ_LIB_FILES) \
  -lmstoolkit

xlink_score_spectrum_SOURCES=./xlink/xhhc_score_peptide_spectrum.cpp
xlink_score_spectrum_CXXFLAGS= \
  $(AM_CFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/external/proteowizard/install/include/ \
  -I$(top_srcdir)/src/c/xlink

xlink_score_spectrum_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/external/MSToolkit \
  -L$(top_builddir)/src/c/xlink \
  -lm -lpthread -ldl \
  -lxlink \
  libcrux.la \
  -L$(top_builddir)/src/external/proteowizard/install/lib $(PWIZ_LIB_FILES) \
  -lmstoolkit

xlink_assign_ions_SOURCES=./xlink/xlink_assign_ions.cpp
xlink_assign_ions_CXXFLAGS= \
  $(AM_CFLAGS) \
  -DMAIN \
  -I$(top_srcdir)/src/external/proteowizard/install/include/ \
  -I$(top_srcdir)/src/c/xlink

xlink_assign_ions_LDADD = \
  $(AM_LDFLAGS) \
  -L$(top_builddir)/src/external/MSToolkit \
  -L$(top_builddir)/src/c/xlink \
  -lm -lpthread -ldl \
  -lxlink \
  libcrux.la \
  -L$(top_builddir)/src/external/proteowizard/install/lib $(PWIZ_LIB_FILES) \
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
