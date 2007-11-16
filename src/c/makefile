#######################################################################
# FILE: makefile
# AUTHOR: William Stafford Noble
# VERSION: $Revision: 1.96 $   #-Werror \
#######################################################################

# The GNU C compiler is preferred.
CC	= gcc
DEBUG   = # -DNDEBUG -DBOUNDS_CHECK 
DEBUG   = -DDEBUG -g
DEBUG   = # -DDEBUG -pg
OFLAGS  = # -O2 
CFLAGS = $(OFLAGS) $(DEBUG) -Wall -W -Werror \
         -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
CPP_FLAGS =  -lstdc++
OS_FLAGS = -DCYGWIN -fno-strict-aliasing

##### Different module directories ####
PERCOLATOR_DIR = ../../../percolator

##### Global variable definitions

# Location of binary files.
BINDIR = ../../bin

# Objects directory
OBJDIR = objs

# Library directory
LIBDIR = ../../lib

# Archive file.
LIBFILE = ${LIBDIR}/libcrux.a
PERCOLATOR_LIB = ${PERCOLATOR_DIR}/src/Release/lib/libpercolator.a
RANLIBFILE = ${LIBDIR}/libcrux.ranlib

# Libraries to use
LIBS = ${LIBFILE} -lm
#LIBS = ${LIBFILE} -lm -lefence

# Objects to compile into the library
LIBOBJS =\
	${OBJDIR}/utils.o \
	${OBJDIR}/peak.o \
	${OBJDIR}/spectrum.o \
	${OBJDIR}/peptide.o \
	${OBJDIR}/protein.o \
	${OBJDIR}/protein_index.o \
	${OBJDIR}/peptide_src.o \
	${OBJDIR}/peptide_constraint.o \
	${OBJDIR}/database.o \
	${OBJDIR}/spectrum_collection.o \
	${OBJDIR}/parse_arguments.o \
	${OBJDIR}/crux-utils.o \
	${OBJDIR}/array.o \
	${OBJDIR}/alphabet.o \
	${OBJDIR}/ion.o \
	${OBJDIR}/ion_series.o \
	${OBJDIR}/carp.o \
	${OBJDIR}/mass.o \
	${OBJDIR}/index.o \
	${OBJDIR}/sorter.o \
	${OBJDIR}/get_ms2_spectrum.o \
	${OBJDIR}/generate_peptides.o \
	${OBJDIR}/generate_peptides_iterator.o \
	${OBJDIR}/create_protein_index.o \
	${OBJDIR}/create_index.o \
	${OBJDIR}/predict_peptide_ions.o \
	${OBJDIR}/scorer.o \
	${OBJDIR}/parameter.o \
	${OBJDIR}/score_peptide_spectrum.o \
	${OBJDIR}/match.o \
	${OBJDIR}/match_collection.o \
	${OBJDIR}/hash.o \
	${OBJDIR}/match_search.o \
	${OBJDIR}/match_analysis.o
#	${OBJDIR}/score_peptide_list_spectrum.o

# These are the names of the programs.
PROGS = \
	${BINDIR}/crux-get-ms2-spectrum \
	${BINDIR}/crux-generate-peptides \
	${BINDIR}/crux-create-index  \
	${BINDIR}/crux-search-for-matches \
	${BINDIR}/crux-analyze-matches \
	${BINDIR}/crux-predict-peptide-ions \
	${BINDIR}/crux-score-peptide-spectrum \
	${BINDIR}/create_psm_files \
	${BINDIR}/create_ion_files \
	${BINDIR}/utils \
	${BINDIR}/create_protein_index \
	${BINDIR}/transform_ms2_spectra
#	${BINDIR}/get_ms2_spectrum \
#	${BINDIR}/generate_peptides \
#	${BINDIR}/create_index  \
#	${BINDIR}/match_search \
#	${BINDIR}/match_analysis \
#	${BINDIR}/predict_peptide_ions \
#	${BINDIR}/score_peptide_spectrum \
#	${BINDIR}/score_peptide_list_spectrum
#	${BINDIR}/charge-vector \
#	${BINDIR}/charge-vector-ms2 \
#	${BINDIR}/precursor

all: srclinks Percolator libs bins binlinks TAGS 
most: srclinks Percolator libs bins binlinks 
libs: ${LIBFILE}(${LIBOBJS})

bins: ${PROGS} libs

python: ${PYTHON}

doc:	*.h *.c
	doxygen Doxyfile

# Remove extraneous files.
clean: 
	rm -rf ${OBJDIR} ${LIBDIR} ${BINDIR} *.bak ${PROGS} ${LIBFILE} ${RANLIBFILE} core html latex

# Build dependencies.  Note that -p argument is quoted for both
# make and shell so ${OBJDIR}/ is inserted literally.
depend: 
	rm -f .depend
	touch .depend
	makedepend -Y -p'$${OBJDIR}/' -p '$(PERCOLATOR_DIR)/src/Release/' -f .depend *.c *.h *.cpp

# Build a tags file for emacs.
TAGS: *.c *.h *.cpp
	etags *.c *.h
	ctags *.c *.h

# Make Percolator and copy over header and library file to CRUX
Percolator:
	# make percolator 
	@cd $(PERCOLATOR_DIR)/src ; make
	@cp $(PERCOLATOR_DIR)/src/PercolatorCInterface.h  . 
	@mkdir -p ${LIBDIR}
	@cp $(PERCOLATOR_DIR)/src/Release/lib/libpercolator.a ${LIBDIR} 

##### Dependencies

${PROGS}: ${RANLIBFILE}

${RANLIBFILE}: ${LIBFILE}

include .depend

##### Rules

# Undefine some built-in rules; especially avoid the archive target
# (and other targets) that removes .o files.
(%): %
(%.o): %.o
(%.o): %.cc
(%.o): %.cpp
(%.o): %.c
%.o: %.c
%: %.o
%: %.c

# Rule for compiling C files into .o files.
${OBJDIR}/%.o: %.c
	@mkdir -p ${OBJDIR}
	${CC} ${CFLAGS} ${OS_FLAGS} -c -o ${OBJDIR}/$*.o $<

#################
${OBJDIR}/%.o: %.cpp
	@mkdir -p ${OBJDIR}
	${CC} ${CFLAGS} ${OS_FLAGS} -c -o ${OBJDIR}/$*.o $<
##################

# Rule for adding .o files to an archive library
(%.o): %.o
	@mkdir -p ${LIBDIR}
	${AR} r $@ $*.o

# Rule for doing ranlib.
%.ranlib: %.a
	ranlib $<
	touch $@

# Rules for compiling programs
${BINDIR}/%: ${OBJDIR}/%.o ${RANLIBFILE}
	@mkdir -p ${BINDIR}
	${CC} -DMAIN ${CFLAGS} ${OS_FLAGS} -o ${BINDIR}/$* $*.c ${LIBS} ${PERCOLATOR_LIB}  ${CPP_FLAGS}

srclinks:
	#make source links
	@if [ ! -e crux-create-index.c ] ; then ln -s create_index.c crux-create-index.c; fi
	@if [ ! -e crux-get-ms2-spectrum.c ] ; then ln -s get_ms2_spectrum.c crux-get-ms2-spectrum.c; fi
	@if [ ! -e crux-generate-peptides.c ]; then ln -s generate_peptides.c crux-generate-peptides.c; fi
	@if [ ! -e crux-search-for-matches.c ]; then ln -s match_search.c crux-search-for-matches.c; fi
	@if [ ! -e crux-analyze-matches.c ]; then ln -s match_analysis.c crux-analyze-matches.c; fi
	@if [ ! -e crux-predict-peptide-ions.c ]; then ln -s predict_peptide_ions.c crux-predict-peptide-ions.c; fi
	@if [ ! -e crux-score-peptide-spectrum.c ]; then ln -s score_peptide_spectrum.c crux-score-peptide-spectrum.c; fi

binlinks: 
	# make links to old executable names
	@rm -f ${BINDIR}/create_index
	@ln -s crux-create-index ${BINDIR}/create_index
	@rm -f ${BINDIR}/get_ms2_spectrum
	@ln -s crux-get-ms2-spectrum ${BINDIR}/get_ms2_spectrum
	@rm -f ${BINDIR}/generate_peptides
	@ln -s crux-generate-peptides ${BINDIR}/generate_peptides
	@rm -f ${BINDIR}/match_search
	@ln -s crux-search-for-matches ${BINDIR}/match_search
	@rm -f ${BINDIR}/match_analysis
	@ln -s crux-analyze-matches ${BINDIR}/match_analysis
	@rm -f ${BINDIR}/predict_peptide_ions
	@ln -s crux-predict-peptide-ions ${BINDIR}/predict_peptide_ions
	@rm -f ${BINDIR}/score_peptide_spectrum
	@ln -s crux-score-peptide-spectrum ${BINDIR}/score_peptide_spectrum
