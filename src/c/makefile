#######################################################################
# FILE: makefile
# AUTHOR: William Stafford Noble
# VERSION: $Revision: 1.33 $
#######################################################################

# The GNU C compiler is preferred.
CC	= gcc
DEBUG   = -DNDEBUG -DBOUNDS_CHECK 
DEBUG   = -DDEBUG -g
OFLAGS  = #-O2 
CFLAGS = $(OFLAGS) $(DEBUG) -Wall -W -Werror \
         -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

##### Global variable definitions

# Location of binary files.
BINDIR = ../../bin

# Objects directory
OBJDIR = objs

# Library directory
LIBDIR = ../../lib

# Archive file.
LIBFILE = ${LIBDIR}/libcrux.a
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
	${OBJDIR}/create_index.o \
	${OBJDIR}/predict_peptide_ions.o

# These are the names of the programs.
PROGS = \
	${BINDIR}/get_ms2_spectrum \
	${BINDIR}/generate_peptides \
	${BINDIR}/create_peptide_index \
	${BINDIR}/create_index \
	${BINDIR}/predict_peptide_ions
#	${BINDIR}/charge-vector \
#	${BINDIR}/charge-vector-ms2 \
#	${BINDIR}/precursor

all: libs bins TAGS 
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
	makedepend -Y -p'$${OBJDIR}/' -f .depend *.c *.h

# Build a tags file for emacs.
TAGS: *.c *.h
	etags *.c *.h
	ctags *.c *.h

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
	${CC} ${CFLAGS} -c -o ${OBJDIR}/$*.o $<

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
	${CC} -DMAIN ${CFLAGS} -o ${BINDIR}/$* $*.c ${LIBS}
