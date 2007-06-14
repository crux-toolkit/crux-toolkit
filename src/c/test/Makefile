CFLAGS= -g -Wall
DEBUG   = -DNDEBUG -DBOUNDS_CHECK 
DEBUG   = -DDEBUG -g
CC	= gcc	
BINDIR=../../../bin/
SRCDIR=../
SRC = \
	../utils.c \
	../peak.c \
	../spectrum.c \
	../peptide.c \
	../protein.c \
	../peptide_src.c \
	../mass.c \
	../alphabet.c \
	../array.c \
	../crux-utils.c \
	../database.c \
	../carp.c \
	../peptide_constraint.c \
	 ../index.c \
	../sorter.c \
	../spectrum_collection.c \
	../protein_index.c \
	../ion.c \
	../ion_series.c \
	../parameter.c \
	../scorer.c \
	../parse_arguments.c \
	../generate_peptides_iterator.c \
	../match.c \
	../match_collection.c
TESTS = \
	check-peak.c \
	check-spectrum.c \
	check-spectrum_collection.c \
	check-ion.c \
	check-ion_series.c \
	check-index.c \
	check-scorer.c \
	check-peptide.c \
	check-protein.c \
	check-database.c \
	check-match.c \
	crux-unit-tests.c


unittests:  $(SRC) $(TESTS)
	$(CC) -o unittests $(CFLAGS) $(SRC) $(TESTS) -lcheck -lm -D_GNU_SOURCE
.PHONY : check clean

check: unittests
	./unittests

clean:
	rm -f *.o unittests

