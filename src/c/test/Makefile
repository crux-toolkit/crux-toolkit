BINDIR  =../../../bin/
SRCDIR  =../
PERCOLATOR_DIR=../../../../percolator/src
CC	= gcc	
CFLAGS  = -g -Wall -I$(PERCOLATOR_DIR)
DEBUG   = -DNDEBUG -DBOUNDS_CHECK 
DEBUG   = -DDEBUG -g

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
#	check-parameter.c \
#	check-modifications.c \
#	check-peptide-modifications.c \
#	check-linked-list.c \
#	check-generate-peptides-iterator.c \
#	check-modified-peptides-iterator.c \

CRUX_LIBS = \
  ../.libs/libcrux.a \
  $(PERCOLATOR_DIR)/Release/lib/libpercolator.a

unittests:  $(SRC) $(TESTS)
	$(CC) -o unittests $(CFLAGS) $(TESTS) $(CRUX_LIBS) -lcheck -lm -D_GNU_SOURCE
.PHONY : check clean

check: unittests
	./unittests

clean:
	rm -f *.o unittests

