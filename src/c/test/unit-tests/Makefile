BINDIR  =../../../../bin/
SRCDIR  =../../
PERCOLATOR_DIR=../../../c/percolator
CC	= g++	
CFLAGS  = -g -Wall -I../../ -I$(PERCOLATOR_DIR) -I$(HOME)/local/include -L$(HOME)/local/lib
DEBUG   = -DNDEBUG -DBOUNDS_CHECK 
DEBUG   = -DDEBUG -g

TESTS = \
	check-peak.c \
	check-spectrum.c \
	check-spectrum_collection.c \
	check-ion.c \
	check-ion_series.c \
	check-index2.c \
	check-scorer.c \
	check-peptide.c \
	check-protein.c \
	check-database.c \
	check-match.c \
	check-match-collection.c \
	check-parameter.c \
	check-modifications.c \
	check-peptide-modifications.c \
	check-linked-list.c \
	check-generate-peptides-iterator.c \
	check-modified-peptides-iterator.c \
	check-serialize.c \
	check-mass.c \
	crux-unit-tests.c

CRUX_LIBS = \
  ../../.libs/libcrux.a \
  $(PERCOLATOR_DIR)/.libs/libpercolator.a \
  ../../MSToolkit/.libs/libmstoolkit.a

unittests:  $(SRC) $(TESTS) $(CRUX_LIBS)
	$(CC) -o unittests $(CFLAGS) $(TESTS) $(CRUX_LIBS) -lcheck -lm -D_GNU_SOURCE
.PHONY : check clean

check: unittests
	./unittests

clean:
	rm -f *.o unittests

