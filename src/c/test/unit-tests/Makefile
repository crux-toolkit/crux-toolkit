# NOTE: When debugging unit tests using the check library you usually
# want to set the environment variable CK_FORK to no. This keeps the
# check library from running the tests on multiple threads which makes
# them hard to debug.
BINDIR  =../../../../bin/
SRCDIR  =../../
PERCOLATOR_DIR=../../../c/percolator
CC	= g++	
INTERNAL_CXXFLAGS  = -g -Wall -I../../ -I$(PERCOLATOR_DIR) -Icheck -Lcheck $(CXXFLAGS)
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
	check-crux-utils.c \
	crux-unit-tests.c

CRUX_LIBS = \
  ../../.libs/libcrux.a \
  $(PERCOLATOR_DIR)/.libs/libpercolator.a \
  ../../MSToolkit/.libs/libmstoolkit.a

unittests:  $(SRC) $(TESTS) $(CRUX_LIBS)
	cd check && $(MAKE)
	$(CC) -o unittests $(INTERNAL_CXXFLAGS) $(TESTS) $(CRUX_LIBS) check/check.ar  -lm -D_GNU_SOURCE
.PHONY : check clean

check: unittests
	./unittests

clean:
	cd check  && $(MAKE) clean
	rm -f *.o unittests

