# NOTE: When debugging unit tests using the check library you usually
# want to set the environment variable CK_FORK to no. This keeps the
# check library from running the tests on multiple threads which makes
# them hard to debug.
BINDIR  =../../bin/
SRCDIR  =../../src/
PERCOLATOR_DIR=../../src/percolator
PWIZ_DIR=../../external/proteowizard/install/
CC	= g++	
INTERNAL_CXXFLAGS  = -g -Wall -I../../src/ -I$(PERCOLATOR_DIR) -I$(PWIZ_DIR)/include -Icheck -Lcheck $(CXXFLAGS)
DEBUG   = -DNDEBUG -DBOUNDS_CHECK 
DEBUG   = -DDEBUG -g

TESTS = \
	check-crux-utils.c \
	check-database.c \
	check-generate-peptides-iterator.c \
	check-index2.c \
	check-ion_series.c \
	check-linked-list.c \
	check-mass.c \
	check-match-collection.c \
	check-modifications.c \
	check-modified-peptides-iterator.c \
	check-parameter.c \
	check-peak.c \
	check-peptide-modifications.c \
	check-peptide.c \
	check-protein.c \
	check-scorer.c \
	crux-unit-tests.c
	#check-match.c \
	#check-serialize.c \
	#check-qvalue.c \

CRUX_LIBS = \
  ../../src/.libs/libcrux.a \
  ../../ext/MSToolkit/.libs/libmstoolkit.a

PWIZ_LIBS = \
  -L$(PWIZ_DIR)/lib \
  -lpwiz_data_msdata \
  -lpwiz_data_common \
  -lpwiz_data_identdata \
  -lpwiz_data_msdata_version \
  -lpwiz_data_misc \
  -lpwiz_utility_misc \
  -lpwiz_utility_minimxml \
  -lpwiz_data_identdata_version \
  -lboost_system-mt-s \
  -lboost_iostreams-mt-s \
  -lz-mt-s-1_2 \
  -lboost_filesystem-mt-s \
  -lboost_regex-mt-s \
  -lboost_thread-mt-s \
  -lpthread



unittests:  $(SRC) $(TESTS) $(CRUX_LIBS)
	cd check && $(MAKE)
	$(CC) -o unittests $(INTERNAL_CXXFLAGS) $(TESTS) $(CRUX_LIBS) $(PWIZ_LIBS) check/check.ar  -lm -lpthread -ldl -D_GNU_SOURCE
.PHONY : check clean

check: unittests
	./unittests

clean:
	cd check  && $(MAKE) clean
	rm -f *.o unittests

