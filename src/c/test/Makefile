CFLAGS= -g -Wall
BINDIR=../../../bin/
SRC =  ../utils.c ../peak.c ../spectrum.c ../peptide.c ../protein.c ../peptide_src.c ../mass.c ../alphabet.c ../array.c  check-index.c ../crux-utils.c ../database.c ../carp.c ../peptide_constraint.c ../index.c crux-unit-tests.c
unittests:  $(SRC)
	$(CC) -o unittests $(CFLAGS) $(SRC) -lcheck -lm -D_GNU_SOURCE
.PHONY : check clean

check: unittests
	./unittests

clean:
	rm -f *.o unittests
