CFLAGS= -g -Wall
SRC =  ../utils.c ../peak.c ../spectrum.c ../peptide.c ../protein.c ../protein_peptide_association.c ../mass.c unit-tests.c check-peptide.c
unittests:  $(SRC)
	$(CC) -o unittests $(CFLAGS) $(SRC) -lcheck -lm -D_GNU_SOURCE

.PHONY : check clean

check: unittests
	./unittests

clean:
	rm -f *.o unittests
