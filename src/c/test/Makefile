CFLAGS= -g -Wall
SRC =  ../utils.c ../peak.c ../spectrum.c unit-tests.c check-spectrum.c
unittests:  $(SRC)
	$(CC) -o unittests $(CFLAGS) $(SRC) -lcheck -lm -D_GNU_SOURCE

.PHONY : check clean

check: unittests
	./unittests

clean:
	rm -f *.o unittests
