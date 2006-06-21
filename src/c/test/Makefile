
SRC = peak.c peak-unit-tests.c check-peak.c utils.c
unittests:  $(SRC)
	$(CC) -o unittests $(CFLAGS) $(SRC) -lcheck -lm  

.PHONY : check clean

check: unittests
	./unittests

clean:
	rm -f *.o unittests
