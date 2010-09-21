CC        = g++  
CFLAGS    = -Icppunit-1.12.1/include -I../..
CRUX_LIB  = ../../.libs/libcrux.a
UNIT_LIB  = cppunit-1.12.1/src/cppunit/.libs/libcppunit.a 
LDFLAGS   = -ldl

TESTS = \
	main.cpp \
	TestModifications.cpp

unittests: $(TESTS) $(CRUX_LIB) $(UNIT_LIB) 
	$(CC) -o unittests $(CFLAGS) $(TESTS) $(CRUX_LIB) $(UNIT_LIB) $(LDFLAGS)

$(UNIT_LIB):
	cd cppunit-1.12.1; ./configure; make

check: unittests
	./unittests

clean:
	rm -f *o unittests
