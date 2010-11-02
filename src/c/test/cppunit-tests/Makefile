CC        = g++  
CFLAGS    = -Icppunit-1.12.1/include -I../..
CRUX_LIB  = ../../.libs/libcrux.a
MSTOOLKIT_LIB = ../../MSToolkit/.libs/libmstoolkit.a
UNIT_LIB  = cppunit-1.12.1/src/cppunit/.libs/libcppunit.a 
LDFLAGS   = -ldl

TESTS = \
	main.cpp \
	TestModifications.cpp \
	TestXml.cpp \
        TestSpectrum.cpp \
        TestMatchFileReader.cpp \
        TestDelimitedFileWriter.cpp \
        TestMatchFileWriter.cpp

unittests: $(TESTS) $(CRUX_LIB) $(MSTOOLKIT_LIB) $(UNIT_LIB)  
	$(CC) -o unittests $(CFLAGS) $(TESTS) $(CRUX_LIB) $(MSTOOLKIT_LIB) $(UNIT_LIB) $(LDFLAGS)

$(UNIT_LIB):
	cd cppunit-1.12.1; ./configure; make

check: unittests
	./unittests

clean:
	rm -f *o unittests
