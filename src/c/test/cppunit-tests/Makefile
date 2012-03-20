CC        = g++  
CFLAGS    = -Icppunit-1.12.1/include -I../.. -I../../qranker-barista
CRUX_LIB  = ../../.libs/libcrux.a
MSTOOLKIT_LIB = ../../../external/MSToolkit/.libs/libmstoolkit.a
BARISTA_LIB = ../../qranker-barista/.libs/libqranker_barista.a
PERCOLATOR_LIB = ../../percolator/.libs/libpercolator.a
PEP_LIB = ../../posterior-error/.libs/libpep.a
ARRAY_LIB = ../../ArrayLibrary/.libs/libArrayLibrary.a

UNIT_LIB  = cppunit-1.12.1/src/cppunit/.libs/libcppunit.a 
LDFLAGS   = -lpthread -ldl

TESTS = \
	main.cpp \
	TestModifications.cpp \
	TestXml.cpp \
        TestSpectrum.cpp \
        TestMatchFileReader.cpp \
        TestDelimitedFileWriter.cpp \
        TestMatchFileWriter.cpp \
	TestProtein.cpp

unittests: $(TESTS) $(CRUX_LIB) $(MSTOOLKIT_LIB) $(UNIT_LIB)  
	$(CC) -o unittests $(CFLAGS) $(TESTS) $(CRUX_LIB) $(MSTOOLKIT_LIB) $(BARISTA_LIB) $(PERCOLATOR_LIB) $(PEP_LIB) $(ARRAY_LIB) $(UNIT_LIB) $(LDFLAGS)

$(UNIT_LIB):
	cd cppunit-1.12.1; ./configure; make

check: unittests
	./unittests

clean:
	rm -f *o unittests
