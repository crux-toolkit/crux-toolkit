CFLAGS= -g -Wall
SRC =	../utils.c ../peak.c ../spectrum.c ../peptide.c ../protein.c ../peptide_src.c ../mass.c ../alphabet.c ../array.c ../crux-utils.c ../database.c
BINDIR=../../../bin

test_generate_peptides: $(SRC) $(BINDIR)/generate_peptides
	( $(BINDIR)/generate_peptides test.fasta > generate_peptides.out) >& generate_peptides.err 
	diff generate_peptides.out standard_generate_peptides.out
	diff generate_peptides.err standard_generate_peptides.err

test_get_ms2_spectrum: $(SRC) $(BINDIR)/get_ms2_spectrum
	rm -f get_ms2_spectrum.out get_ms2_spectrum.err
	( $(BINDIR)/get_ms2_spectrum 2 test.ms2 get_ms2_spectrum.out) >& get_ms2_spectrum.err 
	diff get_ms2_spectrum.out standard_get_ms2_spectrum.out
	diff get_ms2_spectrum.err standard_get_ms2_spectrum.err

smoke_tests: test_generate_peptides test_get_ms2_spectrum

clean:
	rm -f *.o 
