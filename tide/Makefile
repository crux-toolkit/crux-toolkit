# The directory system is as follows:
#
# Source files (.h and .cc) are in src/
# Source files (.proto) are in protobuf/
# Compiled protocol buffers (.pb.h and .pb.cc and _pb2.py) are in protoobj/
# Compiled object files (.o including .pb.o) are in obj/
# Executables compile to the top-level directory

export
ARCH = $(shell uname -m)

ifneq (y, $(DEBUG:Y=y))
# Default values
OPTIONS = -ggdb -O6 -DNDEBUG
OBJDIR = obj-opt-$(ARCH)
BINDIR = bin-opt-$(ARCH)
else
# Debug values
OPTIONS = -ggdb -O0
OBJDIR = obj-deb-$(ARCH)
BINDIR = bin-deb-$(ARCH)
endif

PROTOBUFS=raw_proteins peptides spectrum header results
EXECUTABLES= \
$(BINDIR)/test_spectrum_preprocess \
$(BINDIR)/peptide_sorter \
$(BINDIR)/peptide_downselect \
$(BINDIR)/make_peptides \
$(BINDIR)/filter_peptides \
$(BINDIR)/peptide_counter \
$(BINDIR)/search \
$(BINDIR)/test_read_ms2 \
$(BINDIR)/read_spectrumrecords \
$(BINDIR)/test_records \
$(BINDIR)/read_fasta \
$(BINDIR)/test_theoretical_peaks \
$(BINDIR)/test_fifo_alloc \
$(BINDIR)/read_peptide_diffs \
$(BINDIR)/peptide_mods \
$(BINDIR)/peptide_mods2 \
$(BINDIR)/test_mod_table \
$(BINDIR)/test_peptide_mods \
$(BINDIR)/peptide_peaks \
$(BINDIR)/index \
$(BINDIR)/results
#$(BINDIR)/peptide_mods3

OPTIONS+=-I./protoobj

ifeq (i686, $(ARCH))
LDFLAGS=-lprotobuf -lgflags -lpthread -lm
else
LDFLAGS=-lprotobuf -lgflags -lpthread -lm
endif

ifeq (y, $(STATIC:Y=y))
LDFLAGS := -static $(LDFLAGS)
BINDIR := static-$(BINDIR)
endif

.PHONY: all python_pb clean

all: python_pb $(EXECUTABLES)

python_pb: $(PROTOBUFS:%=protoobj/%_pb2.py)
	@touch protoobj/__init__.py

index: $(BINDIR)/index
$(BINDIR)/index: \
	$(OBJDIR)/index.o \
	$(OBJDIR)/index_settings.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/read_fasta.o \
	$(OBJDIR)/peptide.o \
	$(OBJDIR)/fifo_alloc.o \
	$(OBJDIR)/max_mz.o \
	$(OBJDIR)/make_peptides.o \
	$(OBJDIR)/peptide_mods3.o \
	$(OBJDIR)/peptide_peaks.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/mass_constants.o

read_fasta: $(BINDIR)/read_fasta
$(BINDIR)/read_fasta: \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/read_fasta.o \
	$(OBJDIR)/read_fasta_cmdline.o

test_spectrum_preprocess: $(BINDIR)/test_spectrum_preprocess
$(BINDIR)/test_spectrum_preprocess: \
	$(OBJDIR)/test_spectrum_preprocess.o \
	$(OBJDIR)/spectrum_preprocess2.o \
	$(OBJDIR)/max_mz.o \
	$(OBJDIR)/spectrum_collection.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/spectrum.pb.o

test_theoretical_peaks: $(BINDIR)/test_theoretical_peaks
$(BINDIR)/test_theoretical_peaks: \
	$(OBJDIR)/test_theoretical_peaks.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/peptide.o \
	$(OBJDIR)/fifo_alloc.o \
	$(OBJDIR)/mass_constants.o \
	$(OBJDIR)/max_mz.o \
	$(OBJDIR)/spectrum_preprocess2.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o	

test_fifo_alloc: $(BINDIR)/test_fifo_alloc
$(BINDIR)/test_fifo_alloc: \
	$(OBJDIR)/test_fifo_alloc.o \
	$(OBJDIR)/fifo_alloc.o

test_records: $(BINDIR)/test_records
$(BINDIR)/test_records: \
	$(OBJDIR)/test_records.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o	

make_peptides: $(BINDIR)/make_peptides
$(BINDIR)/make_peptides: \
	$(OBJDIR)/make_peptides.o \
	$(OBJDIR)/make_peptides_cmdline.o \
	$(OBJDIR)/index_settings.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/mass_constants.o

filter_peptides: $(BINDIR)/filter_peptides
$(BINDIR)/filter_peptides: \
	$(OBJDIR)/filter_peptides.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/peptides.pb.o

peptide_sorter: $(BINDIR)/peptide_sorter
$(BINDIR)/peptide_sorter: \
	$(OBJDIR)/peptide_sorter.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/peptides.pb.o

peptide_downselect: $(BINDIR)/peptide_downselect
$(BINDIR)/peptide_downselect: \
	$(OBJDIR)/peptide_downselect.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/peptides.pb.o

peptide_counter: $(BINDIR)/peptide_counter
$(BINDIR)/peptide_counter: \
	$(OBJDIR)/peptide_counter.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/peptides.pb.o

peptide_mods: $(BINDIR)/peptide_mods
$(BINDIR)/peptide_mods: \
	$(OBJDIR)/peptide_mods.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/mass_constants.o

peptide_mods2: $(BINDIR)/peptide_mods2
$(BINDIR)/peptide_mods2: \
	$(OBJDIR)/peptide_mods2.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/mass_constants.o

peptide_mods3: $(BINDIR)/peptide_mods3
$(BINDIR)/peptide_mods3: \
	$(OBJDIR)/peptide_mods3.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/mass_constants.o

test_mod_table: $(BINDIR)/test_mod_table
$(BINDIR)/test_mod_table: \
	$(OBJDIR)/test_mod_table.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/mass_constants.o

test_peptide_mods: $(BINDIR)/test_peptide_mods
$(BINDIR)/test_peptide_mods: \
	$(OBJDIR)/test_peptide_mods.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/mass_constants.o

peptide_peaks: $(BINDIR)/peptide_peaks
$(BINDIR)/peptide_peaks: \
	$(OBJDIR)/peptide_peaks_cmdline.o \
	$(OBJDIR)/peptide_peaks.o \
	$(OBJDIR)/peptide.o \
	$(OBJDIR)/fifo_alloc.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/max_mz.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/mass_constants.o

search: $(BINDIR)/search
$(BINDIR)/search: \
	$(OBJDIR)/search.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/spectrum.pb.o \
	$(OBJDIR)/results.pb.o \
	$(OBJDIR)/max_mz.o \
	$(OBJDIR)/peptide.o \
	$(OBJDIR)/spectrum_preprocess2.o \
	$(OBJDIR)/spectrum_collection.o \
	$(OBJDIR)/active_peptide_queue.o \
	$(OBJDIR)/fifo_alloc.o \
	$(OBJDIR)/mass_constants.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/report.o 

test_read_ms2: $(BINDIR)/test_read_ms2
$(BINDIR)/test_read_ms2: \
	$(OBJDIR)/test_read_ms2.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/spectrum.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/spectrum_collection.o

read_spectrumrecords: $(BINDIR)/read_spectrumrecords
$(BINDIR)/read_spectrumrecords: \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/spectrum.pb.o \
	$(OBJDIR)/spectrum_collection.o \
	$(OBJDIR)/read_spectrumrecords.o

read_peptide_diffs: $(BINDIR)/read_peptide_diffs
$(BINDIR)/read_peptide_diffs: \
	$(OBJDIR)/read_peptide_diffs.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/peptides.pb.o

results: $(BINDIR)/results
$(BINDIR)/results: \
	$(OBJDIR)/results.o \
	$(OBJDIR)/header.pb.o \
	$(OBJDIR)/peptides.pb.o \
	$(OBJDIR)/raw_proteins.pb.o \
	$(OBJDIR)/spectrum.pb.o \
	$(OBJDIR)/spectrum_collection.o \
	$(OBJDIR)/results_writer.o \
	$(OBJDIR)/results.pb.o \
	$(OBJDIR)/abspath.o \
	$(OBJDIR)/fifo_alloc.o \
	$(OBJDIR)/mass_constants.o \
	$(OBJDIR)/max_mz.o \
	$(OBJDIR)/peptide.o \
	$(OBJDIR)/spectrum_preprocess2.o \
	$(OBJDIR)/crux_sp_spectrum.o \
	$(OBJDIR)/sp_scorer.o


# protoc has to run within the protobuf directory, otherwise include
# directives in .proto files can't find the target file.
protoobj/%_pb2.py protoobj/%.pb.cc protoobj/%.pb.h: protobuf/%.proto
	@mkdir -p $(@D)
	cd protobuf && protoc --cpp_out=../protoobj --python_out=../protoobj $(<F)

clean:
	rm -f protoobj/*
	if test -n '$(OBJDIR)';then rm -f $(OBJDIR)/*; fi
	if test -n '$(BINDIR)';then rm -f $(BINDIR)/*; fi

# I'm not sure about the automatic cleaning part with NODEPS, etc.
# below. Right now, .d files persist after make is done, which I think
# is good. The .d files are needed to store the header
# dependencies. These are computed by g++ -MM (see below) and kept in
# a .d file for each .o file. g++ outputs a Makefile dependency
# clause, which we -include below. The .pb.o files don't participate
# in this, as it was tricky to try to make it work, and we know the
# dependencies anyway. So .pb.o has its own (fixed) set of
# dependencies, namely the correponding .pb.h and .pb.cc.

# Add .d to Make's recognized suffixes. (Not sure why/if this is needed.)
SUFFIXES += .d

#We don't need to clean up when we're making these targets
NODEPS:=clean
#Find all the C++ files in the src/ directory
SOURCES:=$(shell find src/ -name "*.cc")
#These are the dependency files, which make will clean up after it creates them
DEPFILES:=$(patsubst src/%.cc,$(OBJDIR)/%.d,$(SOURCES))

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
        #Chances are, these files don't exist.  GMake will create them and
        #clean up automatically afterwards
	-include $(DEPFILES)
endif

#This is the rule for creating the dependency files.  Because g++
#doesn't know at this point that the .pb.h files will be generated in
#the protoobj directory, we use the little inlined perl script in this
#command to prepend protoobj/ to any generated dependency ending in
#.pb.h that has no other specified directory (i.e. no / within).
$(OBJDIR)/%.d: src/%.cc
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -MM -MG -MT '$(patsubst src/%.cc,$(OBJDIR)/%.o,$<)' $< \
	|perl -p -e 's/(\s)([^\/\s]+.pb.h)/\1protoobj\/\2/g' \
	>$@

#These rules do the compilation
$(OBJDIR)/%.pb.o: protoobj/%.pb.cc protoobj/%.pb.h
	@mkdir -p $(@D)
	$(CXX) $(OPTIONS) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: src/%.cc $(OBJDIR)/%.d
	@mkdir -p $(dir $@)
	$(CXX) $(OPTIONS) $(CXXFLAGS) -o $@ -c $<

$(EXECUTABLES):
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(LDFLAGS)
