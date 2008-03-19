#####################################################################
## FILE: Makefile.am
# AUTHOR: Barbara Frewen, Charles E. Grant
# CREATE DATE: 2008
# PROJECT: crux
# COPYRIGHT: 2008, UW
#####################################################################

## Automake uses this file to generate a Makefile

# This Makefile builds the scripts used to test the crux components.

## Call additional Makefiles in child directories
# BF: no subdirs so far
#SUBDIRS = clustalw2fasta draw-mhmm input mhmm mhmms scaffold

check_PROGRAMS = unittests
unittests_SOURCES = \
        ../utils.c \
        ../peak.c \
        ../spectrum.c \
        ../peptide.c \
        ../protein.c \
        ../peptide_src.c \
        ../mass.c \
        ../alphabet.c \
        ../array.c \
        ../crux-utils.c \
        ../database.c \
        ../carp.c \
        ../peptide_constraint.c \
         ../index.c \
        ../sorter.c \
        ../spectrum_collection.c \
        ../protein_index.c \
        ../ion.c \
        ../ion_series.c \
        ../parameter.c \
        ../scorer.c \
        ../parse_arguments.c \
        ../generate_peptides_iterator.c \
        ../match.c \
        ../match_collection.c \
        ../hash.c

#test_evo_models_LDADD = ../src/libmhmm.a -lm
# BF: there used to be ../lib/libcrux.a, now all libs are *.la
#     can we use one instead of all the *c files?
unittests_LDADD =  -lcheck -lm 

# Set up the Meta-MEME testing script.
TESTS = unittests
# none of these?
#check_SCRIPTS = metatest
TESTS_ENVIRONMENT = testdir=$(top_srcdir)/src/c/test

# Update the testing script with the local path to Perl.
#metatest: metatest.pl
#	sed -e 's#!/usr/bin/perl#!$(PERL_SHEBANG)#' $(srcdir)/metatest.pl > metatest
#	chmod ugo+x metatest

# These files should be removed as part of the 'make clean' target.
CLEANFILES = $(check_PROGRAMS) $(check_SCRIPTS)

# These files should be included with the distribution.
#EXTRA_DIST = <insert data files here?>
