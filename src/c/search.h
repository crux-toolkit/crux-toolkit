/**
 * \file match_search.c
 */
/* AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 * DESCRIPTION: Header file for the search-for-matches crux command
 * REVISION: $Revision: 1.3 $
 */
#ifndef SEARCH_H
#define SEARCH_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "carp.h"
#include "peptide.h"
#include "protein.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "crux-utils.h"
#include "scorer.h"
#include "objects.h"
#include "match.h"
#include "match_collection.h"

#define NUM_SEARCH_OPTIONS 16
#define NUM_SEARCH_ARGS 2

/* Private functions */
int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database);
void open_output_files(FILE*** binary_filehandle_array, 
                       FILE** sqt_filehandle,
                       FILE** decoy_sqt_filehandle);

int search_main(int argc, char** argv);

#ifdef __cplusplus
}
#endif

#endif // SEARCH_H
