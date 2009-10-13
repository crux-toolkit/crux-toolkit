/**
 * \file crux-main.h
 */
/*
 AUTHOR: Barbara Frewen
 CREATE DATE: November 24, 2008
 DESCRIPTION: The starting point for what were previously three
 separate protgrams--create-index, search-for-matches, analyze-matches.
 REVISION: $Revision: 1.2 $
*/

#include "carp.h"
#include "utils.h"
#include "crux-utils.h"
#include "create_index.h"
#include "search.h"
#include "q-ranker.h"
#include "q-value.h"
#include "percolator.h"
#include "print-processed-spectra.h"


/* Private data types */
enum _command { INDEX_CMD,      // create-index
                SEARCH_CMD,     // search-for-matches
                SEQUEST_CMD,     // sequest-search
                QVALUE_CMD,     // compute-q-values
                QRANKER_CMD,    // q-rakner
                PERCOLATOR_CMD, // percolator
                PROCESS_SPEC_CMD,// print-processed-spectra
                INVALID_CMD,    // use for errors

                NUM_CMD_TYPES}; // always keep this last
typedef enum _command COMMAND_T;

/* Private functions */
COMMAND_T string_to_command_type(char*);
char*     command_type_to_string(COMMAND_T);
