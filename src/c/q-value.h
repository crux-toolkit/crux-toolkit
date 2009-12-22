#ifndef QVALUE_CMD_H
#define QVALUE_CMD_H

/**
 * \file match_analysis.c
 */
/*
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 * DESCRIPTION: Header file for compute-q-values crux command
 * $Revision: 1.2 $
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "protein.h"
#include "peptide.h"
#include "spectrum.h"
#include "parse_arguments.h" 
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "scorer.h"
#include "match.h"
#include "match_collection.h"
#include "output-files.h"

int qvalue_main(int argc, char** argv);

#endif //QVALUE_CMD_H

