#ifndef PERCOLATOR_CMD_H
#define PERCOLATOR_CMD_H

/**
 * \file percolator.h
 */
/*
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 25, 2008
 * DESCRIPTION: Header file for crux percolator command. To be made
 * obsolte by q-ranker
 *
 * $Revision: 1.2 $
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
#include "PercolatorCInterface.h"
#include "output-files.h"

int percolator_main(int argc, char** argv);

#endif //PERCOLATOR_CMD_H

