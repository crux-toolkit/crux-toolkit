/**
 * \file create_index.h
 */
/*
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 * DESCRIPTION: Header file for the create-index command of crux
 * REVISION: $Revision: 1.4 $
 */

#ifndef CREATE_INDEX_H
#define CREATE_INDEX_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <signal.h>
#include "objects.h"
#include "carp.h"
#include "peptide.h"
#include "peptide_src.h"
#include "protein.h"
#include "database.h"
#include "index.h"
#include "protein_index.h"
#include "parameter.h"

#define NUM_INDEX_OPTIONS 14
#define NUM_INDEX_ARGS 2

int create_index_main(int argc, char** argv);





#endif //CREATE_INDEX_H











