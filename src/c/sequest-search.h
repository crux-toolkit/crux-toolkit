#ifndef SEQUEST_SEARCH_H
#define SEQUEST_SEARCH_H

/**
 * \file sequest-search.h
 */
/*
 * FILE: sequest-search.h
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Oct 2, 2009
 * PROJECT: crux
 * DESCRIPTION: The crux search routine that emulates SEQUEST.  Scores
 * all candidate peptides with Sp, deletes all but the 500 top-scoring
 * candidates, scores remaining 500 with xcorr, sorts results by xcorr
 * and returns the top 5 plus the match with the best Sp score.
 * Writes results to .sqt and .txt.  Does not compute p-values. 
 */

#include <vector>
#include "OutputFiles.h"
#include "SearchProgress.h"
#include "objects.h"
#include "carp.h"
#include "parameter.h"
#include "Protein.h"
#include "peptide.h"
#include "SpectrumCollection.h"
#include "FilteredSpectrumChargeIterator.h"

using namespace std;

int sequest_search_main(int argc, char** argv);










#endif //SEQUEST_SEARCH_H
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
