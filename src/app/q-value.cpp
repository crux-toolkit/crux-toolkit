/*******************************************************************************
  * \file q-value.cpp
  * AUTHOR: Chris Park
  * CREATE DATE: Jan 03 2007
  * \brief  Given as input a directory containing binary psm files,
  * a protein database, and an optional parameter file, analyze the
  * matches (with percolator or q-value) and return scores indicating
  * how good the matches are. 
  *
  * Handles at most 4 files (target and decoy).  Expects psm files to
  * start with <fileroot>.se and 
  * end with the extension '.txt' and decoys to end with
  * '-decoy#.txt'.  Multiple target files in the given directory are
  * concatinated together and presumed to be non-overlaping parts of
  * the same ms2 file. 
  ****************************************************************************/

#include "q-value.h"
#include "io/MatchCollectionParser.h"
#include "PosteriorEstimator.h"
#include "util/FileUtils.h"
#include "util/Params.h"

#include <map>
#include <utility>
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
