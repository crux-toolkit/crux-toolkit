#ifndef QVALUE_H
#define QVALUE_H

/**
 * \file q-value.h
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "io/carp.h"
#include "util/crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "model/Protein.h"
#include "model/Peptide.h"
#include "model/Spectrum.h"
#include "util/parse_arguments.h" 
#include "io/SpectrumCollection.h"
#include "model/Scorer.h"
#include "model/Match.h"
#include "model/MatchCollection.h"
#include "io/OutputFiles.h"



FLOAT_T* compute_decoy_qvalues(
  FLOAT_T* target_scores,
  int      num_targets,
  FLOAT_T* decoy_scores,
  int      num_decoys,
  bool     reverse,
  FLOAT_T  pi_zero);

FLOAT_T* compute_qvalues_from_pvalues(
  FLOAT_T* pvalues, 
  int      num_pvals,
  FLOAT_T  pi_zero);

MatchCollection* run_qvalue(
  const char* psm_result_folder, 
  const char* fasta_file,
  OutputFiles& output );

#endif //QVALUE_H

