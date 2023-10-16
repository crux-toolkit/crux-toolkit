#include <iomanip>

#include "TideIndexApplication.h"
#include "TideLiteMatchSet.h"
#include "TideSearchApplication.h"
#include "util/Params.h"
#include "util/StringUtils.h"

string TideMatchSet::CleavageType;
string TideMatchSet::decoy_prefix_;

char TideMatchSet::match_collection_loc_[] = {0};
char TideMatchSet::decoy_match_collection_loc_[] = {0};

TideMatchSet::TideMatchSet(Arr* matches, double max_mz)
  : matches_(matches), max_mz_(max_mz), exact_pval_search_(false), elution_window_(0), cur_score_function_(XCORR_SCORE) {
}

TideMatchSet::TideMatchSet(Peptide* peptide, double max_mz)
  : peptide_(peptide), max_mz_(max_mz), exact_pval_search_(false), elution_window_(0), cur_score_function_(XCORR_SCORE) {
}

TideMatchSet::~TideMatchSet() {
}


