#ifndef TIDE_LITE_MATCH_SET_H
#define TIDE_LITE_MATCH_SET_H

#define  BOOST_DATE_TIME_NO_LIB
#include <boost/thread.hpp>
#include <vector>
#include "raw_proteins.pb.h"
#include "tide/records.h"
#include "tide/active_peptide_queue.h"  // no include guard
#include "tide/fixed_cap_array.h"
#include "tide/peptide.h"
#include "tide/sp_scorer.h"
#include "tide/spectrum_collection.h"

#include "model/Modification.h"
#include "model/PostProcessProtein.h"

using namespace std;

typedef vector<const pb::Protein*> ProteinVec;

class TideLiteMatchSet {
};

#endif
