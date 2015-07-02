#include <gflags/gflags.h>
#include "header.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "spectrum_preprocess.h"
#include "max_mz.h"
#include "modifications.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_string(seq, "", "Sequence of amino acids");
DEFINE_bool(make_all, false, "Use TheoreticalPeakSetMakeAll");
DEFINE_bool(by_all, false, "Use TheoreticalPeakSetBYAll");
DEFINE_bool(by_sparse, false, "Use TheoreticalPeakSetBYSparse");
DEFINE_bool(sparse, false, "Use TheoreticalPeakSetSparse");
DEFINE_bool(diff, false, "Use TheoreticalPeakSetDiff");
DEFINE_string(mods, "C+57", "mods specification");

typedef vector<const pb::Protein*> ProteinVec;


class ObservedPeakTester {
 public:
  static void Test(TheoreticalPeakSet* workspace) {
    ObservedPeakSet obs;
    
    obs.max_mz_ = MaxMZ::Global();
    obs.cache_end_ = obs.max_mz_.CacheBinEnd() * NUM_PEAK_TYPES;
    obs.cache_ = new int[obs.cache_end_];
    memset(obs.cache_, 0, obs.cache_end_ * sizeof(*obs.cache_));

    for (int charge = 2; charge <= 3; ++charge) {
      cout << "Charge " << charge << endl;
      for (int i = 0; i < obs.max_mz_.BackgroundBinEnd(); ++i) {
	//      cout << i << endl;
	obs.Peak(PeakMain, i) = 1;
	if (i > 0)
	  obs.Peak(PeakMain, i-1) = 0;
	//obs.ShowAll();
	obs.ComputeCache();
	// obs.ShowAll();
	
	int dot = DotProd(workspace, &obs, charge);
	if (dot != 0)
	  cout << "[" << i << "]: " << dot << endl;
      }
    }
  }

 private:
  static int DotProd(TheoreticalPeakSet* workspace, ObservedPeakSet* observed, int charge) {
    TheoreticalPeakArr peaks_charge_1(2000);
    TheoreticalPeakArr peaks_charge_2(2000);
    TheoreticalPeakArr negs_charge_1(2000);
    TheoreticalPeakArr negs_charge_2(2000);
    workspace->GetPeaks(&peaks_charge_1, &negs_charge_1,
			&peaks_charge_2, &negs_charge_2, NULL);
    TheoreticalPeakArr* peaks = (charge >= 3) ? &peaks_charge_2 : &peaks_charge_1;
    TheoreticalPeakArr* negs = (charge >= 3) ? &negs_charge_2 : &negs_charge_1;
    return observed->DotProd(*peaks) - observed->DotProd(*negs);
  }
};

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  FLAGS_max_mz = 2500; // new default

  google::ParseCommandLineFlags(&argc, &argv, true);

  VariableModTable var_mod_table;
  pb::ModTable mod_table;
  CHECK(var_mod_table.Parse(FLAGS_mods.c_str()));
  mod_table.CopyFrom(*var_mod_table.ParsedModTable());
  CHECK(MassConstants::Init(&mod_table));

  MaxMZ::SetGlobalMaxFromFlag();

  CHECK(!FLAGS_seq.empty());
  CHECK(int(FLAGS_make_all) + int(FLAGS_by_all) + int(FLAGS_by_sparse) + int(FLAGS_sparse) + int(FLAGS_diff) == 1);

  // make one "protein" out of FLAGS_seq
  pb::Protein protein;
  protein.set_residues(FLAGS_seq);
  ProteinVec protein_vec;
  protein_vec.push_back(&protein);
  
  // make one peptide be the whole protein
  pb::Peptide pb_peptide;
  pb_peptide.set_length(FLAGS_seq.length());
  pb_peptide.mutable_first_location()->set_protein_id(0);
  pb_peptide.mutable_first_location()->set_pos(0);

  TheoreticalPeakSet* workspace;
  if (FLAGS_make_all)
    workspace = new TheoreticalPeakSetMakeAll(2000);
  else if (FLAGS_by_all)
    workspace = new TheoreticalPeakSetBYAll(2000);
  else if (FLAGS_by_sparse)
    workspace = new TheoreticalPeakSetBYSparse(2000);
  else if (FLAGS_sparse)
    workspace = new TheoreticalPeakSetSparse(2000);
  else if (FLAGS_diff)
    workspace = new TheoreticalPeakSetDiff(2000);

  Peptide peptide(pb_peptide, protein_vec);
  peptide.ComputeTheoreticalPeaks(workspace);

  ObservedPeakTester::Test(workspace);

  delete workspace;

  return 0;
}
