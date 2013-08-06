#ifndef TIDE_MATCH_SET_H
#define TIDE_MATCH_SET_H

#include <vector>
#include "raw_proteins.pb.h"
#include "tide/records.h"
#include "tide/active_peptide_queue.h"  // no include guard
#include "tide/fixed_cap_array.h"
#include "tide/peptide.h"
#include "tide/sp_scorer.h"
#include "tide/spectrum_collection.h"

#include "OutputFiles.h"
#include "PostProcessProtein.h"

using namespace std;

typedef vector<const pb::Protein*> ProteinVec;

class MatchSet {

public:
  typedef pair<int, int> Pair;
  typedef FixedCapacityArray<Pair> Arr;

  // Matches will be an array of pairs, (score, counter), where counter refers
  // to the index within the ActivePeptideQueue, counting from the back.  This
  // slight complication is due to the way the generated machine code fills the
  // counter in the matches buffer by decrementing the counter.
  MatchSet(
    Arr* matches,
    double max_mz
  );

  ~MatchSet();

  /**
   * Write matches to output files
   */
  void report(
    OutputFiles* output_files,  ///< pointer to output handler
    int top_n,  ///< number of matches to report
    const Spectrum* spectrum, ///< spectrum for matches
    int charge, ///< charge for matches
    const ActivePeptideQueue* peptides, ///< peptide queue
    const ProteinVec& proteins, ///< proteins corresponding with peptides
    bool compute_sp ///< whether to compute sp or not
  );

protected:
  Arr* matches_;
  double max_mz_;

  struct less_score : public binary_function<Pair, Pair, bool> {
    // Compare scores, ignore counters.
    bool operator()(Pair x, Pair y) { return x.first < y.first; }
  };

  /**
   * Create a Crux match from Tide data structures
   */
  Crux::Match* getCruxMatch(
    const Peptide* peptide, ///< Tide peptide for match
    const pb::Protein* protein, ///< Tide protein for match
    Crux::Spectrum* crux_spectrum,  ///< Crux spectrum for match
    SpectrumZState& crux_z_state, ///< Crux z state for match
    bool is_decoy,  /// Is the peptide a decoy
    PostProcessProtein** protein_made ///< out parameter for new protein
  );

  /**
   * Returns a poitner to the modification in the list of mods, adding it if it
   * doesn't exist
   */
  const AA_MOD_T* lookUpMod(
    double delta_mass ///< mass of the mod to look up
  );

  /**
   * Create a pb peptide from Tide peptide
   */
  pb::Peptide* getPbPeptide(
    const Peptide& peptide
  );

  /**
   * Gets the flanking AAs for a Tide peptide sequence
   */
  void getFlankingAAs(
    const Peptide* peptide, ///< Tide peptide to get flanking AAs for
    const pb::Protein* protein, ///< Tide protein for the peptide
    string* out_n,  ///< out parameter for n flank
    string* out_c ///< out parameter for c flank
  );

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
