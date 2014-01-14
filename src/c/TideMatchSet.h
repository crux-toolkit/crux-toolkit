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

class TideMatchSet {

public:
  typedef pair<int, int> Pair;
  typedef FixedCapacityArray<Pair> Arr;

  // Matches will be an array of pairs, (score, counter), where counter refers
  // to the index within the ActivePeptideQueue, counting from the back.  This
  // slight complication is due to the way the generated machine code fills the
  // counter in the matches buffer by decrementing the counter.
  TideMatchSet(
    Arr* matches,
    double max_mz
  );

  ~TideMatchSet();

  /**
   * Write matches to output files
   */
  void report(
    ofstream* target_file,  ///< target file to write to
    ofstream* decoy_file, ///< decoy file to write to
    int top_n,  ///< number of matches to report
    const Spectrum* spectrum, ///< spectrum for matches
    int charge, ///< charge for matches
    const ActivePeptideQueue* peptides, ///< peptide queue
    const ProteinVec& proteins, ///< proteins corresponding with peptides
    const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
    bool compute_sp ///< whether to compute sp or not
  );

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
    const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
    bool compute_sp ///< whether to compute sp or not
  );

  static void writeHeaders(
    ofstream* file,
    bool decoyFile,
    bool sp
  );

  static void initModMap(
    const pb::ModTable& modTable
  );

  static void setCleavageType(
    const string& cleavageType
  );

protected:
  Arr* matches_;
  double max_mz_;
  static map<int, double> mod_map_; // unique delta index -> delta
  static ModCoder mod_coder_;
  static string cleavage_type_;

  // For allocation
  static char match_collection_loc_[sizeof(MatchCollection)];
  static char decoy_match_collection_loc_[sizeof(MatchCollection)];

  struct less_score : public binary_function<Pair, Pair, bool> {
    // Compare scores, ignore counters.
    bool operator()(Pair x, Pair y) { return x.first < y.first; }
  };

  /**
   * Helper function for tab delimited report function
   */
  void writeToFile(
    ofstream* file,
    int top_n,
    const vector<Arr::iterator>& vec,
    const Spectrum* spectrum,
    int charge,
    const ActivePeptideQueue* peptides,
    const ProteinVec& proteins,
    const vector<const pb::AuxLocation*>& locations,
    const map<Arr::iterator, FLOAT_T>& delta_cn_map,
    const map<Arr::iterator, pair<const SpScorer::SpScoreData, int> >* sp_map
  );

  /**
   * Helper function for normal report function
   */
  void addCruxMatches(
    MatchCollection* match_collection,
    bool decoys,
    int top_n,
    vector<PostProcessProtein*>* proteins_made,
    const vector<Arr::iterator>& vec,
    Crux::Spectrum& spectrum,
    const ActivePeptideQueue* peptides,
    const ProteinVec& proteins,
    const vector<const pb::AuxLocation*>& locations,
    SpectrumZState& z_state,
    SpScorer* sp_scorer,
    FLOAT_T* lowest_sp_out
  );

  /**
   * Create a Crux match from Tide data structures
   */
  Crux::Match* getCruxMatch(
    const Peptide* peptide, ///< Tide peptide for match
    const ProteinVec& proteins, ///< Tide proteins
    const vector<const pb::AuxLocation*>& locations, /// auxiliary locations
    Crux::Spectrum* crux_spectrum,  ///< Crux spectrum for match
    SpectrumZState& crux_z_state, ///< Crux z state for match
    vector<PostProcessProtein*>* proteins_made ///< out parameter for new proteins
  );

  /**
   * Returns a pointer to the modification in the list of mods, adding it if it
   * doesn't exist
   */
  const AA_MOD_T* lookUpMod(
    double delta_mass ///< mass of the mod to look up
  );

  void gatherTargetsAndDecoys(
    const ActivePeptideQueue* peptides,
    const ProteinVec& proteins,
    vector<Arr::iterator>& targetsOut,
    vector<Arr::iterator>& decoysOut,
    int top_n
  );

  /**
   * Create a pb peptide from Tide peptide
   */
  static pb::Peptide* getPbPeptide(
    const Peptide& peptide
  );

  /**
   * Gets the protein name with the index appended.
   */
  static string getProteinName(
    const pb::Protein& protein,
    int pos
  );

  /**
   * Gets the flanking AAs for a Tide peptide sequence
   */
  static void getFlankingAAs(
    const Peptide* peptide, ///< Tide peptide to get flanking AAs for
    const pb::Protein* protein, ///< Tide protein for the peptide
    int pos,  ///< location of peptide within protein
    string* out_n,  ///< out parameter for n flank
    string* out_c ///< out parameter for c flank
  );

  static void computeDeltaCns(
    const vector<Arr::iterator>& vec, // xcorr*100000000.0, high to low
    map<Arr::iterator, FLOAT_T>* delta_cn_map, // map to add delta cn scores to
    int top_n // number of top matches we will be reporting
  );

  static void computeSpData(
    const vector<Arr::iterator>& vec,
    map<Arr::iterator, pair<const SpScorer::SpScoreData, int> >* sp_rank_map,
    SpScorer* sp_scorer,
    const ActivePeptideQueue* peptides
  );

  struct spGreater {
    inline bool operator() (const pair<Arr::iterator, SpScorer::SpScoreData>& lhs,
                            const pair<Arr::iterator, SpScorer::SpScoreData>& rhs) {
      return lhs.second.sp_score > rhs.second.sp_score;
    }
  };
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
