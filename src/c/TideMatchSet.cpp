#include <iomanip>

#include "TideMatchSet.h"

extern AA_MOD_T* list_of_mods[MAX_AA_MODS]; // list containing all aa mods
extern int num_mods;  // ANY_POSITION mods

MatchSet::MatchSet(
  Arr* matches,
  double max_mz
) {
  matches_ = matches;
  max_mz_ = max_mz;
}

MatchSet::~MatchSet() {
}

/**
 * Write matches to output files
 */
void MatchSet::report(
  OutputFiles* output_files,  ///< pointer to output handler
  int top_n,  ///< number of matches to report
  const Spectrum* spectrum, ///< spectrum for matches
  int charge, ///< charge for matches
  const ActivePeptideQueue* peptides, ///< peptide queue
  const ProteinVec& proteins,  ///< proteins corresponding with peptides
  bool compute_sp ///< whether to compute sp or not
) {
  if (matches_->size() == 0) {
    return;
  } else if (top_n > matches_->size()) {
    top_n = matches_->size();
  }
  carp(CARP_DETAILED_DEBUG, "Tide MatchSet reporting %d matches", top_n);
  getTop(top_n);

  MatchCollection* crux_collection = new MatchCollection();
  vector<PostProcessProtein*> proteins_made;

  // For Sp scoring
  FLOAT_T* lowest_sp = NULL;
  SpScorer* sp_scorer = (compute_sp) ?
    new SpScorer(proteins, *spectrum, charge, max_mz_) : NULL;

  // Create a Crux spectrum and z state
  Crux::Spectrum* crux_spectrum = new Crux::Spectrum(
    spectrum->SpectrumNumber(), spectrum->SpectrumNumber(),
    spectrum->PrecursorMZ(), vector<int>(1, charge), "");
  SpectrumZState z_state;
  z_state.setMZ(crux_spectrum->getPrecursorMz(), charge);

  // Create a Crux match for each match
  for (Arr::iterator i = matches_->end() - 1; top_n != 0; --i, --top_n) {
    // i->second is counter as it was during scoring and corresponds to the
    // index of the peptide in the ActivePeptideQueue, counting from the back.
    // GetPeptide() retrieves the corresponding Peptide.
    const Peptide* peptide = peptides->GetPeptide(i->second);
    const pb::Protein* protein = proteins[peptide->FirstLocProteinId()];

    PostProcessProtein* new_protein;
    Crux::Match* match = getCruxMatch(peptide, protein, crux_spectrum, z_state,
                                      &new_protein);
    proteins_made.push_back(new_protein);
    crux_collection->addMatch(match);
    Crux::Match::freeMatch(match); // so match gets deleted when collection does

    // Set Xcorr score in match
    match->setScore(XCORR, i->first / 100000000.0);

    if (compute_sp) {
      pb::Peptide* pb_peptide = getPbPeptide(*peptide);

      // Score for Sp
      SpScorer::SpScoreData sp_score_data;
      sp_scorer->Score(*pb_peptide, sp_score_data);
      delete pb_peptide;

      FLOAT_T sp = sp_score_data.sp_score;
      if (lowest_sp == NULL) {
        lowest_sp = new FLOAT_T(sp);
      } else if (sp < *lowest_sp) {
        *lowest_sp = sp;
      }

      // Set Sp, B/Y scores in match
      match->setScore(SP, sp);
      match->setBYIonMatched(sp_score_data.matched_ions);
      match->setBYIonPossible(sp_score_data.total_ions);
    }
  }

  // Set MatchCollection variables
  crux_collection->setZState(z_state);
  crux_collection->setExperimentSize(matches_->size());
  crux_collection->populateMatchRank(XCORR);
  crux_collection->forceScoredBy(XCORR);

  if (compute_sp) {
    crux_spectrum->setTotalEnergy(sp_scorer->TotalIonIntensity());
    if (lowest_sp) {
      crux_spectrum->setLowestSp(*lowest_sp);
      delete lowest_sp;
    }
    delete sp_scorer;

    crux_collection->populateMatchRank(SP);
    crux_collection->forceScoredBy(SP);
  }

  // Write matches
  vector<MatchCollection*> dummy;
  output_files->writeMatches(crux_collection, dummy, XCORR, crux_spectrum);

  // Clean up
  delete crux_collection;
  for (vector<PostProcessProtein*>::iterator i = proteins_made.begin();
       i != proteins_made.end();
       ++i) {
    delete *i;
  }
  delete crux_spectrum;
}

void MatchSet::getTop(
  int top_n
) {
  assert(top_n <= matches_->size());
  // move top n elements to end of array, with largest element last
  make_heap(matches_->begin(), matches_->end(), less_score());
  for (int i = 0; i < top_n; ++i) {
    pop_heap(matches_->begin(), matches_->end() - i, less_score());
  }
}

/**
 * Create a Crux match from Tide data structures
 */
Crux::Match* MatchSet::getCruxMatch(
  const Peptide* peptide, ///< Tide peptide for match
  const pb::Protein* protein, ///< Tide protein for match
  Crux::Spectrum* crux_spectrum,  ///< Crux spectrum for match
  SpectrumZState& crux_z_state, ///< Crux z state for match
  PostProcessProtein** protein_made ///< out parameter for new protein
) {
  // Get flanking AAs
  string n_term, c_term;
  getFlankingAAs(peptide, protein, &n_term, &c_term);

  // Create protein
  PostProcessProtein* parent_protein = new PostProcessProtein();
  *protein_made = parent_protein;
  parent_protein->setId(protein->name().c_str());
  string seq = peptide->Seq();
  int start_idx = parent_protein->findStart(seq, n_term, c_term);

  // Create peptide
  Crux::Peptide* crux_peptide = new Crux::Peptide(
    peptide->Len(), peptide->Mass(), parent_protein, start_idx);

  // Set up modifications for peptide
  const ModCoder::Mod* mods;
  int pep_mods = peptide->Mods(&mods);
  MODIFIED_AA_T* mod_seq;
  convert_to_mod_aa_seq(seq.c_str(), &mod_seq);
  for (int i = 0; i < pep_mods; ++i) {
    int mod_index; // 0 based
    double mod_delta;
    MassConstants::DecodeMod(mods[i], &mod_index, &mod_delta);
    // Look up mod and apply it to AA
    const AA_MOD_T* mod = lookUpMod(mod_delta);
    modify_aa(mod_seq + mod_index, mod);
  }
  bool is_decoy = false;
  crux_peptide->setModifiedAASequence(mod_seq, is_decoy);
  free(mod_seq);

  // Create match and return
  Crux::Match* match = new Crux::Match(
    crux_peptide, crux_spectrum, crux_z_state, false);

  return match;
}

/**
 * Returns a pointer to the modification in the list of mods, adding it if it
 * doesn't exist
 */
const AA_MOD_T* MatchSet::lookUpMod(double delta_mass) {
  for (int i = 0; i < num_mods; ++i) {
    const AA_MOD_T* mod = list_of_mods[i];
    if (aa_mod_get_mass_change(mod) == delta_mass) {
      carp(CARP_DETAILED_DEBUG, "Found existing mod (%f)", delta_mass);
      return mod;
    }
  }

  carp(CARP_DEBUG, "Adding new mod (%f)", delta_mass);
  AA_MOD_T* new_mod = new_aa_mod(num_mods);
  aa_mod_set_mass_change(new_mod, delta_mass);
  list_of_mods[num_mods] = new_mod;

  ++num_mods;

  return new_mod;
}

/**
 * Create a pb peptide from Tide peptide
 */
pb::Peptide* MatchSet::getPbPeptide(
  const Peptide& peptide
) {

  pb::Peptide* pb_peptide = new pb::Peptide();
  pb_peptide->set_id(peptide.Id());
  pb_peptide->set_mass(peptide.Mass());
  pb_peptide->set_length(peptide.Len());
  if (peptide.HasAuxLocationsIndex()) {
    pb_peptide->set_aux_locations_index(peptide.AuxLocationsIndex());
  }

  // Copy over all the modifications for this Peptide
  const ModCoder::Mod* mods;
  int pep_mods = peptide.Mods(&mods);
  for (int i = 0; i < pep_mods; ++i) {
    pb_peptide->add_modifications(mods[i]);
  }

  // Copy over the Peptide's first location within the first protein
  pb::Location* first_location = pb_peptide->mutable_first_location();
  first_location->set_protein_id(peptide.FirstLocProteinId());
  first_location->set_pos(peptide.FirstLocPos());

  return pb_peptide;
}

/**
 * Gets the flanking AAs for a Tide peptide sequence
 */
void MatchSet::getFlankingAAs(
  const Peptide* peptide, ///< Tide peptide to get flanking AAs for
  const pb::Protein* protein, ///< Tide protein for the peptide
  string* out_n,  ///< out parameter for n flank
  string* out_c ///< out parameter for c flank
) {
  int idx_n = peptide->FirstLocPos() - 1;
  int idx_c = peptide->FirstLocPos() + peptide->Len();
  const string& seq = protein->residues();

  *out_n = (idx_n >= 0) ?
    seq.substr(idx_n, 1) : "-";
  *out_c = (idx_c < seq.length()) ?
    seq.substr(idx_c, 1) : "-";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
