#include <fstream>
#include <iomanip>

#include "TideIndexApplication.h"
#include "TideMatchSet.h"
#include "TideSearchApplication.h"

extern AA_MOD_T* list_of_mods[MAX_AA_MODS]; // list containing all aa mods
extern int num_mods;  // ANY_POSITION mods

map<int, double> TideMatchSet::mod_map_;
ModCoder TideMatchSet::mod_coder_;
string TideMatchSet::cleavage_type_ = "";
char TideMatchSet::match_collection_loc_[] = {0};
char TideMatchSet::decoy_match_collection_loc_[] = {0};

TideMatchSet::TideMatchSet(
  Arr* matches,
  double max_mz
) :
  matches_(matches), max_mz_(max_mz) {
}

TideMatchSet::~TideMatchSet() {
}

/**
 * Write matches to output files
 * This is for writing tab-delimited only
 */
void TideMatchSet::report(
  ofstream* target_file,  ///< target file to write to
  ofstream* decoy_file, ///< decoy file to write to
  int top_n,  ///< number of matches to report
  const Spectrum* spectrum, ///< spectrum for matches
  int charge, ///< charge for matches
  const ActivePeptideQueue* peptides, ///< peptide queue
  const ProteinVec& proteins,  ///< proteins corresponding with peptides
  const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
  bool compute_sp ///< whether to compute sp or not
) {
  if (matches_->size() == 0) {
    return;
  }

  carp(CARP_DETAILED_DEBUG, "Tide MatchSet reporting top %d of %d matches",
       top_n, matches_->size());

  vector<Arr::iterator> targets, decoys;
  gatherTargetsAndDecoys(peptides, proteins, targets, decoys, top_n);

  map<Arr::iterator, FLOAT_T> delta_cn_map;
  computeDeltaCns(targets, &delta_cn_map, top_n);
  computeDeltaCns(decoys, &delta_cn_map, top_n);

  map<Arr::iterator, pair<const SpScorer::SpScoreData, int> > sp_map;
  if (compute_sp) {
    SpScorer sp_scorer(proteins, *spectrum, charge, max_mz_);
    computeSpData(targets, &sp_map, &sp_scorer, peptides);
    computeSpData(decoys, &sp_map, &sp_scorer, peptides);
  }

  writeToFile(target_file, top_n, targets, spectrum, charge, peptides, proteins,
              locations, delta_cn_map, compute_sp ? &sp_map : NULL);
  writeToFile(decoy_file, top_n, decoys, spectrum, charge, peptides, proteins,
              locations, delta_cn_map, compute_sp ? &sp_map : NULL);
}

/**
 * Helper function for tab delimited report function
 */
void TideMatchSet::writeToFile(
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
) {
  if (!file) {
    return;
  }

  int cur = 0;
  const vector<Arr::iterator>::const_iterator cutoff =
    (vec.size() >= top_n) ? vec.begin() + top_n : vec.end();
  for (vector<Arr::iterator>::const_iterator i = vec.begin(); i != cutoff; ++i) {
    const Peptide* peptide = peptides->GetPeptide((*i)->second);
    const pb::Protein* protein = proteins[peptide->FirstLocProteinId()];
    int pos = peptide->FirstLocPos();
    string proteinNames = getProteinName(*protein,
      (!protein->has_target_pos()) ? pos : protein->target_pos());
    string flankingAAs, n_term, c_term;
    getFlankingAAs(peptide, protein, pos, &n_term, &c_term);
    flankingAAs = n_term + c_term;

    // look for other locations
    if (peptide->HasAuxLocationsIndex()) {
      const pb::AuxLocation* aux = locations[peptide->AuxLocationsIndex()];
      for (int i = 0; i < aux->location_size(); ++i) {
        const pb::Location& location = aux->location(i);
        protein = proteins[location.protein_id()];
        pos = location.pos();
        proteinNames += "," + getProteinName(*protein,
          (!protein->has_target_pos()) ? pos : protein->target_pos());
        getFlankingAAs(peptide, protein, pos, &n_term, &c_term);
        flankingAAs += "," + n_term + c_term;
      }
    }

    map<size_t, double> modMap; // AA index -> mod delta
    const ModCoder::Mod* mods;
    int pep_mods = peptide->Mods(&mods);
    for (int j = 0; j < pep_mods; ++j) {
      int mod_index, mod_delta_index;
      mod_coder_.DecodeMod(mods[j], &mod_index, &mod_delta_index);
      double mod_delta = mod_map_[mod_delta_index];

      map<size_t, double>::iterator lookup = modMap.find(mod_index);
      if (lookup == modMap.end()) {
        modMap[mod_index] = mod_delta;
      } else {
        modMap[mod_index] += mod_delta;
      }
    }
    string seq = peptide->Seq();
    for (size_t j = seq.length() - 1; j >= 0 && !modMap.empty(); --j) {
      map<size_t, double>::iterator lookup = modMap.find(j);
      if (lookup == modMap.end()) {
        continue;
      }
      stringstream ss;
      ss << '[' << lookup->second << ']';
      modMap.erase(lookup);
      seq.insert(j + 1, ss.str());
    }

    const SpScorer::SpScoreData* sp_data = sp_map ? &(sp_map->at(*i).first) : NULL;

    *file << spectrum->SpectrumNumber() << '\t'
          << charge << '\t'
          << spectrum->PrecursorMZ() << '\t'
          << (spectrum->PrecursorMZ() - MASS_PROTON) * charge << '\t'
          << peptide->Mass() << '\t'
          << delta_cn_map.at(*i) << '\t';
    if (sp_map) {
      *file << sp_data->sp_score << '\t'
            << sp_map->at(*i).second << '\t';
    }
    *file << ((*i)->first / 100000000.0) << '\t'
          << ++cur << '\t';
    if (sp_map) {
      *file << sp_data->matched_ions << '\t'
            << sp_data->total_ions << '\t';
    }
    *file << (!peptide->IsDecoy() ? peptides->ActiveTargets() : peptides->ActiveDecoys()) << '\t'
          << seq << '\t'
          << cleavage_type_ << '\t'
          << proteinNames << '\t'
          << flankingAAs;
    if (peptide->IsDecoy() && !OutputFiles::isProteinLevelDecoys()) {
      // write target sequence
      const string& residues = protein->residues();
      *file << '\t'
            << residues.substr(residues.length() - peptide->Len());
    } else if (OutputFiles::isConcat() && !OutputFiles::isProteinLevelDecoys()) {
      *file << '\t'
            << seq;
    }
    *file << endl;
  }
}

/**
 * Write matches to output files
 */
void TideMatchSet::report(
  OutputFiles* output_files,  ///< pointer to output handler
  int top_n,  ///< number of matches to report
  const Spectrum* spectrum, ///< spectrum for matches
  int charge, ///< charge for matches
  const ActivePeptideQueue* peptides, ///< peptide queue
  const ProteinVec& proteins,  ///< proteins corresponding with peptides
  const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
  bool compute_sp ///< whether to compute sp or not
) {
  if (matches_->size() == 0) {
    return;
  }

  carp(CARP_DETAILED_DEBUG, "Tide MatchSet reporting top %d of %d matches",
       top_n, matches_->size());

  vector<Arr::iterator> targets, decoys;
  gatherTargetsAndDecoys(peptides, proteins, targets, decoys, top_n);

  MatchCollection* crux_collection =
    new(match_collection_loc_) MatchCollection();
  MatchCollection* crux_decoy_collection =
    new(decoy_match_collection_loc_) MatchCollection();
  vector<PostProcessProtein*> proteins_made;

  // For Sp scoring
  FLOAT_T lowest_sp = BILLION;
  SpScorer* sp_scorer = (compute_sp) ?
    new SpScorer(proteins, *spectrum, charge, max_mz_) : NULL;

  // Create a Crux spectrum and z state
  Crux::Spectrum crux_spectrum(
    spectrum->SpectrumNumber(), spectrum->SpectrumNumber(),
    spectrum->PrecursorMZ(), vector<int>(1, charge), "");
  SpectrumZState z_state;
  z_state.setMZ(crux_spectrum.getPrecursorMz(), charge);

  addCruxMatches(crux_collection, false, top_n, &proteins_made, targets, crux_spectrum,
                 peptides, proteins, locations, z_state, sp_scorer, &lowest_sp);
  addCruxMatches(crux_decoy_collection, true, top_n, &proteins_made, decoys, crux_spectrum,
                 peptides, proteins, locations, z_state, sp_scorer, &lowest_sp);

  if (sp_scorer) {
    crux_spectrum.setTotalEnergy(sp_scorer->TotalIonIntensity());
    crux_spectrum.setLowestSp(lowest_sp);
    delete sp_scorer;
  }

  // Write matches
  vector<MatchCollection*> decoy_vector;
  if (!OutputFiles::isConcat()) {
    decoy_vector.push_back(crux_decoy_collection);
  }
  output_files->writeMatches(crux_collection, decoy_vector, XCORR, &crux_spectrum);

  // Clean up
  crux_collection->~MatchCollection();
  crux_decoy_collection->~MatchCollection();
  for (vector<PostProcessProtein*>::iterator i = proteins_made.begin();
       i != proteins_made.end();
       ++i) {
    delete *i;
  }
}

/**
 * Helper function for normal report function
 */
void TideMatchSet::addCruxMatches(
  MatchCollection* match_collection,
  bool decoys,
  int top_n,
  vector<PostProcessProtein*>* proteins_made,
  const vector<Arr::iterator>& vec,
  Crux::Spectrum& crux_spectrum,
  const ActivePeptideQueue* peptides,
  const ProteinVec& proteins,
  const vector<const pb::AuxLocation*>& locations,
  SpectrumZState& z_state,
  SpScorer* sp_scorer,
  FLOAT_T* lowest_sp_out
) {

  FLOAT_T lnNumSp = log(!decoys ? peptides->ActiveTargets()
                                : peptides->ActiveDecoys());
  
  // Create a Crux match for each match
  vector<Arr::iterator>::const_iterator endIter = min(vec.begin() + top_n, vec.end());
  for (vector<Arr::iterator>::const_iterator i = vec.begin(); i != endIter; ++i) {
    const Peptide* peptide = peptides->GetPeptide((*i)->second);

    Crux::Match* match = getCruxMatch(peptide, proteins, locations, &crux_spectrum,
                                      z_state, proteins_made);
    match_collection->addMatch(match);
    if (decoys) {
      match->setNullPeptide(true);
    }
    Crux::Match::freeMatch(match); // so match gets deleted when collection does

    // Set Xcorr score in match
    match->setScore(XCORR, (*i)->first / 100000000.0);
    // Set lnNumSp in match
    match->setLnExperimentSize(lnNumSp);

    if (sp_scorer) {
      pb::Peptide* pb_peptide = getPbPeptide(*peptide);

      // Score for Sp
      SpScorer::SpScoreData sp_score_data;
      sp_scorer->Score(*pb_peptide, sp_score_data);
      delete pb_peptide;

      FLOAT_T sp = sp_score_data.sp_score;
      if (sp < *lowest_sp_out) {
        *lowest_sp_out = sp;
      }

      // Set Sp, B/Y scores in match
      match->setScore(SP, sp);
      match->setBYIonMatched(sp_score_data.matched_ions);
      match->setBYIonPossible(sp_score_data.total_ions);
    }
  }
  match_collection->setZState(z_state);
  match_collection->setExperimentSize(!decoys ? peptides->ActiveTargets() : peptides->ActiveDecoys());
  match_collection->populateMatchRank(XCORR);
  match_collection->forceScoredBy(XCORR);
  if (sp_scorer) {
    match_collection->populateMatchRank(SP);
    match_collection->forceScoredBy(SP);
  }
}

/**
 * Write headers for tab delimited file
 */
void TideMatchSet::writeHeaders(
  ofstream* file,
  bool decoyFile,
  bool sp
) {
  if (!file) {
    return;
  }
  const int headers[] = {
    SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, SP_SCORE_COL, SP_RANK_COL, XCORR_SCORE_COL,
    XCORR_RANK_COL, BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL,
    MATCHES_SPECTRUM_COL, SEQUENCE_COL, CLEAVAGE_TYPE_COL, PROTEIN_ID_COL,
    FLANKING_AA_COL, ORIGINAL_TARGET_SEQUENCE_COL
  };
  size_t numHeaders = sizeof(headers) / sizeof(int);
  for (size_t i = 0; i < numHeaders; ++i) {
    int header = headers[i];
    if (!sp &&
        (header == SP_SCORE_COL || header == SP_RANK_COL ||
         header == BY_IONS_MATCHED_COL || header == BY_IONS_TOTAL_COL)) {
      continue;
    } else if (header == ORIGINAL_TARGET_SEQUENCE_COL &&
               (OutputFiles::isProteinLevelDecoys() ||
                (!decoyFile && !OutputFiles::isConcat()))) {
      continue;
    }
    if (i > 0) {
      *file << '\t';
    }
    *file << get_column_header(header);
  }
  *file << endl;
}

void TideMatchSet::initModMap(
  const pb::ModTable& modTable
) {
  mod_map_.clear();
  int numDeltas = modTable.unique_deltas_size();
  for (int i = 0; i < numDeltas; ++i) {
    mod_map_[i] = modTable.unique_deltas(i);
  }
  mod_coder_.Init(numDeltas);
}

void TideMatchSet::setCleavageType(
  const string& cleavageType
) {
  cleavage_type_ = cleavageType;
}

/**
 * Create a Crux match from Tide data structures
 */
Crux::Match* TideMatchSet::getCruxMatch(
  const Peptide* peptide, ///< Tide peptide for match
  const ProteinVec& proteins, ///< Tide proteins
  const vector<const pb::AuxLocation*>& locations, /// auxiliary locations
  Crux::Spectrum* crux_spectrum,  ///< Crux spectrum for match
  SpectrumZState& crux_z_state, ///< Crux z state for match
  vector<PostProcessProtein*>* proteins_made ///< out parameter for new proteins
) {
  const pb::Protein* protein = proteins[peptide->FirstLocProteinId()];
  int pos = peptide->FirstLocPos();

  // Get flanking AAs
  string n_term, c_term;
  getFlankingAAs(peptide, protein, pos, &n_term, &c_term);

  // Create protein
  PostProcessProtein* crux_protein = new PostProcessProtein();
  proteins_made->push_back(crux_protein);

  crux_protein->setId(protein->name().c_str());
  string originalTargetSeq = !peptide->IsDecoy() ? peptide->Seq() :
    protein->residues().substr(protein->residues().length() - peptide->Len() - 1);
  int start_idx = crux_protein->findStart(originalTargetSeq, n_term, c_term);

  // Create peptide
  Crux::Peptide* crux_peptide = new Crux::Peptide(
    peptide->Len(), peptide->Mass(), crux_protein, start_idx);
  crux_peptide->getPeptideSrc()->setStartIdxOriginal(pos + 1);

  // Add other proteins if any
  if (peptide->HasAuxLocationsIndex()) {
    const pb::AuxLocation* aux = locations[peptide->AuxLocationsIndex()];
    for (int i = 0; i < aux->location_size(); ++i) {
      const pb::Location& location = aux->location(i);
      protein = proteins[location.protein_id()];
      pos = location.pos();
      getFlankingAAs(peptide, protein, pos, &n_term, &c_term);

      crux_protein = new PostProcessProtein();
      proteins_made->push_back(crux_protein);
      crux_protein->setId(protein->name().c_str());
      start_idx = crux_protein->findStart(originalTargetSeq, n_term, c_term);
      PeptideSrc* src = new PeptideSrc(NON_SPECIFIC_DIGEST, crux_protein, start_idx);
      src->setStartIdxOriginal(pos + 1);
      crux_peptide->addPeptideSrc(src);
    }
  }

  // Set up modifications for peptide
  const ModCoder::Mod* mods;
  int pep_mods = peptide->Mods(&mods);
  MODIFIED_AA_T* mod_seq;
  convert_to_mod_aa_seq(peptide->Seq().c_str(), &mod_seq);
  for (int i = 0; i < pep_mods; ++i) {
    int mod_index, mod_delta_index;
    mod_coder_.DecodeMod(mods[i], &mod_index, &mod_delta_index);
    double mod_delta = mod_map_[mod_delta_index];
    // Look up mod and apply it to AA
    const AA_MOD_T* mod = lookUpMod(mod_delta);
    modify_aa(mod_seq + mod_index, mod);
  }
  crux_peptide->setModifiedAASequence(mod_seq, peptide->IsDecoy());
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
const AA_MOD_T* TideMatchSet::lookUpMod(double delta_mass) {
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

void TideMatchSet::gatherTargetsAndDecoys(
  const ActivePeptideQueue* peptides,
  const ProteinVec& proteins,
  vector<Arr::iterator>& targetsOut,
  vector<Arr::iterator>& decoysOut,
  int top_n
) {
  make_heap(matches_->begin(), matches_->end(), less_score());
  if (!OutputFiles::isConcat() && TideSearchApplication::hasDecoys()) {
    for (Arr::iterator i = matches_->end(); i != matches_->begin(); ) {
      pop_heap(matches_->begin(), i--, less_score());
      const Peptide& peptide = *(peptides->GetPeptide(i->second));
      const pb::Protein& protein = *(proteins[peptide.FirstLocProteinId()]);
      vector<Arr::iterator>* vec_ptr = !peptide.IsDecoy() ? &targetsOut : &decoysOut;
      if (vec_ptr->size() < top_n + 1) {
        vec_ptr->push_back(i);
      }
    }
  } else {
    int toAdd = min(top_n + 1, matches_->size());
    for (int i = 0; i < toAdd; ) {
      pop_heap(matches_->begin(), matches_->end() - i, less_score());
      targetsOut.push_back(matches_->end() - (++i));
    }
  }
}

/**
 * Create a pb peptide from Tide peptide
 */
pb::Peptide* TideMatchSet::getPbPeptide(
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
 * Gets the protein name with the index appended.
 */
string TideMatchSet::getProteinName(
  const pb::Protein& protein,
  int pos
) {
  stringstream proteinNameStream;
  proteinNameStream << protein.name()
                    << '(' << pos + 1 << ')';
  return proteinNameStream.str();
}

/**
 * Gets the flanking AAs for a Tide peptide sequence
 */
void TideMatchSet::getFlankingAAs(
  const Peptide* peptide, ///< Tide peptide to get flanking AAs for
  const pb::Protein* protein, ///< Tide protein for the peptide
  int pos,  ///< location of peptide within protein
  string* out_n,  ///< out parameter for n flank
  string* out_c ///< out parameter for c flank
) {
  int idx_n = pos - 1;
  int idx_c = pos + peptide->Len();
  const string& seq = protein->residues();

  *out_n = (idx_n >= 0) ? seq.substr(idx_n, 1) : "-";
  *out_c = (idx_c < seq.length()) ? seq.substr(idx_c, 1) : "-";
}

void TideMatchSet::computeDeltaCns(
  const vector<Arr::iterator>& vec, // xcorr*100000000.0, high to low
  map<Arr::iterator, FLOAT_T>* delta_cn_map, // map to add delta cn scores to
  int top_n // number of top matches we will be reporting
) {
  FLOAT_T lastXcorr = BILLION;
  vector<Arr::iterator>::const_reverse_iterator i = (vec.size() > top_n) ?
    vec.rend() - (top_n + 1) : vec.rbegin();
  for (; i != vec.rend(); ++i) {
    const FLOAT_T xcorr = (*i)->first / 100000000.0;
    delta_cn_map->insert(make_pair(*i, (lastXcorr == BILLION) ?
      0 : (xcorr - lastXcorr) / max(xcorr, FLOAT_T(1))));
    lastXcorr = xcorr;
  }
}

void TideMatchSet::computeSpData(
  const vector<Arr::iterator>& vec,
  map<Arr::iterator, pair<const SpScorer::SpScoreData, int> >* sp_rank_map,
  SpScorer* sp_scorer,
  const ActivePeptideQueue* peptides
) {
  vector< pair<Arr::iterator, SpScorer::SpScoreData> > spData;
  spData.reserve(vec.size());
  for (vector<Arr::iterator>::const_iterator i = vec.begin(); i != vec.end(); ++i) {
    spData.push_back(make_pair(*i, SpScorer::SpScoreData()));
    const Peptide& peptide = *(peptides->GetPeptide((*i)->second));
    pb::Peptide* pb_peptide = getPbPeptide(peptide);
    sp_scorer->Score(*pb_peptide, spData.back().second);
    delete pb_peptide;
  }
  sort(spData.begin(), spData.end(), spGreater());
  for (size_t i = 0; i < spData.size(); ++i) {
    sp_rank_map->insert(make_pair(
      spData[i].first, make_pair(spData[i].second, i + 1)));
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
