/*
 * There are two versions of the report function, which writes matches to output
 * files. The first version, which takes ofstreams as arguments, is used when
 * only tab-delimited output is required. It does not perform any object
 * conversions. The second version takes an OutputFiles object as an argument
 * and is used when any non-tab-delimited output is required. It must convert
 * the data from Tide into Crux objects, which increases runtime.
 */

#include <fstream>
#include <iomanip>

#include "TideIndexApplication.h"
#include "TideMatchSet.h"
#include "TideSearchApplication.h"
#include "util/Params.h"
#include "util/StringUtils.h"

string TideMatchSet::CleavageType;
char TideMatchSet::match_collection_loc_[] = {0};
char TideMatchSet::decoy_match_collection_loc_[] = {0};

TideMatchSet::TideMatchSet(Arr* matches, double max_mz)
  : matches_(matches), max_mz_(max_mz), exact_pval_search_(false), elution_window_(0) {
}

TideMatchSet::TideMatchSet(Peptide* peptide, double max_mz)
  : peptide_(peptide), max_mz_(max_mz), exact_pval_search_(false), elution_window_(0) {
}

TideMatchSet::~TideMatchSet() {
}

/**
 * Write peptide centric matches to output files
 * This is for writing tab-delimited only
 */
void TideMatchSet::report(
  ofstream* target_file,  ///< target file to write to
  ofstream* decoy_file, ///< decoy file to write to
  int top_matches,
  const ActivePeptideQueue* peptides, ///< peptide queue
  const ProteinVec& proteins, ///< proteins corresponding with peptides
  const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
  bool compute_sp ///< whether to compute sp or not
) {
  if (peptide_->spectrum_matches_array.size() == 0) {
    return;
  }

  carp(CARP_DETAILED_DEBUG, "TideMatchSet reporting top %d of %d peptide centric matches",
       top_matches, peptide_->spectrum_matches_array.size());

  int charge;
  double score;
  double d_cn = 0.0;
  int nHit = peptide_->spectrum_matches_array.size();

  if (nHit < top_matches) {
      top_matches = nHit;
  }
  if (exact_pval_search_) {
    sort(peptide_->spectrum_matches_array.begin(),
         peptide_->spectrum_matches_array.end(),
         Peptide::spectrum_matches::compPV);
  } else {
    sort(peptide_->spectrum_matches_array.begin(),
         peptide_->spectrum_matches_array.end(),
         Peptide::spectrum_matches::compSC);
  } 
 
  for (int cnt = 0; cnt < nHit; ++cnt) {
    d_cn = 0.0;
    if (exact_pval_search_ == true) {
      score = peptide_->spectrum_matches_array[cnt].score1_;
      if (cnt < nHit-1) {
        d_cn = (double)((log10(peptide_->spectrum_matches_array[cnt+1].score1_)
                       - log10(peptide_->spectrum_matches_array[cnt].score1_))
                       /max((FLOAT_T)(-1*log10(peptide_->spectrum_matches_array[cnt].score1_)), FLOAT_T(1)));
      }
    } else {
      score = (double)(peptide_->spectrum_matches_array[cnt].score1_ / 100000000.0);
      if (cnt < nHit-1) {
        d_cn = (double)( score 
                      - (double)(peptide_->spectrum_matches_array[cnt+1].score1_ / 100000000.0)
                      / (double)max((FLOAT_T)score , FLOAT_T(1)));
      }
    }
    peptide_->spectrum_matches_array[cnt].score1_ = score;
    peptide_->spectrum_matches_array[cnt].d_cn_ = d_cn;
    peptide_->spectrum_matches_array[cnt].score3_ = nHit;
  }
  //smoothing primary scores in the elution window, only in DIA mode.
  if (elution_window_ > 0) {
    sort(peptide_->spectrum_matches_array.begin(), 
         peptide_->spectrum_matches_array.end(),
         Peptide::spectrum_matches::compRT);
    int cnt;
    double mean = 1.0;

    //initialize sliding window
    int flank = (int)((double)((elution_window_)/2) + 1);
    flank = flank > nHit ? nHit : flank;

    for (cnt = 0; cnt < flank; ++cnt) {
      mean *= peptide_->spectrum_matches_array[cnt].score1_;
    }
    int top = flank;
    int bottom = 0;
    for (cnt = 0; cnt < nHit; ++cnt) {
      peptide_->spectrum_matches_array[cnt].elution_score_ = pow(mean, 1.0/(top-bottom));
      if (top < nHit) {
        mean *= peptide_->spectrum_matches_array[top].score1_;
        ++top;
      }
      if (cnt >= flank-1) {
        mean /= peptide_->spectrum_matches_array[bottom].score1_;
        ++bottom;
      }
    }
    //reorder PSMs according to the smoothed p-value
    sort(peptide_->spectrum_matches_array.begin(),
         peptide_->spectrum_matches_array.end(),
         Peptide::spectrum_matches::compES);
  }
  peptide_->spectrum_matches_array.resize(top_matches);
  if (compute_sp) {
    vector<pair<double, int> > spScoreRank;
    spScoreRank.reserve(top_matches);
    for (int cnt = 0; cnt < top_matches; ++cnt) {  
      SpScorer sp_scorer(proteins, *peptide_->spectrum_matches_array[cnt].spectrum_, 
                         peptide_->spectrum_matches_array[cnt].charge_, max_mz_);
      pb::Peptide* pb_peptide = getPbPeptide(*peptide_);
      sp_scorer.Score(*pb_peptide, peptide_->spectrum_matches_array[cnt].spData_);
      spScoreRank.push_back(make_pair(-1*peptide_->spectrum_matches_array[cnt].spData_.sp_score, cnt));
    }
    sort(spScoreRank.begin(), spScoreRank.end());
    for (size_t i = 0; i < spScoreRank.size(); ++i) {
      peptide_->spectrum_matches_array[spScoreRank[i].second].spData_.sp_rank = i;
    }
  }  
  // target peptide or concat search
  ofstream* file =
    (OutputFiles::isConcat() || !peptide_->IsDecoy()) ? target_file : decoy_file;
  writeToFile(file, peptides, proteins, locations, compute_sp);
}

/**
 * Helper function for tab delimited report function for peptide centric search
 */
void TideMatchSet::writeToFile(
  ofstream* file,
  const ActivePeptideQueue* peptides,
  const ProteinVec& proteins,
  const vector<const pb::AuxLocation*>& locations,
  bool compute_sp ///< whether to compute sp or not
) {
  if (!file) {
    return;
  }
  int cur = 0;

  const Peptide* peptide = peptides->GetPeptide(0);
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

  Crux::Peptide cruxPep = getCruxPeptide(peptide);
  for (vector<Peptide::spectrum_matches>::const_iterator 
        i = peptide_->spectrum_matches_array.begin(); 
        i != peptide_->spectrum_matches_array.end(); 
        ++i) {
    Spectrum* spectrum = i->spectrum_;
    
    *file << spectrum->SpectrumNumber() << '\t'
          << i->charge_ << '\t'
          << spectrum->PrecursorMZ() << '\t'
          << (spectrum->PrecursorMZ() - MASS_PROTON) * i->charge_ << '\t'
          << cruxPep.calcModifiedMass() << '\t'
          << i->d_cn_ << '\t';
    SpScorer::SpScoreData spData;
    if (compute_sp) {
      *file << i->spData_.sp_score << '\t'
          << i->spData_.sp_rank << '\t';
    }
    *file << i->score1_ << '\t';
    if (exact_pval_search_) {
      *file << i->score2_<< '\t';
    }
    
    if (elution_window_ ) {
      *file << i->elution_score_ << '\t';
    }

    *file << ++cur << '\t';
    if (compute_sp) {
      *file << i->spData_.matched_ions << '\t'
            << i->spData_.total_ions << '\t';
    }
    *file << i->score3_ << '\t';

    if (OutputFiles::isConcat()) {
      *file << peptides->ActiveTargets() + peptides->ActiveDecoys() << '\t';
    } else {
      *file << (!peptide->IsDecoy() ? peptides->ActiveTargets() : peptides->ActiveDecoys()) << '\t';
    }
    *file << cruxPep.getModifiedSequenceWithMasses() << '\t'
          << cruxPep.getModsString() << '\t'
          << CleavageType << '\t'
          << proteinNames << '\t'
          << flankingAAs;
    if (peptide->IsDecoy() && !OutputFiles::isProteinLevelDecoys()) {
      // write target sequence
      const string& residues = protein->residues();
      *file << '\t'
            << residues.substr(residues.length() - peptide->Len());
    } else if (OutputFiles::isConcat() && !OutputFiles::isProteinLevelDecoys()) {
      *file << '\t'
            << cruxPep.getUnshuffledSequence();
    }
    *file << endl;
  }
}

/**
 * Write matches to output files
 * This is for writing tab-delimited only
 */
void TideMatchSet::report(
  ofstream* target_file,  ///< target file to write to
  ofstream* decoy_file, ///< decoy file to write to
  int top_n,  ///< number of matches to report
  const string& spectrum_filename, ///< name of spectrum file
  const Spectrum* spectrum, ///< spectrum for matches
  int charge, ///< charge for matches
  const ActivePeptideQueue* peptides, ///< peptide queue
  const ProteinVec& proteins,  ///< proteins corresponding with peptides
  const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
  bool compute_sp, ///< whether to compute sp or not
  bool highScoreBest, //< indicates semantics of score magnitude
  boost::mutex * rwlock
) {
  if (matches_->size() == 0) {
    return;
  }

  carp(CARP_DETAILED_DEBUG, "Tide MatchSet reporting top %d of %d matches",
       top_n, matches_->size());

  vector<Arr::iterator> targets, decoys;
  gatherTargetsAndDecoys(peptides, proteins, targets, decoys, top_n, highScoreBest);

  map<Arr::iterator, FLOAT_T> delta_cn_map;
  computeDeltaCns(targets, &delta_cn_map);
  computeDeltaCns(decoys, &delta_cn_map);

  map<Arr::iterator, pair<const SpScorer::SpScoreData, int> > sp_map;
  if (compute_sp) {
    SpScorer sp_scorer(proteins, *spectrum, charge, max_mz_);
    computeSpData(targets, &sp_map, &sp_scorer, peptides);
    computeSpData(decoys, &sp_map, &sp_scorer, peptides);
  }
  writeToFile(target_file, top_n, targets, spectrum_filename, spectrum, charge,
              peptides, proteins, locations, delta_cn_map, compute_sp ? &sp_map : NULL, rwlock);
  writeToFile(decoy_file, top_n, decoys, spectrum_filename, spectrum, charge,
              peptides, proteins, locations, delta_cn_map, compute_sp ? &sp_map : NULL, rwlock);
}

/**
 * Helper function for tab delimited report function
 */
void TideMatchSet::writeToFile(
  ofstream* file,
  int top_n,
  const vector<Arr::iterator>& vec,
  const string& spectrum_filename,
  const Spectrum* spectrum,
  int charge,
  const ActivePeptideQueue* peptides,
  const ProteinVec& proteins,
  const vector<const pb::AuxLocation*>& locations,
  const map<Arr::iterator, FLOAT_T>& delta_cn_map,
  const map<Arr::iterator, pair<const SpScorer::SpScoreData, int> >* sp_map,
  boost::mutex * rwlock
) {
  if (!file) {
    return;
  }

  int massPrecision = Params::GetInt("mass-precision");
  int precision = Params::GetInt("precision");

  int cur = 0;
  int concatDistinctMatches = peptides->ActiveTargets() + peptides->ActiveDecoys();

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

    Crux::Peptide cruxPep = getCruxPeptide(peptide);
    const SpScorer::SpScoreData* sp_data = sp_map ? &(sp_map->at(*i).first) : NULL;

    rwlock->lock();
    if (Params::GetBool("file-column")) {
      *file << spectrum_filename << '\t';
    }
    *file << spectrum->SpectrumNumber() << '\t'
          << charge << '\t'
          << StringUtils::ToString(spectrum->PrecursorMZ(), massPrecision) << '\t'
          << StringUtils::ToString((spectrum->PrecursorMZ() - MASS_PROTON) * charge, massPrecision) << '\t'
          << StringUtils::ToString(cruxPep.calcModifiedMass(), massPrecision) << '\t'
          << delta_cn_map.at(*i) << '\t';
    if (sp_map) {
      *file << StringUtils::ToString(sp_data->sp_score, precision) << '\t'
            << sp_map->at(*i).second << '\t';
    }
    *file << StringUtils::ToString((*i)->first.first, precision) << '\t';
    if (exact_pval_search_) {
      *file << (*i)->first.second << '\t';
    }
    *file << ++cur << '\t';
    if (sp_map) {
      *file << sp_data->matched_ions << '\t'
            << sp_data->total_ions << '\t';
    }

    if (OutputFiles::isConcat()) {
      *file << concatDistinctMatches << '\t';
    } else {
      *file << (!peptide->IsDecoy() ? peptides->ActiveTargets() : peptides->ActiveDecoys()) << '\t';
    }

    *file << cruxPep.getModifiedSequenceWithMasses() << '\t'
          << cruxPep.getModsString() << '\t'
          << CleavageType << '\t'
          << proteinNames << '\t'
          << flankingAAs;
    if (peptide->IsDecoy() && !OutputFiles::isProteinLevelDecoys()) {
      // write target sequence
      const string& residues = protein->residues();
      *file << '\t'
            << residues.substr(residues.length() - peptide->Len());
    } else if (OutputFiles::isConcat() && !OutputFiles::isProteinLevelDecoys()) {
      *file << '\t'
            << cruxPep.getUnshuffledSequence();
    }
    *file << endl;
    rwlock->unlock();
  }
}

/**
 * Write matches to output files
 */
void TideMatchSet::report(
  OutputFiles* output_files,  ///< pointer to output handler
  int top_n,  ///< number of matches to report
  const string& spectrum_filename, ///< name of spectrum file
  const Spectrum* spectrum, ///< spectrum for matches
  int charge, ///< charge for matches
  const ActivePeptideQueue* peptides, ///< peptide queue
  const ProteinVec& proteins,  ///< proteins corresponding with peptides
  const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
  bool compute_sp, ///< whether to compute sp or not
  bool highScoreBest // indicates semantics of score magnitude
) {
  if (matches_->size() == 0) {
    return;
  }

  carp(CARP_DETAILED_DEBUG, "Tide MatchSet reporting top %d of %d matches",
       top_n, matches_->size());

  vector<Arr::iterator> targets, decoys;
  gatherTargetsAndDecoys(peptides, proteins, targets, decoys, top_n, highScoreBest);

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
    spectrum->PrecursorMZ(), vector<int>(1, charge), spectrum_filename);
  SpectrumZState z_state;
  z_state.setMZ(crux_spectrum.getPrecursorMz(), charge);

  crux_collection->exact_pval_search_ = exact_pval_search_;
  crux_decoy_collection->exact_pval_search_ = exact_pval_search_;

  addCruxMatches(crux_collection, false, top_n, &proteins_made, targets, crux_spectrum,
                 peptides, proteins, locations, z_state, sp_scorer, &lowest_sp);
  addCruxMatches(crux_decoy_collection, true, top_n, &proteins_made, decoys, crux_spectrum,
                 peptides, proteins, locations, z_state, sp_scorer, &lowest_sp);
  crux_collection->setFilePath(spectrum_filename, true);
  crux_decoy_collection->setFilePath(spectrum_filename, true);

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
  SCORER_TYPE_T scoreType =
    !Params::GetBool("exact-p-value") ? XCORR : TIDE_SEARCH_EXACT_PVAL;
  output_files->writeMatches(crux_collection, decoy_vector, scoreType, &crux_spectrum);

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
  FLOAT_T lnNumSp = OutputFiles::isConcat()
    ? log((FLOAT_T) (peptides->ActiveTargets() + peptides->ActiveDecoys()))
    : log((FLOAT_T) (!decoys ? peptides->ActiveTargets() : peptides->ActiveDecoys()));

  bool exactPValue = Params::GetBool("exact-p-value");
  
  // Create a Crux match for each match
  vector<Arr::iterator>::const_iterator endIter = 
    (vec.size() >= top_n + 1) ? vec.begin() + top_n + 1 : vec.end();
  for (vector<Arr::iterator>::const_iterator i = vec.begin(); i != endIter; ++i) {
    const Peptide* peptide = peptides->GetPeptide((*i)->second);

    Crux::Match* match = getCruxMatch(peptide, proteins, locations, &crux_spectrum,
                                      z_state, proteins_made);
    match_collection->addMatch(match);
    if (peptide->IsDecoy()) {
      match->setNullPeptide(true);
    }
    Crux::Match::freeMatch(match); // so match gets deleted when collection does

    // Set Xcorr score in match
    if (!exactPValue) {
      match->setScore(XCORR, (*i)->first.first);
    } else {
      match->setScore(TIDE_SEARCH_EXACT_PVAL, (*i)->first.first);
      match->setScore(TIDE_SEARCH_REFACTORED_XCORR, (*i)->first.second);
    }

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
      match->setScore(BY_IONS_MATCHED, sp_score_data.matched_ions);
      match->setScore(BY_IONS_TOTAL, sp_score_data.total_ions);
    }
  }
  match_collection->setZState(z_state);
  if (!OutputFiles::isConcat()) {
    match_collection->setExperimentSize(!decoys ? peptides->ActiveTargets() : peptides->ActiveDecoys());
  } else {
    match_collection->setExperimentSize(peptides->ActiveTargets() + peptides->ActiveDecoys());
  }
  if (sp_scorer) {
    match_collection->setScoredType(SP, true);
    match_collection->setScoredType(BY_IONS_MATCHED, true);
    match_collection->setScoredType(BY_IONS_TOTAL, true);
    match_collection->populateMatchRank(SP);
  }
  if (!exactPValue) {
    match_collection->setScoredType(XCORR, true);
    match_collection->populateMatchRank(XCORR);
  } else {
    match_collection->setScoredType(TIDE_SEARCH_REFACTORED_XCORR, true);
    match_collection->populateMatchRank(TIDE_SEARCH_REFACTORED_XCORR);
    match_collection->setScoredType(TIDE_SEARCH_EXACT_PVAL, true);
    match_collection->populateMatchRank(TIDE_SEARCH_EXACT_PVAL);
  }
}

/**
 * Write headers for tab delimited file
 */
void TideMatchSet::writeHeaders(ofstream* file, bool decoyFile, bool sp) {
  if (!file) {
    return;
  }
  const int headers[] = {
    FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, SP_SCORE_COL, SP_RANK_COL,
    XCORR_SCORE_COL, XCORR_RANK_COL, BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL,
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, CLEAVAGE_TYPE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, ORIGINAL_TARGET_SEQUENCE_COL
  };
  size_t numHeaders = sizeof(headers) / sizeof(int);
  bool writtenHeader = false;
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
    if (writtenHeader) {
      *file << '\t';
    }
    if (header == FILE_COL &&
        (!Params::GetBool("file-column") || Params::GetBool("peptide-centric-search"))) {
      continue;
    }
    if (header == XCORR_SCORE_COL && Params::GetBool("exact-p-value")) {
      *file << get_column_header(EXACT_PVALUE_COL) << '\t'
            << get_column_header(REFACTORED_SCORE_COL);
      if (Params::GetInt("elution-window-size") > 0) {
        *file << '\t' << get_column_header(ELUTION_WINDOW_COL);
      }
      writtenHeader = true;
      continue;
    }
    if (header == DISTINCT_MATCHES_SPECTRUM_COL) {
      if (Params::GetBool("peptide-centric-search")) {
        *file << get_column_header(DISTINCT_MATCHES_PEPTIDE_COL) << '\t';
        *file << get_column_header(DISTINCT_MATCHES_SPECTRUM_COL);
      } else {
        *file << get_column_header(DISTINCT_MATCHES_SPECTRUM_COL);
      }
      writtenHeader = true;
      continue;
    }
    *file << get_column_header(header);
    writtenHeader = true;
  }
  *file << endl;
}

void TideMatchSet::initModMap(const pb::ModTable& modTable, ModPosition position) {
  for (int i = 0; i < modTable.variable_mod_size(); i++) {
    const pb::Modification& mod = modTable.variable_mod(i);
    if (mod.has_delta() && mod.has_amino_acids()) {
      ModificationDefinition::NewVarMod(mod.amino_acids(), mod.delta(), position);
    }
  }
  for (int i = 0; i < modTable.static_mod_size(); i++) {
    const pb::Modification& mod = modTable.static_mod(i);
    if (mod.has_delta() && mod.has_amino_acids()) {
      ModificationDefinition::NewStaticMod(mod.amino_acids(), mod.delta(), position);
    }
  }
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
    protein->residues().substr(protein->residues().length() - peptide->Len());
  int start_idx = crux_protein->findStart(originalTargetSeq, n_term, c_term);

  // Create peptide
  Crux::Peptide* crux_peptide = new Crux::Peptide(peptide->Len(), crux_protein, start_idx);
  crux_peptide->getPeptideSrc()->setStartIdxOriginal(
    ((!protein->has_target_pos()) ? pos : protein->target_pos()) + 1);
  crux_peptide->setMods(getMods(peptide));

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
      src->setStartIdxOriginal(((!protein->has_target_pos()) ? pos : protein->target_pos()) + 1);
      crux_peptide->addPeptideSrc(src);
    }
  }

  return new Crux::Match(crux_peptide, crux_spectrum, crux_z_state, false);
}

Crux::Peptide TideMatchSet::getCruxPeptide(const Peptide* peptide) {
  return Crux::Peptide(peptide->Seq(), getMods(peptide));
}

vector<Crux::Modification> TideMatchSet::getMods(const Peptide* peptide) {
  vector<Crux::Modification> modVector;
  string seq(peptide->Seq());
  const ModCoder::Mod* mods;
  int pep_mods = peptide->Mods(&mods);
  for (int i = 0; i < pep_mods; i++) {
    int mod_index;
    double mod_delta;
    MassConstants::DecodeMod(mods[i], &mod_index, &mod_delta);
    const ModificationDefinition* modDef = ModificationDefinition::Find(mod_delta, false);
    if (modDef == NULL) {
      carp(CARP_ERROR, "Could not find modification with delta %f", mod_delta);
      continue;
    }
    modVector.push_back(Crux::Modification(modDef, mod_index));
  }
  return modVector;
}

void TideMatchSet::gatherTargetsAndDecoys(
  const ActivePeptideQueue* peptides,
  const ProteinVec& proteins,
  vector<Arr::iterator>& targetsOut,
  vector<Arr::iterator>& decoysOut,
  int top_n,
  bool highScoreBest // indicates semantics of score magnitude
) {
  make_heap(matches_->begin(), matches_->end(), highScoreBest ? lessScore : moreScore);
  if (!OutputFiles::isConcat() && TideSearchApplication::hasDecoys()) {
    for (Arr::iterator i = matches_->end(); i != matches_->begin(); ) {
      pop_heap(matches_->begin(), i--, highScoreBest ? lessScore : moreScore);
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
      pop_heap(matches_->begin(), matches_->end() - i, highScoreBest ? lessScore : moreScore);
      targetsOut.push_back(matches_->end() - (++i));
    }
  }
}

/**
 * Create a pb peptide from Tide peptide
 */
pb::Peptide* TideMatchSet::getPbPeptide(const Peptide& peptide) {
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
string TideMatchSet::getProteinName(const pb::Protein& protein, int pos) {
  stringstream proteinNameStream;
  proteinNameStream << protein.name() << '(' << pos + 1 << ')';
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
  map<Arr::iterator, FLOAT_T>* delta_cn_map // map to add delta cn scores to
) {
  vector<FLOAT_T> scores;
  for (vector<Arr::iterator>::const_iterator i = vec.begin(); i != vec.end(); i++) {
    scores.push_back((*i)->first.first);
  }
  vector< pair<FLOAT_T, FLOAT_T> > deltaCns = MatchCollection::calculateDeltaCns(
    scores, !Params::GetBool("exact-p-value") ? XCORR : TIDE_SEARCH_EXACT_PVAL);
  for (int i = 0; i < vec.size(); i++) {
    delta_cn_map->insert(make_pair(vec[i], deltaCns[i].first));
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
