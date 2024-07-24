/*************************************************************************
 * \file MatchFileReader.cpp
 * \brief Object for parsing the tab-delimited files
 *************************************************************************/

#include "MatchColumns.h"
#include "MatchFileReader.h"
#include "DelimitedFile.h"

#include "model/MatchCollection.h"
#include "model/Modification.h"

using namespace std;

/**
 * \returns a blank MatchFileReader object
 */
MatchFileReader::MatchFileReader() : DelimitedFileReader(), PSMReader() {
}

/**
 * \returns a MatchFileReader object and loads the tab-delimited
 * data specified by file_name.
 */
MatchFileReader::MatchFileReader(const char* file_name) : DelimitedFileReader(file_name, true) {
  parseHeader();
}

/**
 * \returns a MatchFileReader object and loads the tab-delimited
 * data specified by file_name.
 */
MatchFileReader::MatchFileReader(const string& file_name)
  : DelimitedFileReader(file_name, true), PSMReader(file_name) {
  parseHeader();
}

MatchFileReader::MatchFileReader(const string& file_name, Database* database, Database* decoy_database)
  : DelimitedFileReader(file_name, true), PSMReader(file_name, database, decoy_database) {
  parseHeader();
}

MatchFileReader::MatchFileReader(istream* iptr) : DelimitedFileReader(iptr, true, '\t') {
  parseHeader();
}

/**
 * Destructor
 */
MatchFileReader::~MatchFileReader() {
}

/**
 * Open a new file from an existing MatchFileReader.
 */
void MatchFileReader::loadData(
  const char* file_name, ///< new file to open
  bool hasHeader
) {
  DelimitedFileReader::loadData(file_name, hasHeader);
  if( hasHeader ) {
    parseHeader();
  }
}

/**
 * Open a new file from an existing MatchFileReader.
 */
void MatchFileReader::loadData(
  const string& file_name, ///< new file to open
  bool hasHeader) {
  DelimitedFileReader::loadData(file_name, hasHeader);
  if( hasHeader ) {
    parseHeader();
  }
}

/**
 * parses the header and builds the internal hash table
 */
void MatchFileReader::parseHeader() {
  for (int idx = 0; idx < NUMBER_MATCH_COLUMNS; idx++) {
    match_indices_[idx] = findColumn(get_column_header(idx));
  }
}

/**
 * \returns the FLOAT_T value of a cell, checks for infinity
 */
FLOAT_T MatchFileReader::getFloat(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getFloat", get_column_header(col_type));
    return -1;
  }
  return DelimitedFileReader::getFloat(idx);
}

/**
 * \returns the double value of a cell, checks for infinity
 */
double MatchFileReader::getDouble(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  carp(CARP_DETAILED_DEBUG, "reading double from column %s", get_column_header(col_type));
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getDouble", get_column_header(col_type));
    return -1;
  }
  return DelimitedFileReader::getDouble(idx);
}


/**
 * \returns the integer value of a cell, checks for infinity.
 */
int MatchFileReader::getInteger(
  MATCH_COLUMNS_T col_type ///< the column name
) {
  carp(CARP_DETAILED_DEBUG, "Reading integer from column %s", get_column_header(col_type));
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getInteger", get_column_header(col_type));
    return -1;
  }
  return DelimitedFileReader::getInteger(idx);
}

/**
 * \returns the string value of a cell
 */
string MatchFileReader::getString(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  carp(CARP_DETAILED_DEBUG, "Getting string from column %s", get_column_header(col_type));
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getString", get_column_header(col_type));
    return "";
  }
  return DelimitedFileReader::getString(idx);
}

bool MatchFileReader::empty(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  int idx = match_indices_[col_type];
  if (idx == -1) {
    return true;
  }
  return DelimitedFileReader::getString(idx).empty();
}

/**
 * Fills in the given vector with a bool value indicating if each
 * MATCH_COLUMN_T type is present in the file being read.
 * \returns Argument vector has NUM_MATCH_COLUMN_T values if a
 * valid file is open and header has been parsed, else vector is empty.
 */
void MatchFileReader::getMatchColumnsPresent(
  std::vector<bool>& col_is_present) {
  col_is_present.clear();

  // has a header been parsed?
  if( column_names_.empty() ) {
    return;
  }
  col_is_present.assign(NUMBER_MATCH_COLUMNS, false);

  for(int col_idx = 0; col_idx < NUMBER_MATCH_COLUMNS; col_idx++) {
    col_is_present[col_idx] = (match_indices_[col_idx] > -1);
  }
}

// FIXME: Need to generalize this to work with Percolator files.
MatchCollection* MatchFileReader::parse(
  const string& file_path,
  Database* database,
  Database* decoy_database) {
  return MatchFileReader(file_path, database, decoy_database).parse();
}

MatchCollection* MatchFileReader::parse() {
  MatchCollection* match_collection = new MatchCollection();
  match_collection->preparePostProcess();
  int maxRank = Params::GetInt("top-match-in");

  while (hasNext()) {
    FLOAT_T ln_experiment_size = 0;
    if (!empty(DISTINCT_MATCHES_SPECTRUM_COL)) {
      match_collection->setHasDistinctMatches(true);
      ln_experiment_size = log(getFloat(DISTINCT_MATCHES_SPECTRUM_COL));
    } else if (!empty(MATCHES_SPECTRUM_COL)) {
      ln_experiment_size = log(getFloat(MATCHES_SPECTRUM_COL));
    }

    match_collection->setScoredType(DELTA_CN, !empty(DELTA_CN_COL));
    match_collection->setScoredType(DELTA_LCN, !empty(DELTA_LCN_COL));
    match_collection->setScoredType(SP, !empty(SP_SCORE_COL));
    match_collection->setScoredType(XCORR, !empty(XCORR_SCORE_COL));
    match_collection->setScoredType(TIDE_SEARCH_EXACT_PVAL, !empty(EXACT_PVALUE_COL));
    match_collection->setScoredType(TIDE_SEARCH_REFACTORED_XCORR, !empty(REFACTORED_SCORE_COL));
    match_collection->setScoredType(RESIDUE_EVIDENCE_PVAL, !empty(RESIDUE_PVALUE_COL)); //Added by Andy Lin for residue evidence
    match_collection->setScoredType(RESIDUE_EVIDENCE_SCORE, !empty(RESIDUE_EVIDENCE_COL)); //Added by Andy Lin for residue evidence
    match_collection->setScoredType(BOTH_PVALUE, !empty(BOTH_PVALUE_COL)); //Added by Andy Lin for residue evidence
    match_collection->setScoredType(EVALUE, !empty(EVALUE_COL));
    match_collection->setScoredType(DECOY_XCORR_QVALUE, !empty(DECOY_XCORR_QVALUE_COL));
    match_collection->setScoredType(PERCOLATOR_QVALUE, !empty(POUT_QVALUE_COL));
    if (empty(POUT_QVALUE_COL)) {
      match_collection->setScoredType(PERCOLATOR_QVALUE, !empty(PERCOLATOR_QVALUE_COL));
    }
    match_collection->setScoredType(PERCOLATOR_SCORE, !empty(POUT_SCORE_COL));
    if (empty(POUT_SCORE_COL)) {
      match_collection->setScoredType(PERCOLATOR_SCORE, !empty(PERCOLATOR_SCORE_COL));
    }
    match_collection->setScoredType(BY_IONS_MATCHED, !empty(BY_IONS_MATCHED_COL));
    match_collection->setScoredType(BY_IONS_TOTAL, !empty(BY_IONS_TOTAL_COL));
    match_collection->setScoredType(BY_ION_FRACTION, !empty(  BY_IONS_FRACTION_COL));
    match_collection->setScoredType(BY_ION_REPEAT_MATCH, !empty(BY_IONS_REPEAT_MATCH_COL));
    match_collection->setScoredType(TAILOR_SCORE, !empty(TAILOR_COL)); //Added for tailor score calibration method by AKF
    match_collection->setScoredType(QVALUE_TDC, !empty(QVALUE_TDC_COL));    

    // DIAmeter related, added by Yang
    match_collection->setScoredType(PRECURSOR_INTENSITY_RANK_M0, !empty(PRECURSOR_INTENSITY_RANK_M0_COL));
    match_collection->setScoredType(PRECURSOR_INTENSITY_RANK_M1, !empty(PRECURSOR_INTENSITY_RANK_M1_COL));
    match_collection->setScoredType(PRECURSOR_INTENSITY_RANK_M2, !empty(PRECURSOR_INTENSITY_RANK_M2_COL));
    match_collection->setScoredType(RT_DIFF, !empty(RT_DIFF_COL));
    match_collection->setScoredType(DYN_FRAGMENT_PVALUE, !empty(DYN_FRAGMENT_PVALUE_COL));
    match_collection->setScoredType(STA_FRAGMENT_PVALUE, !empty(STA_FRAGMENT_PVALUE_COL));
    match_collection->setScoredType(COELUTE_MS1, !empty(COELUTE_MS1_COL));
    match_collection->setScoredType(COELUTE_MS2, !empty(COELUTE_MS2_COL));
    match_collection->setScoredType(COELUTE_MS1_MS2, !empty(COELUTE_MS1_MS2_COL));
    match_collection->setScoredType(ENSEMBLE_SCORE, !empty(ENSEMBLE_SCORE_COL));

    if (!empty(DECOY_INDEX_COL)) {
      match_collection->setHasDecoyIndexes(true);
    }

    // TODO presumably we can do this once instead of one time per scan
    MATCH_COLUMNS_T rank_col;
    if (empty(PERCOLATOR_RANK_COL) == 1) {
      rank_col = PERCOLATOR_RANK_COL;
    } else if (empty(BOTH_PVALUE_RANK) == 1) {
      rank_col = BOTH_PVALUE_RANK;
    } else if (empty(RESIDUE_RANK_COL) == 1) {
      rank_col = RESIDUE_RANK_COL;
    } else if (empty(XCORR_RANK_COL) == 1) {
      rank_col = XCORR_RANK_COL;
    } else if (empty(SP_RANK_COL) == 1) {
      rank_col = SP_RANK_COL;
    } else {
      carp(CARP_FATAL, "Input file does not contain any reconized rank column.");
    }

    // parse match object
    if (maxRank == 0 || getInteger(rank_col) <= maxRank) {
      Crux::Match* match = parseMatch();
      if (match == NULL) {
        carp(CARP_ERROR, "Failed to parse tab-delimited PSM match");
        return NULL;
      }
      // set all spectrum specific features to parsed match
      // zState also gets set by spectrum = parseSpectrum() and spectrum is used
      // as input for match. This code changes zState only if charge and neutral
      // mass column is present. Note that this if statement is required as
      // Percolator PSMId column is parsed by parseSpectrum()
      if (getFloat(SPECTRUM_NEUTRAL_MASS_COL) != -1 &&
          getInteger(CHARGE_COL) != -1) {
        SpectrumZState zState(getFloat(SPECTRUM_NEUTRAL_MASS_COL),
                              getInteger(CHARGE_COL));
        match->setZState(zState);
      }

      // //set all spectrum specific features to parsed match
      // SpectrumZState zState(getFloat(SPECTRUM_NEUTRAL_MASS_COL),
      //                       getInteger(CHARGE_COL));
      // match->setZState(zState);
      match->setLnExperimentSize(ln_experiment_size);
      //add match to match collection.
      match_collection->addMatchToPostMatchCollection(match);
    }
    //increment pointer.
    next();
  }

  return match_collection;
}

/**
 *\returns a match object that is parsed from the tab-delimited result file
 */
Crux::Match* MatchFileReader::parseMatch() {
  // parse spectrum
  Crux::Spectrum* spectrum = parseSpectrum();
  if (spectrum == NULL) {
    carp(CARP_ERROR, "Failed to parse spectrum (tab delimited).");
  }

  // parse peptide
  Crux::Peptide* peptide = parsePeptide();
  if (peptide == NULL) {
    carp(CARP_ERROR, "Failed to parse peptide (tab delimited)");
    // FIXME should this exit or return null. I think sometimes we can get
    // no peptides, which is valid, in which case NULL makes sense.
    // maybe this should be fixed at the output match level however.
    return NULL;
  }

  Crux::Match* match = new Crux::Match(peptide, spectrum, spectrum->getZState(0), false);
  match->setPostProcess(true);
  match->setDatabaseIndexName(getString(INDEX_NAME_COL));
  string path = spectrum->getFullFilename();
  if (!path.empty()) {
    match->setFilePath(path);
  }

  if (!empty(FILE_COL)) {
    match->setFilePath(getString(FILE_COL));
  }

  if (empty(SP_SCORE_COL) || empty(SP_RANK_COL)) {
    match->setScore(SP, NOT_SCORED);
    match->setRank(SP, 0);
  } else {
    match->setScore(SP, getFloat(SP_SCORE_COL));
    match->setRank(SP, getInteger(SP_RANK_COL));
  }

  match->setScore(XCORR, getFloat(XCORR_SCORE_COL));
  match->setRank(XCORR, getInteger(XCORR_RANK_COL));

  if (!empty(DELTA_CN_COL)) {
    match -> setScore(DELTA_CN, getFloat(DELTA_CN_COL));
  }
  if (!empty(DELTA_LCN_COL)) {
    match -> setScore(DELTA_LCN, getFloat(DELTA_LCN_COL));
  }

  if (!empty(EXACT_PVALUE_COL)) {
    match->setScore(TIDE_SEARCH_EXACT_PVAL, getFloat(EXACT_PVALUE_COL));
    match->setScore(TIDE_SEARCH_REFACTORED_XCORR, getFloat(REFACTORED_SCORE_COL));
    match->setRank(TIDE_SEARCH_EXACT_PVAL, getInteger(XCORR_RANK_COL));
  }
  if (!empty(RESIDUE_EVIDENCE_COL)) {
    match->setScore(RESIDUE_EVIDENCE_SCORE, getFloat(RESIDUE_EVIDENCE_COL));
    match->setRank(RESIDUE_EVIDENCE_SCORE, getInteger(RESIDUE_RANK_COL));
  }
  if (!empty(RESIDUE_PVALUE_COL)) {
    match->setScore(RESIDUE_EVIDENCE_SCORE, getFloat(RESIDUE_EVIDENCE_COL));
    match->setScore(RESIDUE_EVIDENCE_PVAL, getFloat(RESIDUE_PVALUE_COL));
    match->setRank(RESIDUE_EVIDENCE_PVAL, getInteger(RESIDUE_RANK_COL));
  }
  if (!empty(BOTH_PVALUE_COL)) {
    match->setScore(BOTH_PVALUE, getFloat(BOTH_PVALUE_COL));
    match->setRank(BOTH_PVALUE, getInteger(BOTH_PVALUE_RANK));
  }
  if (!empty(DECOY_XCORR_QVALUE_COL)) {
    match->setScore(DECOY_XCORR_QVALUE, getFloat(DECOY_XCORR_QVALUE_COL));
  }
  if (!empty(EVALUE_COL)) {
    match->setScore(EVALUE, getFloat(EVALUE_COL));
  }
  if (!empty(PERCOLATOR_QVALUE_COL)) {
    match->setScore(PERCOLATOR_QVALUE, getFloat(PERCOLATOR_QVALUE_COL));
  } else if (!empty(POUT_QVALUE_COL)) {
    match->setScore(PERCOLATOR_QVALUE, getFloat(POUT_QVALUE_COL));
  }
  if (!empty(PERCOLATOR_SCORE_COL)) {
    match->setScore(PERCOLATOR_SCORE, getFloat(PERCOLATOR_SCORE_COL));
    match->setRank(PERCOLATOR_SCORE, getInteger(PERCOLATOR_RANK_COL));
  } else if (!empty(POUT_SCORE_COL)) {
    match->setScore(PERCOLATOR_SCORE, getFloat(POUT_SCORE_COL));
  }
  if (!empty(BY_IONS_MATCHED_COL)) {
    match->setScore(BY_IONS_MATCHED, getInteger(BY_IONS_MATCHED_COL));
  }
  if (!empty(BY_IONS_TOTAL_COL)) {
    match->setScore(BY_IONS_TOTAL, getInteger(BY_IONS_TOTAL_COL));
  }
  if (!empty(BY_IONS_FRACTION_COL)) {
    match->setScore(BY_ION_FRACTION, getFloat(BY_IONS_FRACTION_COL));
  }
  if (!empty(BY_IONS_REPEAT_MATCH_COL)) {
    match->setScore(BY_ION_REPEAT_MATCH, getInteger(BY_IONS_REPEAT_MATCH_COL));
  }
  if (!empty(DECOY_INDEX_COL)) {
    match->setDecoyIndex(getInteger(DECOY_INDEX_COL));
  }
  if (!empty(TAILOR_COL)) { //Added for tailor score calibration method by AKF
    match->setScore(TAILOR_SCORE, getFloat(TAILOR_COL));
  }
  if (!empty(QVALUE_TDC_COL)) {
    match->setScore(QVALUE_TDC, getFloat(QVALUE_TDC_COL));
  }  

  // DIAmeter related, added by Yang
  if (!empty(PRECURSOR_INTENSITY_RANK_M0_COL)) { match->setScore(PRECURSOR_INTENSITY_RANK_M0, getFloat(PRECURSOR_INTENSITY_RANK_M0_COL)); }
  if (!empty(PRECURSOR_INTENSITY_RANK_M1_COL)) { match->setScore(PRECURSOR_INTENSITY_RANK_M1, getFloat(PRECURSOR_INTENSITY_RANK_M1_COL)); }
  if (!empty(PRECURSOR_INTENSITY_RANK_M2_COL)) { match->setScore(PRECURSOR_INTENSITY_RANK_M2, getFloat(PRECURSOR_INTENSITY_RANK_M2_COL)); }
  if (!empty(RT_DIFF_COL)) { match->setScore(RT_DIFF, getFloat(RT_DIFF_COL)); }
  if (!empty(DYN_FRAGMENT_PVALUE_COL)) { match->setScore(DYN_FRAGMENT_PVALUE, getFloat(DYN_FRAGMENT_PVALUE_COL)); }
  if (!empty(STA_FRAGMENT_PVALUE_COL)) { match->setScore(STA_FRAGMENT_PVALUE, getFloat(STA_FRAGMENT_PVALUE_COL)); }
  if (!empty(COELUTE_MS1_COL)) { match->setScore(COELUTE_MS1, getFloat(COELUTE_MS1_COL)); }
  if (!empty(COELUTE_MS2_COL)) { match->setScore(COELUTE_MS2, getFloat(COELUTE_MS2_COL)); }
  if (!empty(COELUTE_MS1_MS2_COL)) { match->setScore(COELUTE_MS1_MS2, getFloat(COELUTE_MS1_MS2_COL)); }
  if (!empty(ENSEMBLE_SCORE_COL)) { match->setScore(ENSEMBLE_SCORE, getFloat(ENSEMBLE_SCORE_COL)); }


  // get experiment size
  int experimentSize = 0;
  if (!empty(DISTINCT_MATCHES_SPECTRUM_COL)) {
    experimentSize = getInteger(DISTINCT_MATCHES_SPECTRUM_COL);
  } else if (!empty(MATCHES_SPECTRUM_COL)) {
    experimentSize = getInteger(MATCHES_SPECTRUM_COL);
  }
  match->setTargetExperimentSize(experimentSize);
  if (experimentSize == 0) {
    match->setLnExperimentSize(0);
  } else {
    match->setLnExperimentSize(log((FLOAT_T) experimentSize));
  }

  if (!empty(PROTEIN_ID_COL) &&
      StringUtils::StartsWith(getString(PROTEIN_ID_COL), Params::GetString("decoy-prefix"))) {
    match->setNullPeptide(true);
  }

  return match;
}

Crux::Peptide* MatchFileReader::parsePeptide() {
  Crux::Peptide* peptide = NULL;
  string seq = getString(SEQUENCE_COL);
  bool is_percolator = false;
  if (seq.empty()) {
    seq = getString(POUT_PERC_PEPTIDE_COL);
    is_percolator = true;
  }
  if (seq.empty()) {
    carp(CARP_FATAL, "No peptide sequence found.");
  }

  // In cases where the sequence is in X.seq.X format, parse out the seq part
  if (seq.length() > 4 && seq[1] == '.' && seq[seq.length() - 2] == '.') {
    seq = seq.substr(2, seq.length() - 4);
  }

  peptide = new Crux::Peptide();

  string unmodSeq = Crux::Peptide::unmodifySequence(seq);
  vector<Crux::Modification> mods;

  // Parse modifications column first! It has more details about modifications.
  string modsString = getString(MODIFICATIONS_COL);
  if (!modsString.empty()) {
    mods = Crux::Modification::Parse(modsString, &unmodSeq);
  } else {
    Crux::Modification::FromSeq(seq, NULL, &mods);
  }

  peptide->setUnmodifiedSequence(unmodSeq);
  peptide->setMods(mods);

  if (!PeptideSrc::parseTabDelimited(peptide, *this, database_, decoy_database_)) {
    carp(CARP_ERROR, "Failed to parse peptide source.");
    delete peptide;
    return NULL;
  }

  return peptide;
}

Crux::Spectrum* MatchFileReader::parseSpectrum() {
  
  if (getInteger(SCAN_COL) != -1) {
    Crux::Spectrum* new_spectrum = new Crux::Spectrum(getInteger(SCAN_COL),
                              getInteger(SCAN_COL),
                              getFloat(SPECTRUM_PRECURSOR_MZ_COL),
                              vector<int>(1, getInteger(CHARGE_COL)),
                              getString(FILE_COL));
    new_spectrum->setRTime(getFloat(RETENTION_TIME_COL));
    return new_spectrum;
  } else {
    // This part allows Percolator output files to be read in
    string intermediate;
    vector <string> tokens;
    stringstream psmid(getString(POUT_PSMID_COL));

    // tokenize by underscore
    while(getline(psmid, intermediate, '_'))
    {
        tokens.push_back(intermediate);
    }
    carp(CARP_INFO, "chrage state: %d", stoi(tokens[3]));
    Crux::Spectrum* new_spectrum = new Crux::Spectrum(stoi(tokens[2]),
                              stoi(tokens[2]),
                              getFloat(SPECTRUM_PRECURSOR_MZ_COL),
                              vector<int>(1, stoi(tokens[3])),
                              getString(FILE_COL));
    new_spectrum->setRTime(getFloat(RETENTION_TIME_COL));
    return new_spectrum;
  }
    
}
