#include "HTMLWriter.h"

using namespace Crux;
using namespace boost;
using namespace std;

void HTMLWriter::writeRow() {
  if( current_row_.empty() ){
    return;
  }

  // make the row as long as the header
  while( current_row_.size() < column_names_.size()){
    current_row_.push_back("");
  }
  // TODO? warning if row is longer than non-empty header?

  // print each value separated by delimiter
  *file_ptr_ << "\t<tr>" << endl
             << "\t\t<td>" << current_row_[0] << "</td>" << endl;
  for(size_t idx = 1; idx < current_row_.size(); idx++){
    *file_ptr_ << "\t\t<td>" << current_row_[idx] << "</td>" << endl;
  }
  // end with newline
  *file_ptr_ << "\t</tr>" << endl;

  // clear the current_row and refill with blanks
  // if there is a header, that is the min length
  // if not, each row can be a different length
  current_row_.assign(column_names_.size(), "");

}

/**
 * Writes the data in a ProteinMatchCollection to the currently open file
 */
void HTMLWriter::write(
  ProteinMatchCollection* collection ///< collection to be written
  ) {
  if (file_ptr_ == NULL) {
    carp(CARP_FATAL, "No file open to write to.");
  }

  // psm-convert: trying to detect whether additional columns are needed.... Is there a better way to do this?
  SpectrumMatchIterator first_spec = collection->spectrumMatchBegin();
  SpectrumMatch* match = *first_spec;
  addScoreColumnIfExists(match, DELTA_CN, DELTA_CN_COL);
  addScoreColumnIfExists(match, SP, SP_SCORE_COL);
  addRankColumnIfExists(match, SP, SP_RANK_COL);
  addScoreColumnIfExists(match, XCORR, XCORR_SCORE_COL);
  addRankColumnIfExists(match, XCORR, XCORR_RANK_COL);
  addScoreColumnIfExists(match, BY_IONS_MATCHED, BY_IONS_MATCHED_COL);
  addScoreColumnIfExists(match, BY_IONS_TOTAL, BY_IONS_TOTAL_COL);

  *file_ptr_ << "<table border=\"1\">" << endl;
  writeHeader();
//  (this->*this->HTMLWriter::write_function_)(collection);
// TODO: somehow make this work for other write functions...
  writePSMs(collection);
  *file_ptr_ << "</table>" << endl;
}

void HTMLWriter::writeHeader() {

  num_columns_ = 0;
  // set file position index for all columns being printed
  for(unsigned int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++){
    if( match_to_print_[col_type] == true ){
      match_indices_[col_type] = num_columns_++;
    } else {
      match_indices_[col_type] = -1;
    }
  }

  // set all the names for which we have match_indices_
  column_names_.assign(num_columns_, "");
  for(unsigned int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++){
    if( match_indices_[col_type] > -1 ){
      if (get_column_header(col_type) == NULL) {
        carp(CARP_FATAL, "Error col type: %d doesn't exist!", col_type);
      }
      DelimitedFileWriter::setColumnName(get_column_header(col_type), 
                                         match_indices_[col_type]);
    }
  }

  if( column_names_.empty() ){
    return;
  }
  
  if( file_ptr_ == NULL || !file_ptr_->is_open() ){
    carp(CARP_FATAL, "Cannot write to NULL delimited file.");
  }
  *file_ptr_ << "\t<tr>" << endl
	     << "\t\t<td><b>" << column_names_[0] << "</b></td>" << endl;
  for(size_t idx = 1; idx < column_names_.size(); idx++){
    if( column_names_[idx].empty() ){
//      *file_ptr_ << delimiter_ << "column_" << (idx+1); 
    } else {
//      *file_ptr_ << delimiter_ << column_names_[idx];
      *file_ptr_ << "\t\t<td><b>" << column_names_[idx] << "</b></td>" << endl;
    }

  }
  
  *file_ptr_ << "\t</tr>" << endl;

  // with a header, each line must be that length
  current_row_.assign(column_names_.size(), "");
}

void HTMLWriter::write(
  MatchCollection* collection,
  string database
  ) {
  ProteinMatchCollection protein_collection(collection);
  write(&protein_collection);
}

/**
 * Write spectrum match file
 */
void HTMLWriter::writePSMs(
  ProteinMatchCollection* collection ///< collection to be written
  ) {

  carp(CARP_DEBUG, "Writing PSMs");

  string cleavage = getCleavageType();
  MASS_FORMAT_T mass_format_type =
    get_mass_format_type_parameter("mod-mass-format");

  const map<pair<int, int>, int>& spectrum_counts = collection->getMatchesSpectrum();

  bool distinct_matches = collection->hasDistinctMatches();

  for (SpectrumMatchIterator iter = collection->spectrumMatchBegin();
       iter != collection->spectrumMatchEnd();
       ++iter) {
    SpectrumMatch* match = *iter;
    setColumnCurrentRow(FILE_IDX_COL, match->getFileIndex());
    setColumnCurrentRow(FILE_COL, match->getFilePath());

    Spectrum* spectrum = match->getSpectrum();
    SpectrumZState& zstate = match->getZState();
    PeptideMatch* pep_match = match->getPeptideMatch();
    Peptide* peptide = pep_match->getPeptide();

    // psm-convert testing
    addScoreIfExists(match, DELTA_CN, DELTA_CN_COL);                          // TODO: Figure out difference between match and pep_match,
    addScoreIfExists(match, BY_IONS_MATCHED, BY_IONS_MATCHED_COL);            // since Percolator was using pep_match...
    addScoreIfExists(match, BY_IONS_TOTAL, BY_IONS_TOTAL_COL); 

    setColumnCurrentRow(SCAN_COL, spectrum->getFirstScan());
    setColumnCurrentRow(CHARGE_COL, zstate.getCharge());
    setColumnCurrentRow(SPECTRUM_PRECURSOR_MZ_COL, zstate.getMZ());
    setColumnCurrentRow(SPECTRUM_NEUTRAL_MASS_COL, zstate.getNeutralMass());
    setColumnCurrentRow(PEPTIDE_MASS_COL, peptide->getPeptideMass());
//    addScoreIfExists(pep_match, DELTA_CN, DELTA_CN_COL);                    // Original, used by percolator
    addScoreIfExists(match, SP, SP_SCORE_COL);
    addRankIfExists(match, SP, SP_RANK_COL);
    addScoreIfExists(match, XCORR, XCORR_SCORE_COL);
    addRankIfExists(match, XCORR, XCORR_RANK_COL);
    addScoreIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_SCORE_COL);
    addRankIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_RANK_COL);
    addScoreIfExists(match, PERCOLATOR_QVALUE, PERCOLATOR_QVALUE_COL);
//    addScoreIfExists(pep_match, BY_IONS_MATCHED, BY_IONS_MATCHED_COL);      // Original, used by percolator
//    addScoreIfExists(pep_match, BY_IONS_TOTAL, BY_IONS_TOTAL_COL);          // Original, used by percolator
    pair<int, int> scan_charge = make_pair(spectrum->getFirstScan(), zstate.getCharge());
    map<pair<int, int>, int>::const_iterator lookup = spectrum_counts.find(scan_charge);
    if (distinct_matches) {
      setColumnCurrentRow(DISTINCT_MATCHES_SPECTRUM_COL,
        (lookup != spectrum_counts.end()) ? lookup->second : 0);
      } else {
      setColumnCurrentRow(MATCHES_SPECTRUM_COL,
        (lookup != spectrum_counts.end()) ? lookup->second : 0);
    }
    MODIFIED_AA_T* mod_seq = peptide->getModifiedAASequence();
    char* seq_with_masses = modified_aa_string_to_string_with_masses(
      mod_seq, peptide->getLength(), mass_format_type);
    free(mod_seq);
    setAndFree(SEQUENCE_COL, seq_with_masses);

    setColumnCurrentRow(CLEAVAGE_TYPE_COL, cleavage);
//    setColumnCurrentRow(PROTEIN_ID_COL, peptide->getProteinIds()); // TODO: Figure out why this was changed ??? For percolator ???
    setColumnCurrentRow(PROTEIN_ID_COL, peptide->getProteinIdsLocations());
    setAndFree(FLANKING_AA_COL, peptide->getFlankingAAs());

    writeRow();
  }
}
