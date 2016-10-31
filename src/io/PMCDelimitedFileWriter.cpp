#include "PMCDelimitedFileWriter.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"

using namespace Crux;
using namespace boost;

/** 
 * Returns an empty PMCDelimitedFileWriter object
 */
PMCDelimitedFileWriter::PMCDelimitedFileWriter() : PSMWriter() {
  write_html_ = false;
  application_ = NULL;
}

/**
 * Destructor
 */
PMCDelimitedFileWriter::~PMCDelimitedFileWriter() {
  closeFile();
}

/**
 * Opens a file, writes it, and closes the file
 */
void PMCDelimitedFileWriter::writeFile(
  CruxApplication* application, ///< application writing the file
  string filename, ///< name of the file to open
  MATCH_FILE_TYPE type, ///< type of file to be written
  ProteinMatchCollection* collection ///< collection to be written
) {
  openFile(application, filename, type);
  write(collection);
  closeFile();
}

/**
 * Closes any open file, then opens a file for the specified type of writing
 */
void PMCDelimitedFileWriter::openFile(
  CruxApplication* application, ///< application writing the file
  string filename, ///< name of the file to open
  MATCH_FILE_TYPE type ///< type of file to be written
  ) {
  closeFile();

  file_ptr_ = FileUtils::GetWriteStream(filename, Params::GetBool("overwrite"));
  application_ = application;
  if (!file_ptr_->is_open()) {
    carp(CARP_FATAL, "Error creating file '%s'.", filename.c_str());
  }

  // reset columns
  num_columns_ = 0;
  for (int i = 0; i < NUMBER_MATCH_COLUMNS; ++i) {
    match_to_print_[i] = false;
    match_precision_[i] = false;
    match_fixed_float_[i] = true;
  }

  switch (type) {
  case PROTEINS:
    setUpProteinsColumns(application);
    setPrecision();
    write_function_ = &PMCDelimitedFileWriter::writeProteins;
    break;
  case PEPTIDES:
    setUpPeptidesColumns(application);
    setPrecision();
    write_function_ = &PMCDelimitedFileWriter::writePeptides;
    break;
  case PSMS:
    setUpPSMsColumns(application);
    setPrecision();
    write_function_ = &PMCDelimitedFileWriter::writePSMs;
    break;
  default:
    carp(CARP_FATAL, "Invalid match file type specified for "
                     "PMCDelimitedFileWriter.");
    break;
  }
}

/**
 * Closes any open file, if any
 */
void PMCDelimitedFileWriter::closeFile() {
  if (file_ptr_) {
    file_ptr_->close();
    delete file_ptr_;
    file_ptr_ = NULL;
  }
}

/**
 * Write all three file types using the file stem as a base filename
 */
void PMCDelimitedFileWriter::writeAll(
  CruxApplication* application, ///< application writing the files
  ProteinMatchCollection* collection, ///< collection to be written
  string stem ///< filestem to be prepended to the filenames
  ) {

  // write proteins
  openFile(application, stem + ".proteins.txt", PROTEINS);
  write(collection);
  closeFile();

  // write peptides
  openFile(application, stem + ".peptides.txt", PEPTIDES);
  write(collection);
  closeFile();

  // write psms
  openFile(application, stem + ".psms.txt", PSMS);
  write(collection);
  closeFile();

}

void PMCDelimitedFileWriter::write(
  MatchCollection* collection,
  string database
  ) {
  ProteinMatchCollection protein_collection(collection);
  write(&protein_collection);
}

/**
 * Writes the data in a ProteinMatchCollection to the currently open file
 */
void PMCDelimitedFileWriter::write(
  ProteinMatchCollection* collection ///< collection to be written
  ) {
  if (file_ptr_ == NULL) {
    carp(CARP_FATAL, "No file open to write to.");
  } else if (collection == NULL) {
    carp(CARP_FATAL, "ProteinMatchCollection was null");
  }

  // psm-convert: trying to detect whether additional columns are needed....
  // Is there a better way to do this or move to setUpPSMSColumns? Seems messy...
  SpectrumMatchIterator first_spec = collection->spectrumMatchBegin();
  SpectrumMatch* match = *first_spec;
  
  if (!match->getFilePath().empty()) {
    addColumnName(FILE_COL);
  }
  addScoreColumnIfExists(match, DELTA_CN, DELTA_CN_COL);
  addScoreColumnIfExists(match, SP, SP_SCORE_COL);
  addRankColumnIfExists(match, SP, SP_RANK_COL);
  addScoreColumnIfExists(match, XCORR, XCORR_SCORE_COL);
  addRankColumnIfExists(match, XCORR, XCORR_RANK_COL);
  addScoreColumnIfExists(match, TIDE_SEARCH_EXACT_PVAL, EXACT_PVALUE_COL);
  addScoreColumnIfExists(match, TIDE_SEARCH_REFACTORED_XCORR, REFACTORED_SCORE_COL);
  addScoreColumnIfExists(match, BY_IONS_MATCHED, BY_IONS_MATCHED_COL);
  addScoreColumnIfExists(match, BY_IONS_TOTAL, BY_IONS_TOTAL_COL);
  if (collection->hasDistinctMatches()) {
    addColumnName(DISTINCT_MATCHES_SPECTRUM_COL);
  } else {
    addColumnName(MATCHES_SPECTRUM_COL);
  }

  if (!write_html_) {
    writeHeader();
  } else {
    writeHTMLHeader();
  }
  (this->*this->PMCDelimitedFileWriter::write_function_)(collection);
}

// adds score to header and columns if exists
void PMCDelimitedFileWriter::addScoreColumnIfExists(
  AbstractMatch* match, ///< match to get score from
  SCORER_TYPE_T scoreType, ///< score type to get
  MATCH_COLUMNS_T column ///< column to add the score to
  ) {
  if (match->hasScore(scoreType)) {
    addColumnName(column);
  }
}

// adds rank to header and columns if exists
void PMCDelimitedFileWriter::addRankColumnIfExists(
  AbstractMatch* match, ///< match to get rank from
  SCORER_TYPE_T scoreType, ///< rank type to get
  MATCH_COLUMNS_T column ///< column to add the rank to
  ) {
  if (match->hasRank(scoreType)) {
    addColumnName(column);
  }
}

/**
 * Set up the columns for a protein match file
 */
void PMCDelimitedFileWriter::setUpProteinsColumns(
  CruxApplication* application ///< application writing the file
) {

  COMMAND_T command = application->getCommand();

  switch(command) {
  default:
    carp(CARP_FATAL, "Command (%s) not yet implemented for writing tab "
                     "delimited files.", application->getName().c_str());
    return;
  case MISC_COMMAND:
  case PROCESS_SPEC_COMMAND:
  case GENERATE_PEPTIDES_COMMAND:
  case GET_MS2_SPECTRUM_COMMAND:
  case PREDICT_PEPTIDE_IONS_COMMAND:
  case VERSION_COMMAND:
  case NUMBER_COMMAND_TYPES:
  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command (%s) for writing tab delimited file.",
                     application->getName().c_str());
    return;
  case PERCOLATOR_COMMAND:
    addColumnName(PERCOLATOR_SCORE_COL);
    addColumnName(PERCOLATOR_RANK_COL);
    addColumnName(PERCOLATOR_QVALUE_COL);
    break;
  }

  addColumnName(PROTEIN_ID_COL);
  addColumnName(PEPTIDES_COL);
}

/**
 * Write protein match file
 */
void PMCDelimitedFileWriter::writeProteins(
  ProteinMatchCollection* collection ///< collection to be written
  ) {

  carp(CARP_DEBUG, "Writing proteins");

  MASS_FORMAT_T mass_format_type =
    get_mass_format_type_parameter("mod-mass-format");

  for (ProteinMatchIterator iter = collection->proteinMatchBegin();
       iter != collection->proteinMatchEnd();
       ++iter) {
    ProteinMatch* match = *iter;
    Protein* protein = match->getProtein();

    addScoreIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_SCORE_COL);
    addRankIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_RANK_COL);
    addScoreIfExists(match, PERCOLATOR_QVALUE, PERCOLATOR_QVALUE_COL);
    setColumnCurrentRow(PROTEIN_ID_COL, protein->getIdPointer());

    // build peptide list
    vector<string> peptide_strings;
    for (PeptideMatchIterator pep_iter = match->peptideMatchBegin();
         pep_iter != match->peptideMatchEnd();
         ++pep_iter) {
      PeptideMatch* pep_match = *pep_iter;
      Peptide* peptide = pep_match->getPeptide();

      MODIFIED_AA_T* mod_seq = peptide->getModifiedAASequence();
      char* seq_with_masses = modified_aa_string_to_string_with_masses(
        mod_seq, peptide->getLength(), mass_format_type);
      string sequence_str(seq_with_masses);
      free(seq_with_masses);
      free(mod_seq);

      for (SpectrumMatchIterator spec_iter = pep_match->spectrumMatchBegin();
           spec_iter != pep_match->spectrumMatchEnd();
           ++spec_iter) {
        SpectrumMatch* spec_match = *spec_iter;
        Spectrum* spectrum = spec_match->getSpectrum();
        SpectrumZState& zstate = spec_match->getZState();

        stringstream ss;
        ss << sequence_str << '-' << spectrum->getFirstScan() << '.'
           << zstate.getCharge();

        peptide_strings.push_back(ss.str());
      }
    }
    setColumnCurrentRow(PEPTIDES_COL, StringUtils::Join(peptide_strings, ','));

    writeRow();
  }
}

/**
 * Set up the columns for a peptide match file
 */
void PMCDelimitedFileWriter::setUpPeptidesColumns(
  CruxApplication* application ///< application writing the file
) {

  COMMAND_T command = application->getCommand();

  switch(command) {
  default:
    carp(CARP_FATAL, "Command (%s) not yet implemented for writing tab "
                     "delimited files.", application->getName().c_str());
    return;
  case MISC_COMMAND:
  case PROCESS_SPEC_COMMAND:
  case GENERATE_PEPTIDES_COMMAND:
  case GET_MS2_SPECTRUM_COMMAND:
  case PREDICT_PEPTIDE_IONS_COMMAND:
  case VERSION_COMMAND:
  case NUMBER_COMMAND_TYPES:
  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command (%s) for writing tab delimited file.",
                     application->getName().c_str());
    return;
  case PERCOLATOR_COMMAND:
    addColumnName(PERCOLATOR_SCORE_COL);
    addColumnName(PERCOLATOR_RANK_COL);
    addColumnName(PERCOLATOR_QVALUE_COL);
    addColumnName(PERCOLATOR_PEP_COL);
    break;
  }

  addColumnName(SCAN_COL);
  addColumnName(CHARGE_COL);
  addColumnName(SPECTRUM_PRECURSOR_MZ_COL);
  addColumnName(SPECTRUM_NEUTRAL_MASS_COL);
  addColumnName(PEPTIDE_MASS_COL);
  //addColumnName(DELTA_CN_COL);
  //addColumnName(SP_SCORE_COL);
  //addColumnName(SP_RANK_COL);
  //addColumnName(XCORR_SCORE_COL);
  //addColumnName(XCORR_RANK_COL);
  //addColumnName(BY_IONS_MATCHED_COL);
  //addColumnName(BY_IONS_TOTAL_COL);
  //addColumnName(MATCHES_SPECTRUM_COL);
  addColumnName(SEQUENCE_COL);
  addColumnName(CLEAVAGE_TYPE_COL);
  addColumnName(PROTEIN_ID_COL);
  addColumnName(FLANKING_AA_COL);
}

/**
 * Write peptide match file
 */
void PMCDelimitedFileWriter::writePeptides(
  ProteinMatchCollection* collection ///< collection to be written
  ) {

  carp(CARP_DEBUG, "Writing peptides");

  string cleavage = getCleavageType();
  MASS_FORMAT_T mass_format_type =
    get_mass_format_type_parameter("mod-mass-format");

  for (PeptideMatchIterator iter = collection->peptideMatchBegin();
       iter != collection->peptideMatchEnd();
       ++iter) {
    PeptideMatch* match = *iter;
    Peptide* peptide = match->getPeptide();

    // collect spectra info
    vector<string> spec_scans;
    vector<string> spec_charges;
    vector<string> spec_precursors;
    vector<string> spec_neutral_masses;
    for (SpectrumMatchIterator spec_iter = match->spectrumMatchBegin();
         spec_iter != match->spectrumMatchEnd();
         ++spec_iter) {
      SpectrumMatch* spec_match = *spec_iter;
      Spectrum* spec = spec_match->getSpectrum();
      SpectrumZState& zstate = spec_match->getZState();

      spec_scans.push_back(lexical_cast<string>(spec->getFirstScan()));
      spec_charges.push_back(lexical_cast<string>(zstate.getCharge()));
      spec_precursors.push_back(lexical_cast<string>(spec->getPrecursorMz()));
      spec_neutral_masses.push_back(lexical_cast<string>(zstate.getNeutralMass()));
    }

    // collect proteins info
    vector<string> protein_ids;
    vector<string> flanking_aas;
    for (ProteinMatchIterator prot_iter = match->proteinMatchBegin();
         prot_iter != match->proteinMatchEnd();
         ++prot_iter) {
      ProteinMatch* prot_match = *prot_iter;
      Protein* protein = prot_match->getProtein();
      protein_ids.push_back(protein->getIdPointer());
      /** TODO figure out why this works and getFlankingAAs() doesn't **/
      if (protein->isPostProcess()) {
        PostProcessProtein* post_process_protein = (PostProcessProtein*)protein;
        string flanking_str = "";
        flanking_str += post_process_protein->getNTermFlankingAA();
        flanking_str += post_process_protein->getCTermFlankingAA();
        flanking_aas.push_back(flanking_str);
      }
      /*****************************************************************/
    }

    setColumnCurrentRow(SCAN_COL, StringUtils::Join(spec_scans, ','));
    setColumnCurrentRow(CHARGE_COL, StringUtils::Join(spec_charges, ','));
    setColumnCurrentRow(SPECTRUM_PRECURSOR_MZ_COL, StringUtils::Join(spec_precursors, ','));
    setColumnCurrentRow(SPECTRUM_NEUTRAL_MASS_COL, StringUtils::Join(spec_neutral_masses, ','));
    setColumnCurrentRow(PEPTIDE_MASS_COL, peptide->calcModifiedMass());
    addScoreIfExists(match, DELTA_CN, DELTA_CN_COL);
    addScoreIfExists(match, SP, SP_SCORE_COL);
    addRankIfExists(match, SP, SP_RANK_COL);
    addScoreIfExists(match, XCORR, XCORR_SCORE_COL);
    addRankIfExists(match, XCORR, XCORR_RANK_COL);
    addScoreIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_SCORE_COL);
    addRankIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_RANK_COL);
    addScoreIfExists(match, PERCOLATOR_QVALUE, PERCOLATOR_QVALUE_COL);
    addScoreIfExists(match, PERCOLATOR_PEP, PERCOLATOR_PEP_COL);
    addScoreIfExists(match, BY_IONS_MATCHED, BY_IONS_MATCHED_COL);
    addScoreIfExists(match, BY_IONS_TOTAL, BY_IONS_TOTAL_COL);
    //addScoreIfExists(match, MATCHES_SPECTRUM, MATCHES_SPECTRUM_COL);

    MODIFIED_AA_T* mod_seq = peptide->getModifiedAASequence();
    char* seq_with_masses = modified_aa_string_to_string_with_masses(
      mod_seq, peptide->getLength(), mass_format_type);
    free(mod_seq);
    setAndFree(SEQUENCE_COL, seq_with_masses);

    setColumnCurrentRow(CLEAVAGE_TYPE_COL, cleavage);
    setColumnCurrentRow(PROTEIN_ID_COL, StringUtils::Join(protein_ids, ','));
    //setAndFree(FLANKING_AA_COL, peptide->getFlankingAAs());
    setColumnCurrentRow(FLANKING_AA_COL, StringUtils::Join(flanking_aas, ','));

    writeRow();
  }
}

/**
 * Set up the columns for a spectrum match file
 */
void PMCDelimitedFileWriter::setUpPSMsColumns(
  CruxApplication* application ///< application writing the file
) {
  switch (application->getCommand()) {
  default:
    carp(CARP_FATAL, "Command (%s) not yet implemented for writing tab "
                     "delimited files.", application->getName().c_str());
    return;
  case MISC_COMMAND:
  case PROCESS_SPEC_COMMAND:
  case GENERATE_PEPTIDES_COMMAND:
  case GET_MS2_SPECTRUM_COMMAND:
  case PREDICT_PEPTIDE_IONS_COMMAND:
  case VERSION_COMMAND:
  case NUMBER_COMMAND_TYPES:
  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command (%s) for writing tab delimited file.",
                     application->getName().c_str());
    return;
  case PSM_CONVERT_COMMAND:
    break;
  case PERCOLATOR_COMMAND:
    addColumnName(FILE_IDX_COL);
    addColumnName(FILE_COL);
    addColumnName(PERCOLATOR_SCORE_COL);
    addColumnName(PERCOLATOR_RANK_COL);
    addColumnName(PERCOLATOR_PEP_COL);
    addColumnName(PERCOLATOR_QVALUE_COL);
    break;
  }
  addColumnName(SCAN_COL);
  addColumnName(CHARGE_COL);
  addColumnName(SPECTRUM_PRECURSOR_MZ_COL);
  addColumnName(SPECTRUM_NEUTRAL_MASS_COL);
  addColumnName(PEPTIDE_MASS_COL);
  addColumnName(SEQUENCE_COL);
  addColumnName(MODIFICATIONS_COL);
  addColumnName(CLEAVAGE_TYPE_COL);
  addColumnName(PROTEIN_ID_COL);
  addColumnName(FLANKING_AA_COL);
}

/**
 * Write spectrum match file
 */
void PMCDelimitedFileWriter::writePSMs(
  ProteinMatchCollection* collection ///< collection to be written
  ) {

  carp(CARP_DEBUG, "Writing PSMs");

  string cleavage = getCleavageType();
  MASS_FORMAT_T mass_format_type =
    get_mass_format_type_parameter("mod-mass-format");

  const map<pair<int, int>, int>& spectrum_counts = collection->getMatchesSpectrum();

  bool distinct_matches = collection->hasDistinctMatches();
  COMMAND_T command = application_->getCommand();

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

    if (command == PSM_CONVERT_COMMAND) {
      addScoreIfExists(match, DELTA_CN, DELTA_CN_COL);                          // TODO: Figure out difference between match and pep_match,
      addScoreIfExists(match, BY_IONS_MATCHED, BY_IONS_MATCHED_COL);            // since Percolator uses pep_match... While tide-search
      addScoreIfExists(match, BY_IONS_TOTAL, BY_IONS_TOTAL_COL);                // uses match. These are different values, are they supposed to be?
                    // When using pep_match, a lot of things are unknown (-1)
      setColumnCurrentRow(PROTEIN_ID_COL, peptide->getProteinIdsLocations());   // This was used by an older PMCDelimitedFileWriter and retains the protein id location
    } else {
      addScoreIfExists(pep_match, DELTA_CN, DELTA_CN_COL);
      addScoreIfExists(pep_match, BY_IONS_MATCHED, BY_IONS_MATCHED_COL);
      addScoreIfExists(pep_match, BY_IONS_TOTAL, BY_IONS_TOTAL_COL);
      setColumnCurrentRow(PROTEIN_ID_COL, StringUtils::Join(peptide->getProteinIds(), ','));
    }

    setColumnCurrentRow(SCAN_COL, spectrum->getFirstScan());
    setColumnCurrentRow(CHARGE_COL, zstate.getCharge());
    setColumnCurrentRow(SPECTRUM_PRECURSOR_MZ_COL, zstate.getMZ());
    setColumnCurrentRow(SPECTRUM_NEUTRAL_MASS_COL, zstate.getNeutralMass());
    setColumnCurrentRow(PEPTIDE_MASS_COL, peptide->calcModifiedMass());
    addScoreIfExists(match, SP, SP_SCORE_COL);
    addRankIfExists(match, SP, SP_RANK_COL);
    addScoreIfExists(match, XCORR, XCORR_SCORE_COL);
    addRankIfExists(match, XCORR, XCORR_RANK_COL);
    addScoreIfExists(match, TIDE_SEARCH_EXACT_PVAL, EXACT_PVALUE_COL);
    addScoreIfExists(match, TIDE_SEARCH_REFACTORED_XCORR, REFACTORED_SCORE_COL);
    addScoreIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_SCORE_COL);
    addRankIfExists(match, PERCOLATOR_SCORE, PERCOLATOR_RANK_COL);
    addScoreIfExists(match, PERCOLATOR_PEP, PERCOLATOR_PEP_COL);
    addScoreIfExists(match, PERCOLATOR_QVALUE, PERCOLATOR_QVALUE_COL);
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
    setColumnCurrentRow(MODIFICATIONS_COL, peptide->getModsString());

    setColumnCurrentRow(CLEAVAGE_TYPE_COL, cleavage);
    setAndFree(FLANKING_AA_COL, peptide->getFlankingAAs());

    if (!write_html_) {
      writeRow();
    } else {
      writeHTMLRow();
    }

  }
  if (write_html_) {
    *file_ptr_ << "</table>" << endl;
  }
}

/**
 * Gets the cleavage type as a string
 */
string PMCDelimitedFileWriter::getCleavageType() {
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  const char* enzyme_string = enzyme_type_to_string(enzyme);
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  const char* digestion_string = digest_type_to_string(digestion);
  string cleavage_str = enzyme_string;
  cleavage_str += "-";
  cleavage_str += digestion_string;
  return cleavage_str;
}

/**
 * Adds the value of a score to the specified column, if it exists in the match
 */
void PMCDelimitedFileWriter::addScoreIfExists(
  AbstractMatch* match, ///< match to get score from
  SCORER_TYPE_T scoreType, ///< score type to get
  MATCH_COLUMNS_T column ///< column to add the score to
  ) {
  if (match->hasScore(scoreType)) {
    setColumnCurrentRow(column, match->getScore(scoreType));
  } else {
    setColumnCurrentRow(column, -1);
  }
}

/**
 * Adds the value of a rank to the specified column, if it exists in the match
 */
void PMCDelimitedFileWriter::addRankIfExists(
  AbstractMatch* match, ///< match to get rank from
  SCORER_TYPE_T scoreType, ///< rank type to get
  MATCH_COLUMNS_T column ///< column to add the rank to
  ) {
  if (match->hasRank(scoreType)) {
    setColumnCurrentRow(column, match->getRank(scoreType));
  } else {
    setColumnCurrentRow(column, -1);
  }
}

/**
 * Sets the value of a char* to the specified column, then frees the pointer
 */
void PMCDelimitedFileWriter::setAndFree(
  MATCH_COLUMNS_T column, ///< column to add the char* to
  char* value ///< char* to set the column value to
  ) {
  if (value != NULL) {
    setColumnCurrentRow(column, value);
    free(value);
  } else {
    carp(CARP_WARNING, "Cannot set value for column %d", column);
  }
}

/**
  * Sets whether we should write our output in tab delimited or HTML format. Default = false (tab delimited)
  */
void PMCDelimitedFileWriter::setWriteHTML(bool write_html) {
  write_html_ = write_html;
}

/**
 * Writes the header in HTML format, this is used for HTMLWriter instead of
 * DelimitedFileWriter's writeHeader()
 */
void PMCDelimitedFileWriter::writeHTMLHeader() {

  num_columns_ = 0;
  // set file position index for all columns being printed
  for(unsigned int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++) {
    if( match_to_print_[col_type] == true ) {
      match_indices_[col_type] = num_columns_++;
    } else {
      match_indices_[col_type] = -1;
    }
  }

  // set all the names for which we have match_indices_
  column_names_.assign(num_columns_, "");
  for(unsigned int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++) {
    if( match_indices_[col_type] > -1 ) {
      if (get_column_header(col_type) == NULL) {
        carp(CARP_FATAL, "Error col type: %d doesn't exist!", col_type);
      }
      DelimitedFileWriter::setColumnName(get_column_header(col_type), 
                                         match_indices_[col_type]);
    }
  }

  if( column_names_.empty() ) {
    return;
  }
  
  if( file_ptr_ == NULL || !file_ptr_->is_open() ) {
    carp(CARP_FATAL, "Cannot write to NULL delimited file.");
  }
  *file_ptr_ << "<table border=\"1\">" << "\t<tr>" << endl
      << "\t\t<td><b>" << column_names_[0] << "</b></td>" << endl;
  for(size_t idx = 1; idx < column_names_.size(); idx++) {
    if( !column_names_[idx].empty() ) {
      *file_ptr_ << "\t\t<td><b>" << column_names_[idx] << "</b></td>" << endl;
    }
  }
  
  *file_ptr_ << "\t</tr>" << endl;

  // with a header, each line must be that length
  current_row_.assign(column_names_.size(), "");
}

/**
 * Writes a row in HTML format, this is used for HTMLWriter instead of
 * DelimitedFileWriter's writeRow() 
 */
void PMCDelimitedFileWriter::writeHTMLRow() {
  if( current_row_.empty() ) {
    return;
  }

  // make the row as long as the header
  while( current_row_.size() < column_names_.size()) {
    current_row_.push_back("");
  }
  // TODO? warning if row is longer than non-empty header?

  // print each value separated by delimiter
  *file_ptr_ << "\t<tr>" << endl
             << "\t\t<td>" << current_row_[0] << "</td>" << endl;
  for(size_t idx = 1; idx < current_row_.size(); idx++) {
    *file_ptr_ << "\t\t<td>" << current_row_[idx] << "</td>" << endl;
  }
  // end with newline
  *file_ptr_ << "\t</tr>" << endl;

  // clear the current_row and refill with blanks
  // if there is a header, that is the min length
  // if not, each row can be a different length
  current_row_.assign(column_names_.size(), "");

}
