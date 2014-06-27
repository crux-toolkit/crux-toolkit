/*************************************************************************//**
 * \file SQTReader.cpp
 * \brief Object for parsing pepxml files
 ****************************************************************************/

#include "SQTReader.h"
#include "LineFileReader.h"

#include <cstdio>
#include <cstring>

#include <iostream>

#include "DelimitedFile.h"
#include "parameter.h"
#include "MatchCollectionParser.h"
#include "SQTReader.h"

using namespace std;
using namespace Crux;


const int spectrum_low_scan_idx = 1;
const int spectrum_high_scan_idx = 2;
const int spectrum_charge_idx = 3;
const int spectrum_observed_mass_idx = 6;
const int spectrum_num_matches_idx = 9;

const int match_xcorr_rank_idx = 1;
const int match_sp_rank_idx = 2;
const int match_calculated_mass_idx = 3;
const int match_delta_cn_idx = 4;
const int match_xcorr_idx = 5;
const int match_sp_idx = 6;
const int match_matched_ions_idx = 7;
const int match_expected_ions_idx = 8;
const int match_sequence_idx = 9;
const int match_validation_idx = 10;

const int locus_protein_id_idx = 1;
const int locus_protein_desc_idx = 2;

/**
 * Initializes the object
 */
void SQTReader::init() {
  last_parsed_ = SQT_LINE_NONE;
  current_spectrum_ = NULL;
  current_match_ = NULL;
  database_ = NULL;
  decoy_database_ = NULL;
}

/**
 * \returns an initialized object
 */
SQTReader::SQTReader() {
  init();
}

/**
 * \returns an object initialized with the file_path
 */
SQTReader::SQTReader(
  const string& file_path ///< the path of the pep.xml file
  ) {
  
  init();
  file_path_ = file_path;
}

/**
 * \returns an object initialized with the xml path, and the target,decoy databases
 */
SQTReader::SQTReader(
  const string& file_path, ///< the path of the pep.xml
  Database* database, ///< the protein database
  Database* decoy_database ///< the decoy protein database (can be null)
  ) {

  init();
  file_path_ = file_path;
  database_ = database;
  decoy_database_ = decoy_database;

}

/**
 * default destructor
 */
SQTReader::~SQTReader() {
}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* SQTReader::parse() {


  LineFileReader* line_reader = new LineFileReader(file_path_);
  current_match_collection_ = new MatchCollection();
  current_match_collection_->preparePostProcess();
  current_match_collection_->setScoredType(XCORR, true);
  current_match_collection_->setScoredType(SP, true);
  last_parsed_ = SQT_LINE_NONE;

  while (line_reader->hasNext()) {

    string line = line_reader->next();

    if (line.length() > 0) {

      char type = line[0];

      switch (type) {
        case 'H':
          parseHeader(line);
          break;
        case 'S':
          parseSpectrum(line);
          break;
        case 'M':
          parseMatch(line);
          break;
        case 'L':
          parseLocus(line);
          break;
        default:
          carp(CARP_ERROR, "Unknown line %d\n%s",line_reader->getCurrentRow(), line.c_str());
      }
    } else {
      carp(CARP_ERROR, "blank line at row %d", line_reader->getCurrentRow());
    }
  }

  delete line_reader;

  //add the last match
  if (current_match_ != NULL) {
    current_match_collection_->addMatchToPostMatchCollection(current_match_);
  }

  return current_match_collection_;
}


void SQTReader::parseHeader(string& line) {
  last_parsed_ = SQT_LINE_HEADER;
  headers_.push_back(line);
}


void SQTReader::parseSpectrum(string& line) {

  vector<string> tokens;

  tokenize(line, tokens, '\t');

  int low_scan;
  from_string(low_scan, tokens[spectrum_low_scan_idx]);
  
  int high_scan;
  from_string(high_scan, tokens[spectrum_high_scan_idx]);

  int charge;
  from_string(charge, tokens[spectrum_charge_idx]);

  double observed_mass;
  from_string(observed_mass, tokens[spectrum_observed_mass_idx]);

  from_string(current_num_matches_, tokens[spectrum_num_matches_idx]);

  current_ln_experiment_size_ = logf((FLOAT_T)current_num_matches_);
  
  last_parsed_ = SQT_LINE_SPECTRUM;

  current_zstate_.setSinglyChargedMass(observed_mass, charge);
  current_spectrum_ = new Spectrum(
    low_scan,
    high_scan,
    current_zstate_.getMZ(),
    vector<int>(1, charge),
    "");
  /*
  cerr << "spectrum line:"<<line<<endl;
  cerr << "low scan:"<<low_scan<<endl;
  cerr << "high scan:"<<high_scan<<endl;
  cerr << "charge:"<<charge<<endl;
  cerr << "observed mass:"<<observed_mass<<endl;
  cerr << "num matches:"<<num_matches<<endl;
  cerr << "======================="<<endl;
 */

}


void SQTReader::parseMatch(string& line) {

  vector<string> tokens;
  tokenize(line, tokens, '\t');

  int xcorr_rank;
  from_string(xcorr_rank, tokens[match_xcorr_rank_idx]);
  
  int sp_rank;
  from_string(sp_rank, tokens[match_sp_rank_idx]);

  double calculated_mass;
  from_string(calculated_mass, tokens[match_calculated_mass_idx]);

  double delta_cn;
  from_string(delta_cn, tokens[match_delta_cn_idx]);

  double xcorr;
  from_string(xcorr, tokens[match_xcorr_idx]);

  double sp;
  from_string(sp, tokens[match_sp_idx]);

  int matched_ions;
  from_string(matched_ions, tokens[match_matched_ions_idx]);

  int expected_ions;
  from_string(expected_ions, tokens[match_expected_ions_idx]);

  string sqt_sequence = tokens[match_sequence_idx];
  vector<string> sequence_tokens;
  tokenize(sqt_sequence, sequence_tokens, '.');

  current_prev_aa_ = sequence_tokens.front();
  current_next_aa_ = sequence_tokens.back();
  sequence_tokens.erase(sequence_tokens.begin());
  sequence_tokens.pop_back();
  current_peptide_sequence_ = DelimitedFile::splice(sequence_tokens, '.');
 /*
  cerr << "Match line:"<<line<<endl;
  cerr << "xcorr rank:"<<xcorr_rank<<endl;
  cerr << "sp rank:"<<sp_rank<<endl;
  cerr << "calculated mass:"<<calculated_mass<<endl;
  cerr << "delta cn:"<<delta_cn<<endl;
  cerr << "xcorr:"<<xcorr<<endl;
  cerr << "sp:"<<sp<<endl;
  cerr << "matched ions:"<<matched_ions<<endl;
  cerr << "expected ions:"<<expected_ions<<endl;
  cerr << "sequence:"<<current_peptide_sequence_<<endl;
  cerr << "prev_aa:"<<current_prev_aa_<<endl;
  cerr << "next_aa:"<<current_next_aa_<<endl;
  cerr << "====================" << endl;
*/
  if (current_match_ != NULL) {
    current_match_collection_->addMatchToPostMatchCollection(current_match_);
  }

  Peptide* peptide = new Peptide();

  peptide->setPeptideMass(calculated_mass - MASS_PROTON);

  MODIFIED_AA_T* mods = NULL;
  convert_to_mod_aa_seq(current_peptide_sequence_.c_str(), &mods);

  char* unmodified_sequence_pointer = unmodify_sequence(current_peptide_sequence_.c_str());
  current_peptide_sequence_ = unmodified_sequence_pointer;
  free(unmodified_sequence_pointer);
  peptide->setLength(current_peptide_sequence_.length());
  peptide->setModifiedAASequence(mods, false);

  current_match_ = new Match(peptide, current_spectrum_, current_zstate_, false);
  current_match_->setScore(XCORR, xcorr);
  current_match_->setScore(SP, sp);
  current_match_->setDeltaCn(delta_cn);
  current_match_->setRank(XCORR, xcorr_rank);
  current_match_->setRank(SP, sp_rank);

  current_match_->setBYIonMatched(matched_ions);
  current_match_->setBYIonPossible(expected_ions);
  current_match_->setBYIonFractionMatched((FLOAT_T)matched_ions / (FLOAT_T)expected_ions);
  current_match_->setLnExperimentSize(current_ln_experiment_size_);

  last_parsed_ = SQT_LINE_MATCH;
}

void SQTReader::parseLocus(string& line) {

  vector<string> tokens;
  tokenize(line, tokens, '\t');

  string protein_id = tokens[locus_protein_id_idx];
  string protein_desc = "";
  if (tokens.size() > 2) {
    protein_desc = tokens[locus_protein_desc_idx];
  }
/*
  cerr << "Locus line:"<<line<<endl;
  cerr << "Protein id:"<<protein_id<<endl;
  cerr << "Protein desc:"<<protein_desc<<endl;
  cerr << "========================="<<endl;
*/
  bool is_decoy;
  Protein* protein = MatchCollectionParser::getProtein(database_,
                                                       decoy_database_,
                                                       protein_id,
                                                       is_decoy);

  int start_idx = protein->findStart(current_peptide_sequence_, current_prev_aa_, current_next_aa_);

  PeptideSrc* peptide_src = new PeptideSrc((DIGEST_T)0, protein, start_idx);
  
  current_match_->getPeptide()->addPeptideSrc(peptide_src);

  if (is_decoy) {
    current_match_->setNullPeptide(true);
  }

  last_parsed_ = SQT_LINE_LOCUS;
}

/**
 * /returns the start position of the peptide sequence within the protein
 */
int SQTReader::findStart(
  Protein* protein,  ///< the protein to find the sequence 
  string peptide_sequence, ///< the peptide sequence to find
  string prev_aa, ///< the amino acid before the sequence in the protein
  string next_aa ///< the next amino acid after the sequence in the protein
  ) {

  if (prev_aa == "-") {
    return 1;
  } else if (next_aa == "-") {
    return protein->getLength() - peptide_sequence.length() + 1;
  } else {
    //use the flanking amino acids to further constrain our search in the sequence
    size_t pos = string::npos; 
    string seq = prev_aa + peptide_sequence + next_aa;
    string protein_seq = protein->getSequencePointer();
    pos = protein_seq.find(seq);
    if (pos == string::npos) {
      carp(CARP_DEBUG, "could not find %s in protein %s\n%s", seq.c_str(), protein->getIdPointer(), protein_seq.c_str());
      //finding the sequence with the flanks failed, try finding without the flanks.
      seq = peptide_sequence;
      pos = protein_seq.find(seq);
      if (pos == string::npos) {
        carp(CARP_FATAL, "could not %s in protein %s\n%s", seq.c_str(), protein->getIdPointer(), protein_seq.c_str());
      }
      return (pos+1);
    }
    return (pos+2);
  }
}

/**
 * sets the target protein database
 */
void SQTReader::setDatabase(
  Database* database ///< the target protein database
  ) {

  database_ = database;
}

/**
 * sets the decoy protein database
 */
void SQTReader::setDecoyDatabase(
  Database* decoy_database ///< sets the decoy protein database
  ) {
  decoy_database_ = decoy_database;
}


/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* SQTReader::parse(
  const char* file_path, ///< path of the xml file
  Database* database, ///< target protein database
  Database* decoy_database ///< decoy protein database (can be null)
  ) {
  SQTReader* reader = new SQTReader(file_path);
  reader->setDatabase(database);
  reader->setDecoyDatabase(decoy_database);
  MatchCollection* collection = reader->parse();

  return collection;


}

//TODO - remove this code after some time of debugging.
#ifdef MAIN
int main(int argc, char** argv) {

  initialize_parameters();

  char* file_path = argv[1];
  char* database_path = argv[2];

  MatchCollection* match_collection = MatchCollectionParser::create(file_path, database_path);

  cerr << "there are "<<match_collection->getMatchTotal()<<" matches read"<<endl;

  MatchIterator* match_iterator = new MatchIterator(match_collection, XCORR, true);

  while(match_iterator->hasNext()) {
    Match* match = match_iterator->next();

    cout << "xcorr:"<<match->getScore(XCORR);
    cout <<" rank:"<<match->getRank(XCORR);
    cout <<" sequence:"<<match->getPeptide()->getSequence();
    cout <<" protein:"<< match->getPeptide()->getProteinIdsLocations()<<endl;

  }



  return 0;
}
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
