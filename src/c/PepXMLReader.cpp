/*************************************************************************//**
 * \file PepXMLReader.cpp
 * \brief Object for parsing pepxml files
 ****************************************************************************/

#include "PepXMLReader.h"
#include "mass.h"
#include "expat.h"

#include "Protein.h"
#include "PostProcessProtein.h"
#include "Peptide.h"


#include <cstdio>
#include <cstring>

#include <iostream>

#include "DelimitedFile.h"
#include "parameter.h"
#include "MatchCollectionParser.h"


using namespace std;
using namespace Crux;

//Buffer for expat's xml reading routines
#define BUFFSIZE 8192
char Buff[BUFFSIZE];

void open_handler(void *data, const char *el, const char **attr) {
  
  PepXMLReader* reader = (PepXMLReader*)data;
  if (strcmp(el, "spectrum_query") == 0) {
    reader->spectrumQueryOpen(attr);
  } else if (strcmp(el, "search_result") == 0) {
    reader->searchResultOpen();
  } else if (strcmp(el, "search_hit") == 0) {
    reader->searchHitOpen(attr);
  } else if (strcmp(el, "modification_info") == 0) {
    reader->modificationInfoOpen(attr);
  } else if (strcmp(el, "mod_aminoacid_mass") == 0) {
    reader->modAminoAcidMassOpen(attr);
  } else if (strcmp(el, "alternative_protein") == 0) {
    reader->alternativeProteinOpen(attr);
  } else if (strcmp(el, "search_score") == 0) {
    reader->searchScoreOpen(attr);
  } else if (strcmp(el, "peptideprophet_result") == 0) {
    reader->peptideProphetResultOpen(attr);
  } else {
    carp(CARP_DEBUG, "Unsupported open tag:%s", el);
  }
}

void close_handler(void *data, const char *el) {

  PepXMLReader* reader = (PepXMLReader*)data;
  if (strcmp(el, "spectrum_query") == 0) {
    reader->spectrumQueryClose();
  } else if (strcmp(el, "search_result") == 0) {
    reader->searchResultClose();
  } else if (strcmp(el, "search_hit") == 0) {
    reader->searchHitClose();
  } else if (strcmp(el, "modification_info") == 0) {
    reader->modificationInfoClose();
  } else if (strcmp(el, "mod_aminoacid_mass") == 0) {
    reader->modAminoAcidMassClose();
  } else if (strcmp(el, "alternative_protein") == 0) {
    reader->alternativeProteinClose();
  } else if (strcmp(el, "search_score") == 0) {
    reader->searchScoreClose();
  } else if (strcmp(el, "peptideprophet_result") == 0) {
    reader->peptideProphetResultClose();
  } else {
    carp(CARP_DEBUG, "Unsupported close tag:%s", el);
  }
}  /* End of end handler */

/**
 * Initializes the object
 */
void PepXMLReader::init() {
  spectrum_query_open_ = false;
  search_result_open_ = false;
  search_hit_open_ = false;
  alternative_protein_open_ = false;
  search_score_open_ = false;
  peptideprophet_result_open_ = false;

  current_spectrum_ = NULL;
  database_ = NULL;
  decoy_database_ = NULL;
}

/**
 * \returns an initialized object
 */
PepXMLReader::PepXMLReader() {
  init();
}

/**
 * \returns an object initialized with the file_path
 */
PepXMLReader::PepXMLReader(
  const string& file_path ///< the path of the pep.xml file
  ) {
  
  init();
  file_path_ = file_path;
}

/**
 * \returns an object initialized with the xml path, and the target,decoy databases
 */
PepXMLReader::PepXMLReader(
  const string& file_path, ///< the path of the pep.xml
  Database* database, ///< the protein database
  Database* decoy_database ///< the decoy protein database (can be null)
  ) {

  file_path_ = file_path;
  database_ = database;
  decoy_database_ = decoy_database;

}

/**
 * default destructor
 */
PepXMLReader::~PepXMLReader() {
}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* PepXMLReader::parse() {

  FILE* file_ptr = fopen(file_path_.c_str(), "r");

  XML_Parser xml_parser = XML_ParserCreate(NULL);

  if (! xml_parser) {
    carp(CARP_FATAL, "Couldn't allocate memory for parser");
  }

  current_match_collection_ = new MatchCollection();
  current_match_collection_->preparePostProcess();
  int total_matches=current_match_collection_->getMatchTotal();
  cerr<<"total matches: "<< total_matches<<endl;
  XML_SetUserData(xml_parser, this);
  XML_SetElementHandler(xml_parser, open_handler, close_handler);

  int done = 0;

  while (done == 0) {

    int len = fread(Buff, 1, BUFFSIZE, file_ptr);
    if (ferror(stdin)) {
      carp(CARP_FATAL, "Read error");
      exit(-1);
    }
    done = feof(file_ptr);

    if (! XML_Parse(xml_parser, Buff, len, done)) {
      carp(CARP_FATAL, "Parse error at line %d:\n%s\n",
           (int)XML_GetCurrentLineNumber(xml_parser),
           XML_ErrorString(XML_GetErrorCode(xml_parser)));
      exit(-1);
    }
  }

  fclose(file_ptr);

  return current_match_collection_;
}

/**
 * Handles the spectrum_query open tag event 
 */
void PepXMLReader::spectrumQueryOpen(
  const char** attr ///< attribute array for element
  ) {

  if (spectrum_query_open_) {
    carp(CARP_FATAL, "spectrum_query not closed before another was opened!");
    exit(-1);
  }
  spectrum_query_open_ = true;

  int first_scan = -1;
  int last_scan = -1;
  double precursor_mass = -1;
  int  charge = -1;

  for (int i = 0; attr[i]; i += 2) {
    if (strcmp(attr[i], "start_scan") == 0) {
      first_scan = atoi(attr[i+1]);
    } else if (strcmp(attr[i], "end_scan") == 0) {
      last_scan = atoi(attr[i+1]);
    } else if (strcmp(attr[i], "precursor_neutral_mass") == 0) {
      precursor_mass = atof(attr[i+1]);
    } else if (strcmp(attr[i], "assumed_charge") == 0) {
      charge = atoi(attr[i+1]);
    }
  }

  double precursor_mz = (precursor_mass + (MASS_PROTON * (double) charge)) / (double)charge;
  vector<int> charge_vec;
  charge_vec.push_back(charge);
  current_spectrum_ = new Spectrum(first_scan, last_scan, precursor_mz, charge_vec, "");
  current_zstate_.setNeutralMass(precursor_mass, charge);
}

/**
 * Handles the spectrum_query close tag event
 */
void PepXMLReader::spectrumQueryClose() {
  spectrum_query_open_ = false;
}

/**
 * Handles the search_result open tag event
 */
void PepXMLReader::searchResultOpen() {
  if (search_result_open_) {
    carp(CARP_FATAL, "search_result not closed before another opened!");
  }
  search_result_open_ = true;
}

/**
 * Handles the search_result close tag event
 */
void PepXMLReader::searchResultClose() {
  search_result_open_ = false;
}

/**
 * /returns the start position of the peptide sequence within the protein
 */
int PepXMLReader::findStart(
  Protein* protein,  ///< the protein to find the sequence 
  string peptide_sequence, ///< the peptide sequence to find
  string prev_aa, ///< the amino acid before the sequence in the protein
  string next_aa ///< the next amino acid after the sequence in the protein
  ) {

  if (protein->getSequencePointer() == NULL) {
    return -1;
  }
  int case_ = -1;

  int ans = -1;

  if (prev_aa == "-") {
    case_ = 0;
    ans = 1;
  } else if (next_aa == "-") {
    case_ = 1;
    ans = protein->getLength() - peptide_sequence.length() + 1;
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
        carp(CARP_FATAL, "could not find %s in protein %s\n%s", seq.c_str(), protein->getIdPointer(), protein_seq.c_str());
      }
      case_ = 2;
      ans = (pos+1);
    } else {
      case_ = 3;
      ans = (pos+2);
    }
  }

  return ans;

}

/**
 * Handles the search_hit open tag event
 */
void PepXMLReader::searchHitOpen(
  const char** attr ///< atttribute array for element
  ) {
  if (search_hit_open_) {
    carp(CARP_FATAL, "Search Hit not closed before another open!");
  }
  search_hit_open_ = true;


  int hit_rank = -1;
  string protein_string;

  string prev_aa;
  string next_aa;
  unsigned  by_ions_matched=0; 
  unsigned by_ions_total=0; 
  unsigned current_num_matches=0; 

  double peptide_mass = 0.0;
  for (int i = 0; attr[i]; i += 2) {
    if (strcmp(attr[i], "hit_rank") == 0) {
      hit_rank = atoi(attr[i+1]);
    } else if (strcmp(attr[i], "peptide") == 0) {
      current_peptide_sequence_ = attr[i+1];
    } else if (strcmp(attr[i], "protein") == 0) {
      protein_string = attr[i+1];
    } else if (strcmp(attr[i], "num_tot_proteins") == 0) {
      ; // do nothing.
    }else if (strcmp(attr[i], "num_matched_ions") == 0) {
      by_ions_matched = atoi(attr[i+1]);  
    }else if (strcmp(attr[i], "tot_num_ions") == 0){
      by_ions_total = atoi(attr[i+1]); 
    }else if (strcmp(attr[i], "calc_neutral_pep_mass") == 0) {
      peptide_mass = atof(attr[i+1]);
    } else if (strcmp(attr[i], "peptide_prev_aa") == 0) {
      prev_aa = attr[i+1];
    } else if (strcmp(attr[i], "peptide_next_aa") == 0) {
      next_aa = attr[i+1];
    } else if( strcmp(attr[i], "num_matched_peptides") == 0){
      current_num_matches=atoi(attr[i+1]);
    }
  }

  
  unsigned char length = current_peptide_sequence_.length();
    

  bool is_decoy;

  Protein* protein = 
    MatchCollectionParser::getProtein(database_, decoy_database_, protein_string, is_decoy);
  int start_idx = protein->findStart(current_peptide_sequence_, prev_aa, next_aa);
  Peptide* peptide = new Peptide(length, peptide_mass, protein, start_idx);

  current_match_ = new Match(peptide, current_spectrum_, current_zstate_, is_decoy);
  if (is_decoy) {
    current_match_->setNullPeptide(true);
  }

  if ((hit_rank > 0) && (current_match_->getRank(XCORR) == 0)) {
    current_match_->setRank(XCORR, hit_rank);
  }
  if(by_ions_total>0){
    current_match_->setBYIonMatched(by_ions_matched);
    current_match_->setBYIonPossible(by_ions_total);
    current_match_-> setBYIonFractionMatched((FLOAT_T)by_ions_matched/(FLOAT_T)by_ions_total); 
  }else 
    current_match_->setBYIonFractionMatched(0);
  if(current_num_matches>0)
    current_match_->setLnExperimentSize(logf(current_num_matches));
  else 
    current_match_->setLnExperimentSize(0);
}

/**
 * Handles the search_hit close tag event
 */
void PepXMLReader::searchHitClose() {
  search_hit_open_ = false;
  
  //We should have all the information needed to add the match object.
  current_match_collection_->addMatch(current_match_);

}

/**
 * Handles the modification_info open tag event
 */
void PepXMLReader::modificationInfoOpen(
  const char** attr ///< attribute array for element
  ) {
  
  modification_info_open_ = true;
  for (int idx = 0; attr[idx]; idx += 2) {
    if (strcmp(attr[idx], "mod_nterm_mass") == 0) {

      FLOAT_T mod_nterm_mass = atof(attr[idx+1]);
      //set nterm mod mass
      const AA_MOD_T* aa_mod = get_aa_mod_from_mass(mod_nterm_mass);
      MODIFIED_AA_T* mod_seq = current_match_->getPeptide()->getModifiedAASequence();
      modify_aa(&mod_seq[0], aa_mod);
      current_match_->getPeptide()->setModifiedAASequence(mod_seq, true);
      free(mod_seq);
    } 
  }
}

/**
 * Handles the modification_info close tag event
 */
void PepXMLReader::modificationInfoClose() {
  modification_info_open_ = false;
}

/**
 * Handles the mod_aminoacid_mass open tag event
 */
void PepXMLReader::modAminoAcidMassOpen(
  const char** attr ///< attribute array for element
  ) {
  
  mod_aminoacid_mass_open_ = true;
  int position = -1;
  FLOAT_T mod_mass = 0;
  bool have_mod_mass = false;
  
  for (int idx = 0; attr[idx]; idx += 2) {
    if (strcmp(attr[idx], "position") == 0) {
      position = atoi(attr[idx+1]);
    } else if (strcmp(attr[idx],"mass") == 0) {
      mod_mass = atof(attr[idx+1]);
      have_mod_mass = true;
    }
  }
  
  if (position > 0 && have_mod_mass) {
    //set aminoacid modification.
    const AA_MOD_T* aa_mod = get_aa_mod_from_mass(mod_mass);
    MODIFIED_AA_T* mod_seq = current_match_->getPeptide()->getModifiedAASequence();
    modify_aa(&mod_seq[position-1], aa_mod);
    current_match_->getPeptide()->setModifiedAASequence(mod_seq, true);
    free(mod_seq);
  } else {
    carp(CARP_WARNING, "mod_aminoacid_mass error");
  }
}

/**
 * Handles the mod_aminoacid_mass close tag event
 */
void PepXMLReader::modAminoAcidMassClose() {
  mod_aminoacid_mass_open_ = false;
}

/**
 * Handles the alternative_protein open tag event
 */
void PepXMLReader::alternativeProteinOpen(
  const char** attr ///< atttribute array for element
  ) {


  alternative_protein_open_ = true;

  string protein_string;
  string prev_aa;
  string next_aa;

  for (int i = 0; attr[i]; i += 2) {
    if (strcmp(attr[i], "protein") == 0) {
      protein_string = attr[i+1];
    } else if (strcmp(attr[i], "peptide_prev_aa") == 0) {
      prev_aa = attr[i+1];
    } else if (strcmp(attr[i], "peptide_next_aa") == 0) {
      next_aa = attr[i+1];
    }
  }
  
  bool is_decoy;

  Protein* protein = 
    MatchCollectionParser::getProtein(database_, decoy_database_, protein_string, is_decoy);
  if (is_decoy) {
    current_match_->setNullPeptide(true);
  }
  int start_idx = protein->findStart(current_peptide_sequence_, prev_aa, next_aa);

  PeptideSrc* src = new PeptideSrc((DIGEST_T)0, protein, start_idx);
  current_match_->getPeptide()->addPeptideSrc(src);

}

/**
 * Handles the alternative_protein close tag event
 */
void PepXMLReader::alternativeProteinClose() {

  alternative_protein_open_ = false;

}

/**
 * Handles the search_score open tag event
 */
void PepXMLReader::searchScoreOpen(
  const char** attr ///< attribute array for element
  ) {

  search_score_open_ = true;
  string name;
  double value = 0.0;

  for (int i = 0; attr[i]; i += 2) {
    if (strcmp(attr[i], "name") == 0) {
      name = attr[i+1];
    } else if (strcmp(attr[i], "value") == 0) {
      value = atof(attr[i+1]);
    }
  }

  if (name == "xcorr_score" || name == "xcorr") {
    current_match_collection_->setScoredType(XCORR, true);
    current_match_->setScore(XCORR, value);
  } else if (name == "xcorr_rank") {
    current_match_->setRank(XCORR, value);
  } else if (name == "expect") {
    current_match_->setScore(EVALUE, value);
  } else if (name == "delta_cn" || name == "deltacn") {
    current_match_->setDeltaCn(value);
  } else if (name == "sp" || name == "spscore") {
    current_match_collection_->setScoredType(SP, true);
    current_match_->setScore(SP, value);
  } else if (name == "sp_rank" || name == "sprank") {
    current_match_->setRank(SP, value);
  } else if (name == "percolator_score") {
    current_match_collection_->setScoredType(PERCOLATOR_SCORE, true);
    current_match_->setScore(PERCOLATOR_SCORE, value);
  } else if (name == "percolator_qvalue") {
    current_match_collection_->setScoredType(PERCOLATOR_QVALUE, true);
    current_match_->setScore(PERCOLATOR_QVALUE, value);
  } else if (name == "percolator_PEP") {
    current_match_collection_->setScoredType(PERCOLATOR_PEP, true);
    current_match_->setScore(PERCOLATOR_PEP, value);
  } else if (name == "qranker_score") {
    current_match_collection_->setScoredType(QRANKER_SCORE, true);
    current_match_->setScore(QRANKER_SCORE, value);
  } else if (name == "qranker_qvalue") {
    current_match_collection_->setScoredType(QRANKER_QVALUE, true);
    current_match_->setScore(QRANKER_QVALUE, value);
  } else if (name == "qranker_PEP") {
    current_match_collection_->setScoredType(QRANKER_PEP, true);
    current_match_->setScore(QRANKER_PEP, value);
  } 
  //set the custom score
  current_match_->setCustomScore(name, value);

}

/**
 * Handles the search_score close tag event
 */
void PepXMLReader::searchScoreClose() {
  search_score_open_ = false;
}


/**
 * Handles the peptideprophet_result open tag event
 */
void PepXMLReader::peptideProphetResultOpen(
  const char** attr ///< attribute array for element
  ) {


  if (peptideprophet_result_open_) {
    carp(CARP_FATAL, "peptideprophet_result_open_ before close!");
  }  

  peptideprophet_result_open_ = true;


  FLOAT_T probability = 0;
  bool probability_parsed = false;
   for (int idx = 0; attr[idx]; idx += 2) {
    if (strcmp(attr[idx], "probability") == 0) {
      probability = atof(attr[idx+1]);
      probability_parsed = true;
    }
  }
  
  //Place the peptideprophet probability score as a custom score.
  if (probability_parsed) {
    current_match_->setCustomScore("peptideprophet", probability);
  } else {
    carp_once(CARP_WARNING, "Couldn't parse probability from peptideprophet_result");
  }
}


/**
 * Handles the peptideprophet_result close tag event
 */
void PepXMLReader::peptideProphetResultClose() {

  peptideprophet_result_open_ = false;
}


/**
 * sets the target protein database
 */
void PepXMLReader::setDatabase(
  Database* database ///< the target protein database
  ) {

  database_ = database;
}

/**
 * sets the decoy protein database
 */
void PepXMLReader::setDecoyDatabase(
  Database* decoy_database ///< sets the decoy protein database
  ) {
  decoy_database_ = decoy_database;
}


/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* PepXMLReader::parse(
  const char* file_path, ///< path of the xml file
  Database* database, ///< target protein database
  Database* decoy_database ///< decoy protein database (can be null)
  ) {

  PepXMLReader* reader = new PepXMLReader(file_path);
  reader->setDatabase(database);
  reader->setDecoyDatabase(decoy_database);

  MatchCollection* collection = reader->parse();

  return collection;


}

//TODO - remove this code after some time of debugging.
#ifdef MAIN
int main(int argc, char** argv) {

  initialize_parameters();
  set_verbosity_level(CARP_INFO);
  char* file_path = argv[1];
  char* database_path = argv[2];

  Database* database;
  Database* decoy_database;
 
  //MatchCollectionFactory::loadDatabase(database_path, database, decoy_database);
  database = new Database();
  decoy_database = new Database();


  PepXMLReader* reader = new PepXMLReader(file_path);
  reader->setDatabase(database);
  reader->setDecoyDatabase(decoy_database);

  MatchCollection* match_collection = reader->parse();

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
