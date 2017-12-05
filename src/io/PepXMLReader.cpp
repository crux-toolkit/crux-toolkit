/*************************************************************************
 * \file PepXMLReader.cpp
 * \brief Object for parsing pepxml files
 *************************************************************************/

#include "expat.h"
#include "PepXMLReader.h"
#include "util/mass.h"
#include "util/AminoAcidUtil.h"
#include "util/MathUtil.h"
#include "model/Protein.h"
#include "model/PostProcessProtein.h"
#include "model/Peptide.h"
#include "DelimitedFile.h"
#include "parameter.h"
#include "MatchCollectionParser.h"

#include <cstdio>
#include <cstring>
#include <iostream>

using namespace std;
using namespace Crux;

void open_handler(void *data, const char *el, const char **attr) {
  PepXMLReader* reader = (PepXMLReader*)data;
  if (strcmp(el, "aminoacid_modification") == 0) {
    reader->aminoacidModificationOpen(attr);
  } else if (strcmp(el, "spectrum_query") == 0) {
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
  if (strcmp(el, "aminoacid_modification") == 0) {
    reader->aminoacidModificationClose();
  } else if (strcmp(el, "spectrum_query") == 0) {
    reader->spectrumQueryClose();
  } else if (strcmp(el, "search_result") == 0) {
    reader->searchResultClose();
  } else if (strcmp(el, "search_hit") == 0) {
    reader->searchHitClose();
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
  maxRank_ = Params::GetInt("top-match-in");
  aminoacid_modification_open_ = false;
  spectrum_query_open_ = false;
  search_result_open_ = false;
  search_hit_open_ = false;
  alternative_protein_open_ = false;
  search_score_open_ = false;
  peptideprophet_result_open_ = false;

  current_spectrum_ = NULL;
}

/**
 * \returns an initialized object
 */
PepXMLReader::PepXMLReader() : PSMReader() {
  init();
}

/**
 * \returns an object initialized with the file_path
 */
PepXMLReader::PepXMLReader(
  const string& file_path ///< the path of the pep.xml file
  ) : PSMReader(file_path) {
  init();
}

/**
 * \returns an object initialized with the xml path, and the target,decoy databases
 */
PepXMLReader::PepXMLReader(
  const string& file_path, ///< the path of the pep.xml
  Database* database, ///< the protein database
  Database* decoy_database ///< the decoy protein database (can be null)
  ) : PSMReader(file_path, database, decoy_database) {
  init();
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
  if (file_ptr == NULL) {
    carp(CARP_FATAL, "Opening %s or reading failed", file_path_.c_str());
  }

  XML_Parser xml_parser = XML_ParserCreate(NULL);

  if (!xml_parser) {
    carp(CARP_FATAL, "Couldn't allocate memory for parser");
  }

  current_match_collection_ = new MatchCollection();
  current_match_collection_->preparePostProcess();
  XML_SetUserData(xml_parser, this);
  XML_SetElementHandler(xml_parser, open_handler, close_handler);

  int done = 0;
  while (!done) {
    char buf[8192];
    int len = fread(buf, 1, sizeof(buf), file_ptr);
    if (ferror(stdin)) {
      carp(CARP_FATAL, "Read error");
    }
    done = feof(file_ptr);
    if (!XML_Parse(xml_parser, buf, len, done)) {
      carp(CARP_FATAL, "Parse error at line %d:\n%s\n",
           (int)XML_GetCurrentLineNumber(xml_parser),
           XML_ErrorString(XML_GetErrorCode(xml_parser)));
    }
  }
  fclose(file_ptr);
  return current_match_collection_;
}

void PepXMLReader::aminoacidModificationOpen(const char** attr) {
  if (aminoacid_modification_open_) {
    carp(CARP_FATAL, "aminoacid_modification not closed before another was opened!");
  }
  aminoacid_modification_open_ = true;

  ModPosition pos = ANY;
  char aa = '\0';
  double mass = 0, dMass = 0;
  bool variable = false;
  char symbol = '\0';

  for (int i = 0; attr[i]; i += 2) {
    if (strcmp(attr[i], "aminoacid") == 0) {
      aa = attr[i + 1][0];
    } else if (strcmp(attr[i], "mass") == 0) {
      mass = atof(attr[i + 1]);
    } else if (strcmp(attr[i], "massdiff") == 0) {
      dMass = atof(attr[i + 1]);
    } else if (strcmp(attr[i], "variable") == 0) {
      variable = strcmp(attr[i + 1], "Y") == 0;
    } else if (strcmp(attr[i], "peptide_terminus") == 0) {
      // TODO: is this protein terminus?
      if (strcmp(attr[i + 1], "n") == 0) {
        pos = PEPTIDE_N;
      } else if (strcmp(attr[i + 1], "c") == 0) {
        pos = PEPTIDE_C;
      } else if (strcmp(attr[i + 1], "nc") == 0) {
        //pos = PEPTIDE_NC; TODO
      }
    } else if (strcmp(attr[i], "symbol") == 0) {
      symbol = attr[i + 1][0];
    }
  }

  ModificationDefinition::New(string(1, aa), dMass, pos, !variable);
}

void PepXMLReader::aminoacidModificationClose() {
  aminoacid_modification_open_ = false;
}

/**
 * Handles the spectrum_query open tag event
 */
void PepXMLReader::spectrumQueryOpen(
  const char** attr ///< attribute array for element
  ) {
  if (spectrum_query_open_) {
    carp(CARP_FATAL, "spectrum_query not closed before another was opened!");
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
  current_spectrum_ = new Crux::Spectrum(first_scan, last_scan, precursor_mz, charge_vec, "");
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
      carp(CARP_DEBUG, "could not find %s in protein %s\n%s", seq.c_str(), protein->getIdPointer().c_str(), protein_seq.c_str());
      //finding the sequence with the flanks failed, try finding without the flanks.
      seq = peptide_sequence;
      pos = protein_seq.find(seq);
      if (pos == string::npos) {
        carp(CARP_FATAL, "could not find %s in protein %s\n%s",
             seq.c_str(), protein->getIdPointer().c_str(), protein_seq.c_str());
      }
      case_ = 2;
      ans = (pos+1);
    } else {
      case_ = 3;
      ans = pos+2;
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
  string protein_string, prev_aa, next_aa;
  unsigned by_ions_matched = 0;
  unsigned by_ions_total = 0;
  unsigned current_num_matches = 0;

  for (int i = 0; attr[i]; i += 2) {
    if (strcmp(attr[i], "hit_rank") == 0) {
      hit_rank = atoi(attr[i+1]);
    } else if (strcmp(attr[i], "peptide") == 0) {
      current_peptide_sequence_ = attr[i+1];
    } else if (strcmp(attr[i], "protein") == 0) {
      protein_string = attr[i+1];
    } else if (strcmp(attr[i], "num_tot_proteins") == 0) {
      ; // do nothing.
    } else if (strcmp(attr[i], "num_matched_ions") == 0) {
      current_match_collection_->setScoredType(BY_IONS_MATCHED, true);
      by_ions_matched = atoi(attr[i+1]);
    } else if (strcmp(attr[i], "tot_num_ions") == 0) {
      current_match_collection_->setScoredType(BY_IONS_TOTAL, true);
      by_ions_total = atoi(attr[i+1]);
    } else if (strcmp(attr[i], "peptide_prev_aa") == 0) {
      prev_aa = attr[i+1];
    } else if (strcmp(attr[i], "peptide_next_aa") == 0) {
      next_aa = attr[i+1];
    } else if (strcmp(attr[i], "num_matched_peptides") == 0) {
      current_num_matches = atoi(attr[i+1]);
    }
  }

  unsigned char length = current_peptide_sequence_.length();
  bool is_decoy;

  Protein* protein =
    MatchCollectionParser::getProtein(database_, decoy_database_, protein_string, is_decoy);
  int start_idx = protein->findStart(current_peptide_sequence_, prev_aa, next_aa);
  Peptide* peptide = new Crux::Peptide(length, protein, start_idx);

  current_match_ = new Match(peptide, current_spectrum_, current_zstate_, is_decoy);
  if (is_decoy) {
    current_match_->setNullPeptide(true);
  }

  if (hit_rank > 0 && current_match_->getRank(XCORR) == 0) {
    current_match_->setRank(XCORR, hit_rank);
  }
  if (by_ions_total > 0) {
    current_match_->setScore(BY_IONS_MATCHED, by_ions_matched);
    current_match_->setScore(BY_IONS_TOTAL, by_ions_total);
  }

  if (current_num_matches > 0) {
    current_match_->setLnExperimentSize(logf(current_num_matches));
    current_match_collection_->setHasDistinctMatches(true);
  } else {
    current_match_->setLnExperimentSize(0);
  }
}

/**
 * Handles the search_hit close tag event
 */
void PepXMLReader::searchHitClose() {
  search_hit_open_ = false;
  //We should have all the information needed to add the match object.
  if (maxRank_ == 0 || current_match_->getRank(XCORR) <= maxRank_) {
    current_match_collection_->addMatch(current_match_);
  }
}

void PepXMLReader::modificationInfoOpen(const char** attr) {
  double modMass = 0.0;
  ModPosition modPosition = UNKNOWN;
  for (int i = 0; attr[i]; i += 2) {
    if (strcmp(attr[i], "mod_nterm_mass") == 0) {
      // mass of modification + n-terminus
      modPosition = PEPTIDE_N;
      modMass = atof(attr[i+1]) - MASS_H;
    } else if (strcmp(attr[i], "mod_cterm_mass") == 0) {
      // mass of modification + c-terminus
      modPosition = PEPTIDE_C;
      modMass = atof(attr[i+1]) - MASS_O - MASS_H;
    }

    if (modPosition == UNKNOWN || MathUtil::Round(modMass, Params::GetInt("mod-precision")) == 0.0) {
      continue;
    }

    char* seq = current_match_->getPeptide()->getSequence();
    unsigned char position = 0;
    if (modPosition == PEPTIDE_C) {
      position = strlen(seq) - 1;
    }
    char aa = seq[position];
    free(seq);

    // look for a variable mod with this mass
    const ModificationDefinition* mod = ModificationDefinition::Find(modMass, false);
    if (mod == NULL) {
      // no variable mod with this mass, try finding a static one
      mod = ModificationDefinition::Find(modMass, true);
    }
    if (mod != NULL) {
      if (!mod->Static()) {
        current_match_->getPeptide()->addMod(mod, position);
      }
    } else {
      // mod was not defined at top of file, add it as a variable modification
      modMass = MathUtil::Round(modMass, Params::GetInt("mod-precision"));
      mod = ModificationDefinition::NewVarMod(string(1, aa), modMass, modPosition);
      current_match_->getPeptide()->addMod(mod, position);
    }
  }
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
    } else if (strcmp(attr[idx], "mass") == 0) {
      mod_mass = atof(attr[idx+1]);
      have_mod_mass = true;
    }
  }

  if (position > 0 && have_mod_mass) {
    char* seq = current_match_->getPeptide()->getSequence();
    ModPosition mod_position = ANY;
    if (position <= 1) {
      mod_position = PEPTIDE_N;
    } else if (position >= strlen(seq)) {
      mod_position = PEPTIDE_C;
    }
    char aa = seq[position - 1];
    free(seq);
    double staticDeltaMass = ModificationDefinition::DeltaMass(aa, mod_position);
    mod_mass -= staticDeltaMass;
    // mass includes amino acid mass, subtract it
    mod_mass -= AminoAcidUtil::GetMass(aa, Params::GetString("isotopic-mass") == "mono");
    if (MathUtil::Round(mod_mass, Params::GetInt("mod-precision")) != 0.0) {
      // look for a variable mod with this mass
      const ModificationDefinition* mod = ModificationDefinition::Find(mod_mass, false);
      if (mod == NULL) {
        // no variable mod with this mass, try finding a static one
        mod = ModificationDefinition::Find(staticDeltaMass, true);
      }
      if (mod != NULL) {
        if (!mod->Static()) {
          current_match_->getPeptide()->addMod(mod, position - 1);
        }
      } else {
        // mod was not defined at top of file, add it as a variable modification
        mod_mass = MathUtil::Round(mod_mass, Params::GetInt("mod-precision"));
        mod = ModificationDefinition::NewVarMod(string(1, aa), mod_mass, UNKNOWN);
        current_match_->getPeptide()->addMod(mod, position - 1);
      }
    }
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

  PeptideSrc* src = new PeptideSrc(INVALID_DIGEST, protein, start_idx);
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
    current_match_->setScore(DELTA_CN, value);
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
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* PepXMLReader::parse(
  const string& file_path, ///< path of the xml file
  Database* database, ///< target protein database
  Database* decoy_database ///< decoy protein database (can be null)
  ) {
  PepXMLReader reader(file_path);
  reader.setDatabase(database);
  reader.setDecoyDatabase(decoy_database);
  return reader.parse();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
