/*************************************************************************
 * \file SQTReader.cpp
 * \brief Object for parsing sqt files
 *************************************************************************/

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include "parameter.h"
#include "LineFileReader.h"
#include "MatchCollectionParser.h"
#include "SQTReader.h"
#include "util/StringUtils.h"

using namespace std;
using namespace Crux;

/**
 * Initializes the object
 */
void SQTReader::init() {
  maxRank_ = Params::GetInt("top-match-in");
  last_parsed_ = SQT_LINE_NONE;
  current_spectrum_ = NULL;
  current_match_ = NULL;
}

/**
 * \returns an initialized object
 */
SQTReader::SQTReader() : PSMReader() {
  init();
}

/**
 * \returns an object initialized with the file_path
 */
SQTReader::SQTReader(
  const string& file_path ///< the path of the pep.xml file
  ) : PSMReader(file_path) {
  init();
}

/**
 * \returns an object initialized with the xml path, and the target,decoy databases
 */
SQTReader::SQTReader(
  const string& file_path, ///< the path of the sqt
  Database* database, ///< the protein database
  Database* decoy_database ///< the decoy protein database (can be null)
  ) : PSMReader(file_path, database, decoy_database) {
  init();
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
  current_match_collection_ = new MatchCollection();
  current_match_collection_->preparePostProcess();
  current_match_collection_->setScoredType(XCORR, true);
  current_match_collection_->setScoredType(SP, true);
  last_parsed_ = SQT_LINE_NONE;

  LineFileReader line_reader(file_path_);
  while (line_reader.hasNext()) {
    string line = line_reader.next();
    if (!line.empty()) {
      switch (line[0]) {
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
          carp(CARP_ERROR, "Unknown line %d\n%s", line_reader.getCurrentRow(), line.c_str());
      }
    }
  }

  //add the last match
  if (current_match_ != NULL &&
      (maxRank_ == 0 || current_match_->getRank(XCORR) <= maxRank_)) {
    current_match_collection_->addMatchToPostMatchCollection(current_match_);
  }
  return current_match_collection_;
}


void SQTReader::parseHeader(const string& line) {
  last_parsed_ = SQT_LINE_HEADER;
  string content = StringUtils::Trim(line.substr(2));
  if (StringUtils::StartsWith(content, "StaticMod")) {
    content = StringUtils::Trim(content.substr(10));
    vector<string> tokens = StringUtils::Split(content, '=');
    if (tokens.size() == 2) {
      ModificationDefinition::NewStaticMod(
        tokens[0], StringUtils::FromString<double>(tokens[1]), UNKNOWN /*TODO*/);
    }
  } else if (StringUtils::StartsWith(content, "DiffMod")) {
    content = StringUtils::Trim(content.substr(8));
    vector<string> tokens = StringUtils::Split(content, '=');
    if (tokens.size() == 2) {
      string t0 = StringUtils::Trim(tokens[0]);
      string aminoAcids = t0.substr(0, t0.length() - 1);
      char symbol = t0[t0.length() - 1];
      double deltaMass = StringUtils::FromString<double>(tokens[1]);
      ModificationDefinition::NewVarMod(
        aminoAcids, deltaMass, UNKNOWN /*TODO*/, false, false, symbol);
    }
  }
}

void SQTReader::parseSpectrum(const string& line) {
  const int spectrum_low_scan_idx = 1;
  const int spectrum_high_scan_idx = 2;
  const int spectrum_charge_idx = 3;
  const int spectrum_observed_mass_idx = 6;
  const int spectrum_total_ion_intensity_idx = 7;
  const int spectrum_lowest_sp_idx = 8;
  const int spectrum_num_matches_idx = 9;

  vector<string> tokens = StringUtils::Split(line, '\t');

  int low_scan = -1;
  int high_scan = -1;
  int charge = -1;
  current_ln_experiment_size_ = 0.0;
  for (int i = 0; i < tokens.size(); i++) {
    switch (i) {
      case spectrum_low_scan_idx:
        low_scan = StringUtils::FromString<int>(tokens[i]);
        break;
      case spectrum_high_scan_idx:
        high_scan = StringUtils::FromString<int>(tokens[i]);
        break;
      case spectrum_charge_idx:
        charge = StringUtils::FromString<int>(tokens[i]);
        break;
      case spectrum_observed_mass_idx: {
        double observed_mass = StringUtils::FromString<double>(tokens[i]);
        current_zstate_.setSinglyChargedMass(observed_mass, charge);
        vector<int> chargeVec(1, charge);
        current_spectrum_ = new Spectrum(low_scan, high_scan,
                                         current_zstate_.getMZ(), chargeVec, "");
        break;
      }
      case spectrum_num_matches_idx: {
        int current_num_matches = StringUtils::FromString<int>(tokens[i]);
        if (current_num_matches > 0) {
          current_match_collection_->setHasDistinctMatches(true);
          current_ln_experiment_size_ = logf((FLOAT_T)current_num_matches);
        }
        break;
      }
      case spectrum_total_ion_intensity_idx:
        if (!tokens[i].empty()) {
          current_spectrum_->setHasTotalEnergy(true);
          current_spectrum_->setTotalEnergy(StringUtils::FromString<double>(tokens[i]));
        }
        break;
      case spectrum_lowest_sp_idx:
        if (!tokens[i].empty()) {
          current_spectrum_->setHasLowestSp(true);
          current_spectrum_->setLowestSp(StringUtils::FromString<double>(tokens[i]));
        }
        break;
    }
  }
  last_parsed_ = SQT_LINE_SPECTRUM;
}

void SQTReader::parseMatch(const string& line) {
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

  vector<string> tokens = StringUtils::Split(line, '\t');

  int xcorr_rank = StringUtils::FromString<int>(tokens[match_xcorr_rank_idx]);

  if (maxRank_ == 0 || xcorr_rank <= maxRank_) {
    int sp_rank = StringUtils::FromString<int>(tokens[match_sp_rank_idx]);
    double calculated_mass = StringUtils::FromString<double>(tokens[match_calculated_mass_idx]);
    double delta_cn = StringUtils::FromString<double>(tokens[match_delta_cn_idx]);
    double xcorr = StringUtils::FromString<double>(tokens[match_xcorr_idx]);
    double sp = StringUtils::FromString<double>(tokens[match_sp_idx]);
    int matched_ions = StringUtils::FromString<int>(tokens[match_matched_ions_idx]);
    int expected_ions = StringUtils::FromString<int>(tokens[match_expected_ions_idx]);

    string sqt_sequence = tokens[match_sequence_idx];
    vector<string> sequence_tokens = StringUtils::Split(sqt_sequence, '.');

    current_prev_aa_ = sequence_tokens.front();
    current_next_aa_ = sequence_tokens.back();
    sequence_tokens.erase(sequence_tokens.begin());
    sequence_tokens.pop_back();
    current_peptide_sequence_ = StringUtils::Join(sequence_tokens, '.');
    if (current_match_ != NULL) {
      current_match_collection_->addMatchToPostMatchCollection(current_match_);
    }

    string unmodifiedSequence;
    vector<Modification> mods;
    for (size_t i = 0; i < current_peptide_sequence_.size(); i++) {
      char c = current_peptide_sequence_[i];
      if ('A' <= c && c <= 'Z') {
        unmodifiedSequence.push_back(c);
      } else {
        const ModificationDefinition* modDef = ModificationDefinition::Find(c);
        if (modDef != NULL) {
          mods.push_back(Modification(modDef, (i > 0) ? i - 1 : 0));
        } else {
          carp(CARP_ERROR, "Unknown modification in %s: %c",
               current_peptide_sequence_.c_str(), c);
        }

      }
    }
    Peptide* peptide = new Peptide(unmodifiedSequence, mods);
    current_peptide_sequence_ = unmodifiedSequence;

    current_match_ = new Match(peptide, current_spectrum_, current_zstate_, false);
    current_match_->setScore(XCORR, xcorr);
    current_match_->setScore(SP, sp);
    current_match_->setScore(DELTA_CN, delta_cn);
    current_match_->setRank(XCORR, xcorr_rank);
    current_match_->setRank(SP, sp_rank);

    current_match_->setScore(BY_IONS_MATCHED, matched_ions);
    current_match_->setScore(BY_IONS_TOTAL, expected_ions);
    current_match_->setLnExperimentSize(current_ln_experiment_size_);
  }

  last_parsed_ = SQT_LINE_MATCH;
}

void SQTReader::parseLocus(const string& line) {
  const int locus_protein_id_idx = 1;
  const int locus_protein_desc_idx = 2;

  vector<string> tokens = StringUtils::Split(line, '\t');

  string protein_id = tokens[locus_protein_id_idx];
  string protein_desc = "";
  if (tokens.size() > 2) {
    protein_desc = tokens[locus_protein_desc_idx];
  }
  bool is_decoy;
  Protein* protein =
    MatchCollectionParser::getProtein(database_, decoy_database_, protein_id, is_decoy);

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
      carp(CARP_DEBUG, "could not find %s in protein %s\n%s", seq.c_str(), protein->getIdPointer().c_str(), protein_seq.c_str());
      //finding the sequence with the flanks failed, try finding without the flanks.
      seq = peptide_sequence;
      pos = protein_seq.find(seq);
      if (pos == string::npos) {
        carp(CARP_FATAL, "could not %s in protein %s\n%s", seq.c_str(), protein->getIdPointer().c_str(), protein_seq.c_str());
      }
      return pos+1;
    }
    return pos+2;
  }
}

/**
 * \returns the MatchCollection resulting from the parsed xml file
 */
MatchCollection* SQTReader::parse(
  const string& file_path, ///< path of the xml file
  Database* database, ///< target protein database
  Database* decoy_database ///< decoy protein database (can be null)
  ) {
  SQTReader* reader = new SQTReader(file_path);
  reader->setDatabase(database);
  reader->setDecoyDatabase(decoy_database);
  return reader->parse();
}

/**
 * \use the modification symbols in this SQT file
 * looks for DiffMod <residues><symbol>=<delta mass.
 */
void SQTReader::readSymbols(const string& file, bool append) {
  ifstream sqt(file.c_str());
  if (!sqt) {
    carp(CARP_FATAL, "Could not open '%s' for reading.", file.c_str());
  }

  bool doneReset = false;

  AA_MOD_T** list_of_mods;
  int num_mods = get_all_aa_mod_list(&list_of_mods);

  string line;
  carp(CARP_DEBUG, "Reading modifications from %s", file.c_str());
  while (std::getline(sqt, line)) {
    size_t idx;
    if ((idx = line.find_first_not_of(" \t")) != string::npos &&
        line[idx] == 'H' &&
        (idx = line.find("DiffMod", idx + 1)) != string::npos &&
        (idx = line.find_first_of(" \t", idx + 1)) != string::npos) {
      if (num_mods == MAX_AA_MODS) {
        carp(CARP_FATAL, "Too many modifications in file '%s'", file.c_str());
      }
      if ((idx = line.find_first_not_of(" \t", idx + 1)) == string::npos) {
        continue;
      }
      char residue = toupper(line[idx]);
      if (++idx >= line.length()) {
        continue;
      }
      char symbol = line[idx];
      if ((idx = line.find('=', idx + 1)) == string::npos ||
          idx >= line.length() - 1) {
        continue;
      }
      FLOAT_T mass = StringUtils::FromString<FLOAT_T>(line.substr(idx + 1));

      if (residue < 'A' || residue > 'Z') {
        carp(CARP_FATAL, "The DiffMod residue '%c' in file '%s' is invalid",
                         residue, file.c_str());
      }

      if (!doneReset && !append) {
        resetMods();
        doneReset = true;
      }

      AA_MOD_T* mod = NULL;
      for (int i = 0; i < num_mods; ++i) {
        if (list_of_mods[i]->getMassChange() == mass) {
          mod = list_of_mods[i];
          mod->getAAList()[residue - 'A'] = true;
          break;
        }
      }
      if (!mod) {
        mod = new AA_MOD_T(num_mods);
        mod->setMassChange(mass);
        mod->setSymbol(symbol);
        mod->getAAList()[residue - 'A'] = true;
        list_of_mods[num_mods] = mod;
        incrementNumMods();
        ++num_mods;
      }

      carp(CARP_DETAILED_DEBUG, "Read mod %c (symbol %c) with mass %f",
           residue, symbol, mass);
    }
  }
  carp(CARP_DEBUG, "Read %d modifications from %s", num_mods, file.c_str());
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
