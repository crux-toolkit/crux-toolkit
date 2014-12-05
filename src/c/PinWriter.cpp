/**
 * \file PinWriter.cpp
 * \brief Writes search results in the .pin format.
 */
#include <iostream>
#include "PinWriter.h"
#include "MatchCollection.h"
#include "crux-utils.h"
#include "parameter.h"
#include "MatchCollectionParser.h"
#include <cctype>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <vector>
#include "Spectrum.h"
#include "SQTReader.h"
#include "SpectrumZState.h"
#include <fstream>
#include <limits>
#include <boost/filesystem.hpp>

using namespace std;
using namespace Crux; 

PinWriter::PinWriter():
  output_file_(NULL),
  enzyme_(get_enzyme_type_parameter("enzyme")),
  isotopic_mass_(get_mass_type_parameter("isotopic-mass")),
  precision_(get_int_parameter("precision")),
  mass_precision_(get_int_parameter("mass-precision")),
  scan_number_(-1),
  is_sp_(get_boolean_parameter("compute-sp") || get_boolean_parameter("sqt-output")),
  decoy_prefix_(get_string_parameter("decoy-prefix"))
{
}

PinWriter::~PinWriter(){ 
  closeFile(); 
}

/**
 * Open a file of the given name.  Replace an existing file if
 * overwrite is true, else exit if an existing file is found.
 */
void PinWriter::openFile(const char* filename, const char* output_dir, bool overwrite) {
  output_file_ = create_file_in_path(filename, output_dir, overwrite);
  if (output_file_ == NULL) {
    carp(CARP_FATAL, "Can't open file '%s'", filename);
  }
}

/**
 * Close the file, if open.
 */
void PinWriter::closeFile() {
  if (output_file_) {
    fclose(output_file_);
    output_file_ = NULL;
  }
}

void PinWriter::write(
  MatchIterator* iterator,
  Spectrum* spectrum,
  int top_rank) {
  vector<Match*> matches; 
  while (iterator->hasNext()) {
    matches.push_back(iterator->next());
  }
  calculateDeltaCN(matches);
  for (vector<Match*>::const_iterator i = matches.begin(); i != matches.end(); ++i) {
    int xcorrRank = (*i)->getRank(XCORR);
    if (xcorrRank <= top_rank)
      printPSM(*i, spectrum);
    else 
      return;
  }
}

void PinWriter::calculateDeltaCN(map<pair<int, int>, vector<Match*> >& scan_charge_to_matches) {
  for (map< pair<int, int>, vector<Match*> >::iterator i = scan_charge_to_matches.begin();
       i != scan_charge_to_matches.end();
       ++i) {
    calculateDeltaCN(i->second);
  }
}

void PinWriter::calculateDeltaCN(vector<Match*>& collection) {
  // We have to use an explicit 'new' here because
  // the MatchCollection includes a couple of huge 'C' style
  // arrays, which are guaranteed to blow the stack.
  MatchCollection* tmp_matches = new MatchCollection();

  tmp_matches->setScoredType(XCORR, true);
  for (vector<Match*>::const_iterator i = collection.begin();
    i != collection.end();
    ++i) {
    tmp_matches->addMatch(*i);
  }
  tmp_matches->calculateDeltaCn();
  
  delete tmp_matches;
}

void PinWriter::calculateDeltaCN(
  MatchCollection* target_collection, 
  vector<MatchCollection*>& decoys
) {
  calculateDeltaCN(target_collection);
  for (vector<MatchCollection*>::iterator i = decoys.begin(); 
       i != decoys.end(); 
       ++i) {
    calculateDeltaCN(*i);
  }
}

void PinWriter::calculateDeltaCN(MatchCollection* matchCollection) {
  map< pair<int, int>, vector<Match*> > scanChargeToMatches;
  MatchIterator matchIterator(matchCollection);  
  while (matchIterator.hasNext()) {
    Match* match = matchIterator.next();
    pair<int, int> scanCharge = make_pair(
      match->getSpectrum()->getFirstScan(), match->getZState().getCharge());
    if (scanChargeToMatches.find(scanCharge) == scanChargeToMatches.end()) {
      scanChargeToMatches[scanCharge] = vector<Match*>();
    }
    scanChargeToMatches[scanCharge].push_back(match);
  }
  calculateDeltaCN(scanChargeToMatches);
}

void PinWriter::write( 
  MatchCollection* target_collection,
  vector<MatchCollection*>& decoys,
  Spectrum* spectrum,
  int top_rank
){
  calculateDeltaCN(target_collection, decoys);

  if (scan_number_ != spectrum->getFirstScan()) { 
    scan_number_ = spectrum->getFirstScan();
  }

  MatchIterator target_match_iterator(target_collection);
  write(&target_match_iterator, spectrum, top_rank);

  // Allow for multiple sets of decoys
  for (vector<MatchCollection*>::iterator i = decoys.begin();
       i != decoys.end(); 
       i++) {
    MatchIterator decoy_match_iterator(*i);  
    write(&decoy_match_iterator, spectrum, top_rank);
  }
}

/*creates a pin file from two sqt, txt ,or pep.xml files*/
void PinWriter::write( 
  MatchCollection* target_collection,
  vector<MatchCollection*>& decoys,
  int top_rank
){
  calculateDeltaCN(target_collection, decoys);
  is_sp_ = target_collection->getScoredType(SP);

  map< int, vector<Match*> > scan_to_matches;
  MatchIterator target_match_iterator(target_collection);  
  while (target_match_iterator.hasNext()) {
    Match* match = target_match_iterator.next();
    int scan = match->getSpectrum()->getFirstScan();
    if (scan_to_matches.find(scan) == scan_to_matches.end()) {
      scan_to_matches[scan] = vector<Match*>();
    }
    scan_to_matches[scan].push_back(match);
    charges_.insert(match->getZState().getCharge());
  }
    
  // Allow for multiple sets of decoys
  for (vector<MatchCollection*>::iterator i = decoys.begin();
       i != decoys.end();
       i++) {
    MatchIterator decoy_match_iterator(*i);  
    while (decoy_match_iterator.hasNext()) {
      Match* match = decoy_match_iterator.next();
      int scan = match->getSpectrum()->getFirstScan();
      if (scan_to_matches.find(scan) == scan_to_matches.end()) {
        vector<Match*> matches;
        scan_to_matches[scan] = matches;
      }
      scan_to_matches[scan].push_back(match); 
      charges_.insert(match->getZState().getCharge());
    }
  }

  printHeader();
  for (map<int, vector<Match*> >::iterator iter = scan_to_matches.begin();
    iter != scan_to_matches.end();
    ++iter) {
    vector<Match*> matches = iter->second;
    if (scan_number_ != iter->first) {
      scan_number_ = iter->first;
    }
    for (size_t i = 0; i < matches.size(); ++i) {   
      int xcorrRank = matches[i]->getRank(XCORR);
      if (xcorrRank <= top_rank) {
        printPSM(matches[i], matches[i]->getSpectrum());
      }
    }
  }
}

bool PinWriter::isInfinite(FLOAT_T x) {
  return x == numeric_limits<FLOAT_T>::infinity();
}

void PinWriter::printHeader() {
  stringstream features;
  if (is_sp_) {
    features << "lnrSp\t";
  }
  features << "deltLCn\tdeltCn\tXcorr\t";
  if (is_sp_) {
    features << "Sp\tIonFrac\t";
  }
  features << "Mass\tPepLen\t";
  for (set<int>::const_iterator i = charges_.begin(); i != charges_.end(); ++i) {
    features << "Charge" << *i << '\t';
  }
  features << "enzN\tenzC\tenzInt\tlnNumSP\tdm\tabsdM\t";

  /*if (po->calcPTMs) 
    push_backFeatureDescription("ptm");
  if (po->pngasef) 
    push_backFeatureDescription("PNGaseF");
  if (po->calcAAFrequencies)
    for (std::string::const_iterator it = aaAlphabet.begin(); it != aaAlphabet.end(); it++)
    {
      std::string temp = boost::lexical_cast<std::string>(*it)+"-Freq";
      push_backFeatureDescription(temp.c_str());
    }*/
  fprintf(output_file_,
    "SpecId\tLabel\tScanNr\t"
    "%s"
    "Peptide\tProteins\n",
    features.str().c_str()
  );
}

void PinWriter::printPSM(
  Match* match,
  Spectrum* spectrum
){ 
  Peptide* peptide = match->getPeptide();

  int charge = match->getCharge();
  bool enzC = false;
  bool enzN = false;
  FLOAT_T obsMass = match->getZState().getSinglyChargedMass();
  FLOAT_T calcMass = peptide->calcMass(isotopic_mass_) + calcMassOfMods(peptide);
  FLOAT_T dM = (obsMass - calcMass) / charge;

  char* sequence = peptide->getSequence();
  int missedCleavages = get_num_internal_cleavage(sequence, enzyme_);
  get_terminal_cleavages(sequence, peptide->getNTermFlankingAA(),
                         peptide->getCTermFlankingAA(), enzyme_, enzN, enzC);
  free(sequence);
 
  fprintf(output_file_, "%s\t%d\t%d\t",
    getId(match, spectrum->getFirstScan()).c_str(), // SpecId
    match->getNullPeptide() ? -1 : 1, // Label
    spectrum->getFirstScan() // ScanNr
  );
  if (is_sp_) {
    fprintf(output_file_, "%.*f\t",
      precision_, match->getRank(SP) > 0 ? log((double) match->getRank(SP)) : 0 // lnrSp
    );
  }
  fprintf(output_file_, "%.*f\t%.*f\t%.*f\t",
    precision_, isInfinite(fabs(match->getDeltaLCn())) ? 0 : match->getDeltaLCn(), // deltLCn
    precision_, isInfinite(fabs(match->getDeltaCn())) ? 0 : match->getDeltaCn(), // deltCn
    precision_, match->getScore(XCORR) // XCorr
  );
  if (is_sp_) {
    fprintf(output_file_, "%.*f\t%.*f\t",
      precision_, match->getScore(SP), // Sp
      precision_, isnan(match->getBYIonFractionMatched()) ? 0 : match->getBYIonFractionMatched() // IonFrac
    );
  }
  fprintf(output_file_, "%.*f\t%u\t",
    mass_precision_, obsMass, // Mass
    peptide->getLength() // PepLen
  );
  for (set<int>::const_iterator i = charges_.begin(); i != charges_.end(); ++i) {
    fprintf(output_file_, "%u\t",
      charge == *i ? 1 : 0 // ChargeN
    );
  }
  fprintf(output_file_, "%u\t%u\t%u\t%.*f\t%.*f\t%.*f\t%s\t%s",
    enzN ? 1 : 0, // enzN
    enzC ? 1 : 0, // enzC
    missedCleavages, // enzInt
    precision_, match->getLnExperimentSize(), // lnNumSP
    precision_, dM, // dM
    precision_, fabs(dM), // absdM
    getPeptide(peptide).c_str(), // Peptide
    getProteins(peptide).c_str() // Proteins
  );

  fprintf(output_file_, "\n");
}

string PinWriter::getPeptide(Peptide* peptide) {
  stringstream sequence;

  sequence << peptide->getNTermFlankingAA() << '.';

  char* unmodified_sequence = peptide->getModifiedSequenceWithMasses(MOD_MASS_ONLY);
  sequence << unmodified_sequence;
  free(unmodified_sequence);

  sequence << '.' << peptide->getCTermFlankingAA();

  return sequence.str();
}

string PinWriter::getProteins(
  Peptide* peptide
) {
  vector<string> proteinIds;
  vector<string> proteinDescriptions; 
  int numProteins = peptide->getProteinInfo(proteinIds, proteinDescriptions);
  if (numProteins < 1) {
    return "";
  }
  vector<string>::const_iterator i = proteinIds.begin();
  string proteinIdStr = *(i++);
  for (; i != proteinIds.end(); ++i) {
    proteinIdStr += ',' + *i;
  }
  return proteinIdStr;
}

string PinWriter::getId(
  Match* match,
  int scan_number
){
  string prefix;
  if (get_boolean_parameter("filestem-prefixes")) {
    string filename = match->getFilePath();
    if (!filename.empty()) {
      boost::filesystem::path path(filename);
      prefix = path.stem().string();
    }
  }

  stringstream psm_id; 
  if (prefix.empty()) {
    psm_id << (match->getNullPeptide() ? "decoy" : "target")
           << '_' << match->getFileIndex();
  } else {
    psm_id << prefix;
  }
  psm_id << '_' << scan_number << '_' << match->getCharge() << '_'
         << match->getRank(XCORR);
  return psm_id.str();   
}

FLOAT_T PinWriter::calcMassOfMods(Peptide* peptide) {
  FLOAT_T total_mass_delta = 0;
  MODIFIED_AA_T* mod_seq = peptide->getModifiedAASequence();
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  for (size_t i = 0; mod_seq[i] != MOD_SEQ_NULL; ++i) {
    for (size_t j = 0; j < total_mods; ++j) {
      if (is_aa_modified(mod_seq[i], mod_list[j])) {
        total_mass_delta += aa_mod_get_mass_change(mod_list[j]);
      }
    }
  }
  free(mod_seq);
  return total_mass_delta;
}

