/**
 * \file PinWriter.cpp
 * \brief Writes search results in the .pin format.
 */
#include <iostream>
#include "PinWriter.h"
#include "model/MatchCollection.h"
#include "util/crux-utils.h"
#include "parameter.h"
#include "MatchCollectionParser.h"
#include <cctype>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <vector>
#include "model/Spectrum.h"
#include "SQTReader.h"
#include "model/SpectrumZState.h"
#include <fstream>
#include <limits>
#include <boost/filesystem.hpp>
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include "util/MathUtil.h"
#include "MassHandler.h"
#include <boost/foreach.hpp>

using namespace std;
using namespace Crux; 

#define MAX_LOG_P 100.0 // Needed for handling p-value = 0, which should not happen, but alas does.

PinWriter::PinWriter():
  out_(NULL),
  enzyme_(get_enzyme_type_parameter("enzyme")),
  precision_(Params::GetInt("precision")),
  mass_precision_(Params::GetInt("mass-precision")) {
    
  features_.push_back(make_pair("SpecId", true));
  features_.push_back(make_pair("Label", true));
  features_.push_back(make_pair("filename", true));
  features_.push_back(make_pair("ScanNr", true));
  features_.push_back(make_pair("ExpMass", true));
  features_.push_back(make_pair("CalcMass", true));
  features_.push_back(make_pair("lnrSp", false));
  features_.push_back(make_pair("deltLCn", true));
  features_.push_back(make_pair("deltCn", true));
  features_.push_back(make_pair("XCorr", true));
  features_.push_back(make_pair("TailorScore", true));  
  features_.push_back(make_pair("byIonsMatched", true));  
  features_.push_back(make_pair("byIonsTotal", true));  
  features_.push_back(make_pair("byIonsFraction", true));  
  features_.push_back(make_pair("byIonsRepeatMatch", true));  

  features_.push_back(make_pair("Sp", false));
  features_.push_back(make_pair("IonFrac", false));
  features_.push_back(make_pair("RefactoredXCorr", false));
  features_.push_back(make_pair("NegLog10PValue", false));
  features_.push_back(make_pair("NegLog10ResEvPValue", false));
  features_.push_back(make_pair("NegLog10CombinePValue", false));
  features_.push_back(make_pair("PepLen", true));
  for (int i = '1'; i <= '9'; i++) {
    features_.push_back(make_pair("Charge" + string(1, i), false));
  }
  features_.push_back(make_pair("enzN", true));
  features_.push_back(make_pair("enzC", true));
  features_.push_back(make_pair("enzInt", true));
  features_.push_back(make_pair("lnNumSP", true));
  features_.push_back(make_pair("lnNumDSP", false));
  features_.push_back(make_pair("dM", true));
  features_.push_back(make_pair("absdM", true));

  // DIAmeter related, added by Yang
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-precursor"), 0)) {
     features_.push_back(make_pair("PrecursorIntRankM0", true));
     features_.push_back(make_pair("PrecursorIntRankM1", true));
     features_.push_back(make_pair("PrecursorIntRankM2", true));
  }
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-fragment"), 0)) {
	 features_.push_back(make_pair("DynFragPVal", true));
	 features_.push_back(make_pair("StaFragPVal", true));
  }
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-rtdiff"), 0)) {
	 features_.push_back(make_pair("RTDiff", true));
  }
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-elution"), 0)) {
	 features_.push_back(make_pair("CoeluteMS1", true));
	 features_.push_back(make_pair("CoeluteMS2", true));
	 features_.push_back(make_pair("CoeluteMS1MS2", true));
  }
  features_.push_back(make_pair("EnsembleScore", true));


  features_.push_back(make_pair("Peptide", true));
  features_.push_back(make_pair("Proteins", true));

  scannr_cnt_ = 0;
}

PinWriter::~PinWriter() { 
  closeFile(); 
}

/**
 * Open a file of the given name.  Replace an existing file if
 * overwrite is true, else exit if an existing file is found.
 */
void PinWriter::openFile(const string& filename, const string& output_dir, bool overwrite) {
  if (!(out_ = create_stream_in_path(filename.c_str(), output_dir.c_str(), overwrite))) {
    carp(CARP_FATAL, "Can't open file '%s'", filename.c_str());
  }
}

void PinWriter::openFile(CruxApplication* application, string filename, MATCH_FILE_TYPE type) {
  openFile(filename, "", Params::GetBool("overwrite"));
}

bool PinWriter::getEnabledStatus(const string& name) const {
  for (vector< pair<string, bool> >::const_iterator i = features_.begin(); i != features_.end(); i++) {
    if (i->first == name) {
      return i->second;
    }
  }
  return false;
}

void PinWriter::setEnabledStatus(const string& name, bool enabled) {
  IsFeature finder(name);
  vector< pair<string, bool> >::iterator i =
    find_if(features_.begin(), features_.end(), finder);
  if (i != features_.end()) {
    i->second = enabled;
  } else {
    carp(CARP_WARNING, "setEnabledStatus: feature '%s' not found", name.c_str());
  }
}

/**
 * Close the file, if open.
 */
void PinWriter::closeFile() {
  if (out_) {
    delete out_;
    out_ = NULL;
  }
}

void PinWriter::write( 
  MatchCollection* target_collection,
  const vector<MatchCollection*>& decoys,
  int top_rank
) {
  // added by Yang

  vector<MatchCollection*> collections(1, target_collection);
  collections.insert(collections.end(), decoys.begin(), decoys.end());
  for (vector<MatchCollection*>::iterator i = collections.begin();
       i != collections.end();
       i++) {
    MatchIterator match_iterator(*i);  
    while (match_iterator.hasNext()) {
      Match* match = match_iterator.next();
      if (match->getRank(XCORR) <= top_rank) {
        printPSM(match);
      }
    }
  }
}

void PinWriter::write(MatchCollection* collection, string database) {
  bool sp = collection->getScoredType(SP);
  bool xcorr = collection->getScoredType(XCORR);
  bool exact_p = collection->getScoredType(TIDE_SEARCH_REFACTORED_XCORR);
  bool combine_p = collection->getScoredType(BOTH_PVALUE);
  bool tailor = collection->getScoredType(TAILOR_SCORE);

  if (combine_p) {
    exact_p = true;
  }

  setEnabledStatus("lnrSp", sp);
  setEnabledStatus("deltLCn", collection->getScoredType(DELTA_LCN));
  setEnabledStatus("deltCn", collection->getScoredType(DELTA_CN));
  setEnabledStatus("XCorr", xcorr);
  setEnabledStatus("Sp", sp);
  setEnabledStatus("IonFrac", sp);
  setEnabledStatus("RefactoredXCorr", exact_p);
  setEnabledStatus("NegLog10PValue", exact_p);
  setEnabledStatus("NegLog10ResEvPValue", combine_p);
  setEnabledStatus("NegLog10CombinePValue", combine_p);
  setEnabledStatus("TailorScore", tailor);

  // DIAmeter related, added by Yang
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-precursor"), 0)) {
     setEnabledStatus("PrecursorIntRankM0", collection->getScoredType(PRECURSOR_INTENSITY_RANK_M0));
     setEnabledStatus("PrecursorIntRankM1", collection->getScoredType(PRECURSOR_INTENSITY_RANK_M1));
     setEnabledStatus("PrecursorIntRankM2", collection->getScoredType(PRECURSOR_INTENSITY_RANK_M2));
  }
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-fragment"), 0)) {
     setEnabledStatus("DynFragPVal", collection->getScoredType(DYN_FRAGMENT_PVALUE));
     setEnabledStatus("StaFragPVal", collection->getScoredType(STA_FRAGMENT_PVALUE));
  }
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-rtdiff"), 0)) {
     setEnabledStatus("RTDiff", collection->getScoredType(RT_DIFF));
  }
  if (!MathUtil::AlmostEqual(Params::GetDouble("coeff-elution"), 0)) {
     setEnabledStatus("CoeluteMS1", collection->getScoredType(COELUTE_MS1));
     setEnabledStatus("CoeluteMS2", collection->getScoredType(COELUTE_MS2));
     setEnabledStatus("CoeluteMS1MS2", collection->getScoredType(COELUTE_MS1_MS2));
  }
  setEnabledStatus("EnsembleScore", collection->getScoredType(ENSEMBLE_SCORE));


  int max_charge = 0;
  for (MatchIterator i = MatchIterator(collection); i.hasNext();) {
    max_charge = max(i.next()->getCharge(), max_charge);
  }
  for (int i = 1; i <= max_charge; i++) {
    setEnabledStatus("Charge" + StringUtils::ToString(i), true);
  }

  vector<MatchCollection*> decoyvec;
  int top_match = Params::GetInt("top-match");
  printHeader();
  write(collection, decoyvec, top_match); // TODO: When top match is greater than default (5) in a given PSM File?
}

bool PinWriter::isInfinite(FLOAT_T x) {
  return x == numeric_limits<FLOAT_T>::infinity() || -x == numeric_limits<FLOAT_T>::infinity();
}

void PinWriter::printHeader() {
  enabledFeatures_.clear();
  FeatureCopy copier(&enabledFeatures_);
  for_each(features_.begin(), features_.end(), copier);

  for (vector< pair<string, bool> >::const_iterator i = features_.begin();
       i != features_.end();
       i++) {
    carp(CARP_DEBUG, "PIN feature: '%s'%s",
         i->first.c_str(), !i->second ? " (disabled)" : "");
  }
  *out_ << StringUtils::Join(enabledFeatures_, '\t') << endl;
}

void PinWriter::printPSM(
  Match* match
) { 
  Peptide* peptide = match->getPeptide();
  Spectrum* spectrum = match->getSpectrum();
  int charge = match->getCharge();
  bool enzC = false;
  bool enzN = false;
  double obsMass = match->getZState().getSinglyChargedMass();
  FLOAT_T calcMass = peptide->calcModifiedMass() + MASS_PROTON;
  FLOAT_T dM = MassHandler::massDiff(obsMass, calcMass, charge);

  const char* sequence = peptide->getSequence();
  int missedCleavages = get_num_internal_cleavage(sequence, enzyme_);
  get_terminal_cleavages(sequence, peptide->getNTermFlankingAA(),
                         peptide->getCTermFlankingAA(), enzyme_, enzN, enzC);

  vector<string> fields;
  BOOST_FOREACH(const std::string& feature, enabledFeatures_) {
    if (feature == "SpecId") {
      fields.push_back(getId(match, spectrum->getFirstScan()));
    } else if (feature == "Label") {
      fields.push_back(match->getNullPeptide() ? "-1" : "1");
    } else if (feature == "filename") {
      fields.push_back(match->getFilePath());
    } else if (feature == "ScanNr") {
      // modified by Yang
      if (Params::GetBool("unique-scannr")) { fields.push_back(StringUtils::ToString(++scannr_cnt_)); }
      else { fields.push_back(StringUtils::ToString(spectrum->getFirstScan()));  }

    } else if (feature == "ExpMass") {
      fields.push_back(StringUtils::ToString(obsMass, mass_precision_));
    } else if (feature == "CalcMass") {
      fields.push_back(StringUtils::ToString(calcMass, mass_precision_));
    } else if (feature == "lnrSp") {
      double sp = match->getRank(SP);
      fields.push_back(StringUtils::ToString(log(sp + 1.0), precision_));
    } else if (feature == "deltLCn") {
      FLOAT_T delta_lcn = match->getScore(DELTA_LCN);
      if (isInfinite(delta_lcn) || isnan(delta_lcn)) {
        delta_lcn = 0;
      }
      fields.push_back(StringUtils::ToString(delta_lcn, precision_));
    } else if (feature == "deltCn") {
      FLOAT_T delta_cn = match->getScore(DELTA_CN);
      if (isInfinite(delta_cn) || isnan(delta_cn)) {
        delta_cn = 0;
      }
      fields.push_back(StringUtils::ToString(delta_cn, precision_));
    } else if (feature == "XCorr") {
      fields.push_back(StringUtils::ToString(match->getScore(XCORR), precision_));
    } else if (feature == "Sp") {
      fields.push_back(StringUtils::ToString(match->getScore(SP), precision_));
    } else if (feature == "IonFrac") {
      FLOAT_T ionMatch = match->getScore(BY_IONS_MATCHED);
      FLOAT_T ionTotal = match->getScore(BY_IONS_TOTAL);
      FLOAT_T ionFrac = (ionMatch != NOT_SCORED && ionTotal != NOT_SCORED && ionTotal > 0)
        ? ionMatch / ionTotal
        : 0;
      fields.push_back(StringUtils::ToString(ionFrac, precision_));
    } else if (feature == "RefactoredXCorr") {
      fields.push_back(
        StringUtils::ToString(match->getScore(TIDE_SEARCH_REFACTORED_XCORR), precision_));
    } else if (feature == "NegLog10PValue") {
      FLOAT_T logP = -log10(match->getScore(TIDE_SEARCH_EXACT_PVAL));
      fields.push_back(StringUtils::ToString(isInfinite(logP) ? MAX_LOG_P : logP, precision_));
    } else if (feature == "NegLog10ResEvPValue") {
      FLOAT_T logP = -log10(match->getScore(RESIDUE_EVIDENCE_PVAL));
      fields.push_back(StringUtils::ToString(isInfinite(logP) ? MAX_LOG_P : logP, precision_));
    } else if (feature == "NegLog10CombinePValue" ) {
      FLOAT_T logP = -log10(match->getScore(BOTH_PVALUE));
      fields.push_back(StringUtils::ToString(isInfinite(logP) ? MAX_LOG_P : logP, precision_));
    } else if (feature == "TailorScore" ) {
      FLOAT_T tailor = match->getScore(TAILOR_SCORE);
      fields.push_back(StringUtils::ToString(tailor));
    } else if (feature == "byIonsMatched" ) {
      FLOAT_T byIonsMatched = match->getScore(BY_IONS_MATCHED);
      fields.push_back(StringUtils::ToString(byIonsMatched));
    } else if (feature == "byIonsTotal" ) {
      FLOAT_T byIonsTotal = match->getScore(BY_IONS_TOTAL);
      fields.push_back(StringUtils::ToString(byIonsTotal));
    } else if (feature == "byIonsFraction" ) {
      FLOAT_T byIonsFraction = match->getScore(BY_ION_FRACTION);
      fields.push_back(StringUtils::ToString(byIonsFraction));
    } else if (feature == "byIonsRepeatMatch" ) {
      FLOAT_T byIonsRepeatMatch = match->getScore(BY_ION_REPEAT_MATCH);
      fields.push_back(StringUtils::ToString(byIonsRepeatMatch));
    } else if (feature == "PepLen") {
      fields.push_back(StringUtils::ToString((unsigned) peptide->getLength()));
    } else if (StringUtils::StartsWith(feature, "Charge")) {
      int chargeFeature = StringUtils::FromString<int>(feature.substr(6));
      fields.push_back(charge == chargeFeature ? "1" : "0"); 
    } else if (feature == "enzN") {
      fields.push_back(enzN ? "1" : "0");
    } else if (feature == "enzC") {
      fields.push_back(enzC ? "1" : "0");
    } else if (feature == "enzInt") {
      fields.push_back(StringUtils::ToString(missedCleavages));
    } else if (feature == "lnNumSP" || feature == "lnNumDSP") {
      fields.push_back(StringUtils::ToString(match->getLnExperimentSize(), precision_));
    } else if (feature == "dM") {
      fields.push_back(StringUtils::ToString(dM, precision_));
    } else if (feature == "absdM") {
      fields.push_back(StringUtils::ToString(fabs(dM), precision_));
    } else if (feature == "Peptide") {
      fields.push_back(getPeptide(peptide));
    } else if (feature == "Proteins") {
      fields.push_back(StringUtils::Join(peptide->getProteinIds(), '\t'));
    }
    // DIAmeter related, added by Yang
    else if (feature == "PrecursorIntRankM0" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-precursor"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(PRECURSOR_INTENSITY_RANK_M0))); }
    else if (feature == "PrecursorIntRankM1" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-precursor"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(PRECURSOR_INTENSITY_RANK_M1))); }
    else if (feature == "PrecursorIntRankM2" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-precursor"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(PRECURSOR_INTENSITY_RANK_M2))); }
    else if (feature == "RTDiff" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-rtdiff"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(RT_DIFF))); }
    else if (feature == "DynFragPVal" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-fragment"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(DYN_FRAGMENT_PVALUE))); }
    else if (feature == "StaFragPVal" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-fragment"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(STA_FRAGMENT_PVALUE))); }
    else if (feature == "CoeluteMS1" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-elution"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(COELUTE_MS1))); }
    else if (feature == "CoeluteMS2" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-elution"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(COELUTE_MS2))); }
    else if (feature == "CoeluteMS1MS2" && !MathUtil::AlmostEqual(Params::GetDouble("coeff-elution"), 0)) { fields.push_back(StringUtils::ToString(match->getScore(COELUTE_MS1_MS2))); }
    else if (feature == "EnsembleScore") { fields.push_back(StringUtils::ToString(match->getScore(ENSEMBLE_SCORE))); }

    else {
      carp(CARP_FATAL, "Unknown feature: '%s'", feature.c_str());
    }
  }
  *out_ << StringUtils::Join(fields, '\t') << endl;
}

string PinWriter::getPeptide(Peptide* pep) {
  stringstream sequence;
  sequence << pep->getNTermFlankingAA() << '.'
           << (Params::GetBool("mod-symbols")
                ? pep->getModifiedSequenceWithSymbols()
                : pep->getModifiedSequenceWithMasses())
           << '.' << pep->getCTermFlankingAA();
  return sequence.str();
}

string PinWriter::getId(Match* match, int scan_number) {
  string prefix = Params::GetBool("filestem-prefixes")
    ? FileUtils::Stem(match->getFilePath())
    : "";

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

