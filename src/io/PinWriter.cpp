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
#include <boost/foreach.hpp>

using namespace std;
using namespace Crux; 

PinWriter::PinWriter():
  out_(NULL),
  enzyme_(get_enzyme_type_parameter("enzyme")),
  precision_(Params::GetInt("precision")),
  mass_precision_(Params::GetInt("mass-precision"))
{
  features_.push_back(make_pair("SpecId", true));
  features_.push_back(make_pair("Label", true));
  features_.push_back(make_pair("ScanNr", true));
  features_.push_back(make_pair("ExpMass", true));
  features_.push_back(make_pair("CalcMass", true));
  features_.push_back(make_pair("lnrSp", false));
  features_.push_back(make_pair("deltLCn", true));
  features_.push_back(make_pair("deltCn", true));
  features_.push_back(make_pair("XCorr", true));
  features_.push_back(make_pair("Sp", false));
  features_.push_back(make_pair("IonFrac", false));
  features_.push_back(make_pair("RefactoredXCorr", false));
  features_.push_back(make_pair("NegLog10PValue", false));
  features_.push_back(make_pair("PepLen", true));
  for (int i = '1'; i <= '9'; i++) {
    features_.push_back(make_pair("Charge" + string(1, i), false));
  }
  features_.push_back(make_pair("enzN", true));
  features_.push_back(make_pair("enzC", true));
  features_.push_back(make_pair("enzInt", true));
  features_.push_back(make_pair("lnNumSP", true));
  features_.push_back(make_pair("dM", true));
  features_.push_back(make_pair("absdM", true));
  features_.push_back(make_pair("Peptide", true));
  features_.push_back(make_pair("Proteins", true));
}

PinWriter::~PinWriter(){ 
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
){ 
  Peptide* peptide = match->getPeptide();
  Spectrum* spectrum = match->getSpectrum();
  int charge = match->getCharge();
  bool enzC = false;
  bool enzN = false;
  FLOAT_T obsMass = match->getZState().getSinglyChargedMass();
  FLOAT_T calcMass = peptide->getPeptideMass() + MASS_PROTON;
  FLOAT_T dM = (obsMass - calcMass) / charge;

  char* sequence = peptide->getSequence();
  int missedCleavages = get_num_internal_cleavage(sequence, enzyme_);
  get_terminal_cleavages(sequence, peptide->getNTermFlankingAA(),
                         peptide->getCTermFlankingAA(), enzyme_, enzN, enzC);
  free(sequence);

  vector<string> fields;
  BOOST_FOREACH(const std::string& feature, enabledFeatures_) {
    if (feature == "SpecId") {
      fields.push_back(getId(match, spectrum->getFirstScan()));
    } else if (feature == "Label") {
      fields.push_back(match->getNullPeptide() ? "-1" : "1");
    } else if (feature == "ScanNr") {
      fields.push_back(StringUtils::ToString(spectrum->getFirstScan()));
    } else if (feature == "ExpMass") {
      fields.push_back(StringUtils::ToString(obsMass, mass_precision_));
    } else if (feature == "CalcMass") {
      fields.push_back(StringUtils::ToString(calcMass, mass_precision_));
    } else if (feature == "lnrSp") {
      double sp = match->getRank(SP);
      fields.push_back(StringUtils::ToString(sp > 0 ? log(sp) : 0, precision_));
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
      FLOAT_T ion_frac = match->getBYIonFractionMatched();
      fields.push_back(StringUtils::ToString(!isnan(ion_frac) ? ion_frac : 0, precision_));
    } else if (feature == "RefactoredXCorr") {
      fields.push_back(
        StringUtils::ToString(match->getScore(TIDE_SEARCH_REFACTORED_XCORR), precision_));
    } else if (feature == "NegLog10PValue") {
      fields.push_back(StringUtils::ToString(
        -log10(match->getScore(TIDE_SEARCH_EXACT_PVAL)), precision_));
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
    } else if (feature == "lnNumSP") {
      fields.push_back(StringUtils::ToString(match->getLnExperimentSize(), precision_));
    } else if (feature == "dM") {
      fields.push_back(StringUtils::ToString(dM, precision_));
    } else if (feature == "absdM") {
      fields.push_back(StringUtils::ToString(fabs(dM), precision_));
    } else if (feature == "Peptide") {
      fields.push_back(getPeptide(peptide));
    } else if (feature == "Proteins") {
      fields.push_back(StringUtils::Join(peptide->getProteinIds(), '\t'));
    } else {
      carp(CARP_FATAL, "Unknown feature: '%s'", feature.c_str());
    }
  }
  *out_ << StringUtils::Join(fields, '\t') << endl;
}

string PinWriter::getPeptide(Peptide* peptide) {
  stringstream sequence;
  sequence << peptide->getNTermFlankingAA() << '.';

  char* modified_sequence = Params::GetBool("mod-symbols")
    ? peptide->getModifiedSequenceWithSymbols()
    : peptide->getModifiedSequenceWithMasses(
        get_mass_format_type_parameter("mod-mass-format"));
  sequence << modified_sequence;
  free(modified_sequence);

  sequence << '.' << peptide->getCTermFlankingAA();

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

