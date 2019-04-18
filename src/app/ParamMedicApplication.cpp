#include "ParamMedicApplication.h"
#include "ParamMedicInference.h"
#include "io/carp.h"
#include "io/SpectrumCollectionFactory.h"
#include "parameter.h"
#include "util/mass.h"
#include "util/Params.h"

#include <cmath>
#include <fstream>
#include <numeric>

using namespace Crux;
using namespace std;

const double SQRT_2_PI = 2.50662827463;

// maximum proportion of precursor delta-masses that can be 0, otherwise we give up
const double MAX_PROPORTION_PRECURSOR_DELTAS_ZERO = 0.5;

// maximum peaks to use to fit a mixed distribution
const int MAX_PEAKPAIRS = 100000;

// multipliers to transform standard error values into algorithm parameters
const double PRECURSOR_SIGMA_MULTIPLIER = 37.404067;
const double FRAGMENT_SIGMA_MULTIPLIER = 0.004274;

// separation between averagine peaks used for binning spectra
const double AVERAGINE_PEAK_SEPARATION = 1.0005079;

// minimum allowed values for sigma of the estimated normal
const double MIN_SIGMA_PPM = 0.01;
const double MIN_SIGMA_TH = 0.00001;

// if more than this proportion of mass bins have more than one peak,
// we might be looking at profile-mode data
const double PROPORTION_MASSBINS_MULTIPEAK_PROFILE = 0.5;

ParamMedicApplication::ParamMedicApplication() {
}

ParamMedicApplication::~ParamMedicApplication() {
}

int ParamMedicApplication::main(int argc, char** argv) {
  bool overwrite = Params::GetBool("overwrite");
  ofstream* out = create_stream_in_path(make_file_path("param-medic.txt").c_str(), NULL, overwrite);

  ParamMedic::RunAttributeResult errorCalcResult;
  vector<ParamMedic::RunAttributeResult> modsResult;
  processFiles(Params::GetStrings("spectrum-file"), true, true, &errorCalcResult, &modsResult);

  vector<string> columns;
  columns.push_back("file");
  columns.push_back("precursor_prediction_ppm");
  columns.push_back("precursor_sigma_ppm");
  columns.push_back("fragment_prediction_th");
  columns.push_back("fragment_sigma_ppm");
  columns.push_back("SILAC_4Da_present");
  columns.push_back("SILAC_4Da_statistic");
  columns.push_back("SILAC_6Da_present");
  columns.push_back("SILAC_6Da_statistic");
  columns.push_back("SILAC_8Da_present");
  columns.push_back("SILAC_8Da_statistic");
  columns.push_back("SILAC_10Da_present");
  columns.push_back("SILAC_10Da_statistic");
  columns.push_back("iTRAQ_8plex_present");
  columns.push_back("iTRAQ_8plex_statistic");
  columns.push_back("iTRAQ_4plex_present");
  columns.push_back("iTRAQ_4plex_statistic");
  columns.push_back("TMT_6plex_present");
  columns.push_back("TMT_6plex_statistic");
  columns.push_back("TMT_10plex_present");
  columns.push_back("TMT_10plex_statistic");
  columns.push_back("TMT_2plex_present");
  columns.push_back("TMT_2plex_statistic");
  columns.push_back("phospho_present");
  columns.push_back("phospho_statistic");

  *out << StringUtils::Join(columns, '\t') << endl;

  columns.clear();

  const string ERROR_VAL = "ERROR";
  columns.push_back(Params::GetStrings("spectrum-file").front());
  string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_FAILURE);
  columns.push_back(fail.empty()
    ? StringUtils::ToString(errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_PREDICTION))
    : ERROR_VAL);
  columns.push_back(fail.empty()
    ? StringUtils::ToString(errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_SIGMA))
    : ERROR_VAL);
  fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_FAILURE);
  columns.push_back(fail.empty()
    ? StringUtils::ToString(errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_PREDICTION))
    : ERROR_VAL);
  columns.push_back(fail.empty()
    ? StringUtils::ToString(errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_SIGMA))
    : ERROR_VAL);
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_4Da_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_4Da_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_6Da_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_6Da_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_8Da_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_8Da_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_10Da_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "SILAC_10Da_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "iTRAQ_8plex_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "iTRAQ_8plex_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "iTRAQ_4plex_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "iTRAQ_4plex_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "TMT_6plex_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "TMT_6plex_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "TMT_10plex_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "TMT_10plex_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "TMT_2plex_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "TMT_2plex_statistic"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "phospho_present"));
  columns.push_back(ParamMedic::RunAttributeResult::getFirstValue(modsResult, "phospho_statistic"));

  *out << StringUtils::Join(columns, '\t') << endl;

  delete out;
  return 0;
}

void ParamMedicApplication::processFiles(
  const vector<string>& files,
  bool errorCalcEnabled,
  bool modsEnabled,
  ParamMedic::RunAttributeResult* errorCalcResult,
  vector<ParamMedic::RunAttributeResult>* modsResult
) {
  using ParamMedic::Modification;
  using ParamMedic::RunAttributeDetector;
  using ParamMedic::RunAttributeResult;

  ParamMedic::ErrorCalc errorCalc;
  ParamMedic::SilacDetector silacDetector;
  ParamMedic::ReporterIonProportionCalc reporterDetector;
  ParamMedic::PhosphoLossProportionCalc phosphoDetector;

  vector<RunAttributeDetector*> modDetectors;
  if (modsEnabled) {
    modDetectors.push_back(&silacDetector);
    modDetectors.push_back(&reporterDetector);
    modDetectors.push_back(&phosphoDetector);
  }

  vector<RunAttributeDetector*> detectors;
  if (errorCalcEnabled) {
    detectors.push_back(&errorCalc);
  }
  for (vector<RunAttributeDetector*>::const_iterator i = modDetectors.begin();
       i != modDetectors.end();
       i++) {
    detectors.push_back(*i);
  }

  int numSpectraProcessed = ParamMedic::processSpectra(files, detectors);
  if (numSpectraProcessed == 0) {
    carp(CARP_FATAL, "No spectra found! Quitting.");
  }

  carp(CARP_DEBUG, "Processed all spectra. Summarizing...");
  carp(CARP_INFO, "");
  carp(CARP_INFO, "Search parameter recommendations:");
  if (errorCalcEnabled) {
    RunAttributeResult result = errorCalc.summarize();
    if (errorCalcResult) {
      *errorCalcResult = result;
    }
    string errorCalcMessagesAll = result.getValue("search_param_messages");
    vector<string> errorCalcMessages = StringUtils::Split(errorCalcMessagesAll, '\n');
    for (vector<string>::const_iterator i = errorCalcMessages.begin();
         i != errorCalcMessages.end();
         i++) {
      carp(CARP_INFO, "%s", i->c_str());
    }
  }

  vector<ParamMedic::RunAttributeResult> modsResultTmp;
  if (modsEnabled) {
    for (vector<RunAttributeDetector*>::const_iterator i = modDetectors.begin();
         i != modDetectors.end();
         i++) {
      modsResultTmp.push_back((*i)->summarize());
    }
    if (modsResult) {
      *modsResult = modsResultTmp;
    }

    vector<Modification> mods = Modification::GetFromResults(modsResultTmp);
    if (!mods.empty()) {
      for (vector<Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
        carp(CARP_INFO, "%s", i->str().c_str());
      }
    } else {
      carp(CARP_INFO, "No modifications detected requiring search parameter changes.");
    }
  }
}

string ParamMedicApplication::getName() const {
  return "param-medic";
}

string ParamMedicApplication::getDescription() const {
  return
    "[[html:<p>]]Examine the spectra in a file to estimate the best precursor "
    "and fragment error tolerances for database search.[[html:</p>]]";
}

bool ParamMedicApplication::hidden() const {
  return false;
}

vector<string> ParamMedicApplication::getArgs() const {
  string arr[] = {
    "spectrum-file+"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> ParamMedicApplication::getOptions() const {
  string arr[] = {
    "verbosity",
    "spectrum-parser",
    "pm-min-precursor-mz",
    "pm-max-precursor-mz",
    "pm-min-frag-mz",
    "pm-max-frag-mz",
    "pm-min-scan-frag-peaks",
    "pm-max-precursor-delta-ppm",
    "pm-charges",
    "pm-top-n-frag-peaks",
    "pm-pair-top-n-frag-peaks",
    "pm-min-common-frag-peaks",
    "pm-max-scan-separation",
    "pm-min-peak-pairs"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > ParamMedicApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("stdout",
    "the estimated parameter values for precursor mass tolerance (in ppm) and "
    "fragment bin size (in Th), as well as the standard deviations of the "
    "estimated error distributions for precursor and fragment masses."));
  return outputs;
}

bool ParamMedicApplication::needsOutputDirectory() const {
  return true;
}

namespace ParamMedic {

Modification::Modification():
  location_(""), massDiff_(0.0), variable_(false) {
}

Modification::Modification(const string& location, double massDiff, bool isVariable):
  location_(location), massDiff_(massDiff), variable_(isVariable) {
}

Modification::Modification(const Modification& other):
  location_(other.location_), massDiff_(other.massDiff_), variable_(other.variable_) {
}

Modification::~Modification() {
}

void swap(Modification& x, Modification& y) {
  using std::swap;
  swap(x.location_, y.location_);
  swap(x.massDiff_, y.massDiff_);
  swap(x.variable_, y.variable_);
}

Modification& Modification::operator=(Modification rhs) {
  swap(*this, rhs);
  return *this;
}

string Modification::getLocation() const {
  return location_;
}

double Modification::getMassDiff() const {
  return massDiff_;
}

bool Modification::getVariable() const {
  return variable_;
}

string Modification::str() const {
  stringstream ss;
  ss << (variable_ ? "Variable" : "Static") << " modification of "
     << massDiff_ << "Da on ";
  if (location_ == LOCATION_NTERM) {
    ss << "N terminus";
  } else if (location_ == LOCATION_CTERM) {
    ss << "C terminus";
  } else {
    ss << location_;
  }
  return ss.str();
}

vector<Modification> Modification::GetFromResults(const vector<RunAttributeResult>& results) {
  vector<ParamMedic::Modification> mods;
  for (vector<ParamMedic::RunAttributeResult>::const_iterator i = results.begin(); i != results.end(); i++) {
    const vector<ParamMedic::Modification>& resultMods = i->getMods();
    mods.insert(mods.end(), resultMods.begin(), resultMods.end());
  }
  return mods;
}

const string Modification::LOCATION_NTERM = "N terminus";
const string Modification::LOCATION_CTERM = "C terminus";

RunAttributeResult::RunAttributeResult() {
}

RunAttributeResult::RunAttributeResult(const RunAttributeResult& other):
  modifications_(other.modifications_) {
  nameValuePairs_.insert(other.nameValuePairs_.begin(), other.nameValuePairs_.end());
}

RunAttributeResult::~RunAttributeResult() {
}

void swap(RunAttributeResult& x, RunAttributeResult& y) {
  using std::swap;
  swap(x.modifications_, y.modifications_);
  swap(x.nameValuePairs_, y.nameValuePairs_);
}

RunAttributeResult& RunAttributeResult::operator=(RunAttributeResult rhs) {
  swap(*this, rhs);
  return *this;
}

string RunAttributeResult::getValue(const string& name) const {
  map<string, string>::const_iterator i = nameValuePairs_.find(name);
  if (i != nameValuePairs_.end()) {
    return i->second;
  }
  return "";
}

string RunAttributeResult::getFirstValue(const vector<RunAttributeResult>& results, const string& name) {
  for (vector<RunAttributeResult>::const_iterator i = results.begin(); i != results.end(); i++) {
    string v = i->getValue(name);
    if (!v.empty()) {
      return v;
    }
  }
  return "";
}

void RunAttributeResult::setValue(const string& name, const string& value) {
  nameValuePairs_[name] = value;
}

void RunAttributeResult::addMod(const string& location, double massDiff, bool isVariable) {
  modifications_.push_back(Modification(location, massDiff, isVariable));
}

const vector<Modification>& RunAttributeResult::getMods() const {
  return modifications_;
}

const string RunAttributeResult::ERROR = "ERROR";

ErrorCalc::ErrorCalc():
  numTotalSpectra_(0) {
  vector<string> chargeStrings = StringUtils::Split(Params::GetString("pm-charges"), ',');
  for (vector<string>::const_iterator i = chargeStrings.begin(); i != chargeStrings.end(); i++) {
    int charge = StringUtils::FromString<int>(StringUtils::Trim(*i));
    calcs_[charge] = new PerChargeErrorCalc(charge);
  }
}

ErrorCalc::~ErrorCalc() {
  for (map<int, PerChargeErrorCalc*>::const_iterator i = calcs_.begin(); i != calcs_.end(); i++) {
    delete i->second;
  }
}

void ErrorCalc::processSpectrum(
  const Spectrum* spectrum,
  const vector<double>& binnedSpectrum
) {
  for (map<int, PerChargeErrorCalc*>::iterator i = calcs_.begin(); i != calcs_.end(); i++) {
    i->second->processSpectrum(spectrum, binnedSpectrum);
    numTotalSpectra_++;
  }
}

void ErrorCalc::nextFile() {
  for (map<int, PerChargeErrorCalc*>::iterator i = calcs_.begin(); i != calcs_.end(); i++) {
    i->second->nextFile();
  }
}

void ErrorCalc::calcMassErrorDist(
  string* precursorFailure,
  string* fragmentFailure,
  double* precursorSigmaPpm,
  double* fragmentSigmaPpm,
  double* precursorPredictionPpm,
  double* fragmentPredictionTh
) const {
  precursorFailure->clear();
  fragmentFailure->clear();
  *precursorSigmaPpm = numeric_limits<double>::quiet_NaN();
  *precursorPredictionPpm = numeric_limits<double>::quiet_NaN();
  *fragmentSigmaPpm = numeric_limits<double>::quiet_NaN();
  *fragmentPredictionTh = numeric_limits<double>::quiet_NaN();

  carp(CARP_INFO, "Processed %d total spectra", numTotalSpectra_);
  bool hasNonZeroChargePairs = false;
  for (map<int, PerChargeErrorCalc*>::const_iterator i = calcs_.begin(); i != calcs_.end(); i++) {
    const int charge = i->first;
    const PerChargeErrorCalc* calc = i->second;
    carp(CARP_INFO, "  charge %d", charge);
    carp(CARP_INFO, "Processed %d qualifying spectra", calc->getNumPassingSpectra());

    size_t precursorPairs = calc->getPairedPrecursorMzs().size();
    size_t fragmentPairs = calc->getPairedFragmentPeaks().size();
    carp(CARP_INFO, "Precursor pairs: %d", precursorPairs);
    carp(CARP_INFO, "Fragment pairs: %d", fragmentPairs);
    if (charge > 0 && precursorPairs > 0) {
      hasNonZeroChargePairs = true;
    }

    carp(CARP_INFO, "Total spectra in the same bin as another: %d",
         calc->getNumSpectraSameBin());
    carp(CARP_INFO, "Total spectra in the same bin as another and within m/z tol: %d",
         calc->getNumSpectraWithinPpm());
    carp(CARP_INFO, "Total spectra in the same bin as another and within m/z tol and within scan range: %d",
         calc->getNumSpectraWithinPpmAndScans());
    // check the proportion of mass bins, in the whole file, that have multiple fragments.
    // If that's high, we might be looking at profile-mode data.
    const int single = calc->getNumSingleFragBins();
    const int multiple = calc->getNumMultipleFragBins();
    double proportionMultipleFrags;
    if (single + multiple != 0) {
      proportionMultipleFrags = (double)multiple / (double)(single + multiple);
      carp(CARP_INFO, "Proportion of bins with multiple fragments: %.02f", proportionMultipleFrags);
      if (proportionMultipleFrags > PROPORTION_MASSBINS_MULTIPEAK_PROFILE) {
        carp(CARP_WARNING, "Is this profile-mode data? Proportion of mass bins with multiple "
                           "peaks is quite high (%.02f)", proportionMultipleFrags);
        carp(CARP_WARNING, "Param-Medic will not perform well on profile-mode data.");
      }
    } else {
      proportionMultipleFrags = 0;
      carp(CARP_DEBUG, "No values in any bin!");
    }
  }

  vector< pair<double, double> > pairedPrecursorMzs;
  vector< pair<Peak, Peak> > pairedFragmentPeaks;
  if (hasNonZeroChargePairs) {
    carp(CARP_INFO, "Found paired spectra from known charges, so using those.");
    for (map<int, PerChargeErrorCalc*>::const_iterator i = calcs_.begin(); i != calcs_.end(); i++) {
      if (i->first < 1) {
        continue;
      }
      const vector< pair<double, double> >& calcPairedPrecursors =
        i->second->getPairedPrecursorMzs();
      const vector< pair<Peak, Peak> >& calcPairedFragments =
        i->second->getPairedFragmentPeaks();
      pairedPrecursorMzs.insert(pairedPrecursorMzs.end(),
                                calcPairedPrecursors.begin(), calcPairedPrecursors.end());
      pairedFragmentPeaks.insert(pairedFragmentPeaks.end(),
                                 calcPairedFragments.begin(), calcPairedFragments.end());
    }
  } else {
    carp(CARP_INFO, "Did not find spectra from known charges, so looking for "
                    "unknown-charge spectra.");
    map<int, PerChargeErrorCalc*>::const_iterator lookup = calcs_.find(0);
    pairedPrecursorMzs = lookup->second->getPairedPrecursorMzs();
    pairedFragmentPeaks = lookup->second->getPairedFragmentPeaks();
  }

  if (pairedPrecursorMzs.size() > MAX_PEAKPAIRS) {
    carp(CARP_INFO, "Reducing %d to %d peak pairs for precursor...",
         pairedPrecursorMzs.size(), MAX_PEAKPAIRS);
    random_shuffle(pairedPrecursorMzs.begin(), pairedPrecursorMzs.end(), myrandom_limit);
    pairedPrecursorMzs.resize(MAX_PEAKPAIRS);
  }

  vector<double> precursorDistancesPpm;
  int numZeroPrecursorDeltas = 0;
  for (vector< pair<double, double> >::const_iterator i = pairedPrecursorMzs.begin();
       i != pairedPrecursorMzs.end();
       i++) {
    double diffTh = i->first - i->second;
    if (diffTh == 0.0) {
      ++numZeroPrecursorDeltas;
    }
    precursorDistancesPpm.push_back(diffTh * MILLION / i->first);
  }

  // check for conditions that would cause us to bomb out
  if (precursorDistancesPpm.size() < Params::GetInt("pm-min-peak-pairs")) {
    stringstream msg;
    msg << "Need >= " << Params::GetString("pm-min-peak-pairs")
        << " peak pairs to fit mixed distribution. Got only "
        << precursorDistancesPpm.size() << "." << endl << "Details:" << endl;
    for (map<int, PerChargeErrorCalc*>::const_iterator i = calcs_.begin(); i != calcs_.end(); i++) {
      const PerChargeErrorCalc* calc = i->second;
      msg << "  Charge " << i->first << endl
          << "Spectra in same averagine bin as another: " << calc->getNumSpectraSameBin() << endl
          << "    ... and also within m/z tolerance: " << calc->getNumSpectraWithinPpm() << endl
          << "    ... and also within scan range: " << calc->getNumSpectraWithinPpmAndScans() << endl
          << "    ... and also with sufficient in-common fragments: " << precursorDistancesPpm.size() << endl;
    }
    *precursorFailure = msg.str();
  }
  if (precursorFailure->empty()) {
    double proportionPrecursorMzsZero = (double)numZeroPrecursorDeltas / pairedPrecursorMzs.size();
    carp(CARP_DEBUG, "proportion zero: %f", proportionPrecursorMzsZero);
    if (proportionPrecursorMzsZero > MAX_PROPORTION_PRECURSOR_DELTAS_ZERO) {
      *precursorFailure =
        "Too high a proportion of precursor mass differences (" +
        StringUtils::ToString(proportionPrecursorMzsZero) + ") are exactly 0. "
        "Some processing has been done on this run that param-medic can't handle. "
        "You should investigate what that processing might be.";
    }
  }

  double precursorMuPpm2Measures = numeric_limits<double>::quiet_NaN();
  double precursorSigmaPpm2Measures = numeric_limits<double>::quiet_NaN();
  if (precursorFailure->empty()) {
    carp(CARP_INFO, "Using %d peak pairs for precursor error estimation.",
         precursorDistancesPpm.size());
    estimateMuSigma(precursorDistancesPpm, MIN_SIGMA_PPM,
                    &precursorMuPpm2Measures, &precursorSigmaPpm2Measures);
  }

  if (pairedFragmentPeaks.size() < Params::GetInt("pm-min-peak-pairs")) {
    stringstream msg;
    msg << "Need >= " + Params::GetString("pm-min-peak-pairs")
        << " peak pairs to fit mixed distribution. Got only " 
        << pairedFragmentPeaks.size() << "." << endl << "Details:" << endl;
    for (map<int, PerChargeErrorCalc*>::const_iterator i = calcs_.begin(); i != calcs_.end(); i++) {
      const PerChargeErrorCalc* calc = i->second;
      msg << "  Charge " << i->first << endl
          << "Spectra in same averagine bin as another: " << calc->getNumSpectraSameBin() << endl
          << "    ... and also within m/z tolerance: " << calc->getNumSpectraWithinPpm() << endl
          << "    ... and also within scan range: " << calc->getNumSpectraWithinPpmAndScans() << endl
          << "    ... and also with sufficient in-common fragments: " << pairedFragmentPeaks.size() << endl;
    }
    *fragmentFailure = msg.str();
  }

  double fragmentMuPpm2Measures = numeric_limits<double>::quiet_NaN();
  double fragmentSigmaPpm2Measures = numeric_limits<double>::quiet_NaN();
  if (fragmentFailure->empty()) {
    if (pairedFragmentPeaks.size() > MAX_PEAKPAIRS) {
      carp(CARP_DEBUG, "Reducing %d to %d peak pairs for fragment...",
           pairedFragmentPeaks.size(), MAX_PEAKPAIRS);
      random_shuffle(pairedFragmentPeaks.begin(), pairedFragmentPeaks.end(), myrandom_limit);
      pairedFragmentPeaks.resize(MAX_PEAKPAIRS);
    }
    vector<double> fragmentDistancesTh;
    vector<double> fragmentDistancesPpm;
    for (vector< pair<Peak, Peak> >::const_iterator i = pairedFragmentPeaks.begin();
         i != pairedFragmentPeaks.end();
         i++) {
      double diffTh = i->first.getLocation() - i->second.getLocation();
      fragmentDistancesTh.push_back(diffTh);
      fragmentDistancesPpm.push_back(diffTh * MILLION / i->first.getLocation());
    }
    // estimate the parameters of the component distributions for each of the mixed distributions
    carp(CARP_INFO, "Using %d peak pairs for fragment error estimation.",
         pairedFragmentPeaks.size());
    estimateMuSigma(fragmentDistancesPpm, MIN_SIGMA_PPM,
                    &fragmentMuPpm2Measures, &fragmentSigmaPpm2Measures);
  }

  if (!precursorFailure->empty()) {
    carp(CARP_DEBUG, "Failed precursor! %s", precursorFailure->c_str());
  } else {
    carp(CARP_DEBUG, "precursor_mu_ppm_2measures: %f", precursorMuPpm2Measures);
    carp(CARP_DEBUG, "precursor_sigma_ppm_2measures: %f", precursorSigmaPpm2Measures);
  }

  if (!fragmentFailure->empty()) {
    carp(CARP_DEBUG, "Failed fragment! %s", fragmentFailure->c_str());
  } else {
    carp(CARP_DEBUG, "fragment_mu_ppm_2measures: %f", fragmentMuPpm2Measures);
    carp(CARP_DEBUG, "fragment_sigma_ppm_2measures: %f", fragmentSigmaPpm2Measures);
  }

  // what we have now measured, in the fit Gaussians, is the distribution of the difference
  // of two values drawn from the distribution of error values.
  // Assuming the error values are normally distributed with mean 0 and variance s^2, the
  // differences are normally distributed with mean 0 and variance 2*s^2:
  // http://mathworld.wolfram.com/NormalDifferenceDistribution.html
  // i.e., differences are normally distributed with mean=0 and sd=sqrt(2)*s
  // hence, if differences have sd=diff_sigma, then errors have sd diff_sigma/sqrt(2)
  //
  // incidentally, this transformation doesn't matter one bit, practically, since we're
  // inferring a multiplier for this value empirically. But it lets us report something
  // with an easily-interpretable meaning as an intermediate value
  if (precursorFailure->empty()) {
    *precursorSigmaPpm = precursorSigmaPpm2Measures / sqrt(2);
    // generate prediction by multiplying by empirically-derived value
    *precursorPredictionPpm = PRECURSOR_SIGMA_MULTIPLIER * *precursorSigmaPpm;
  }
  if (fragmentFailure->empty()) {
    *fragmentSigmaPpm = fragmentSigmaPpm2Measures / sqrt(2);
    // generate prediction by multiplying by empirically-derived value
    *fragmentPredictionTh = FRAGMENT_SIGMA_MULTIPLIER * *fragmentSigmaPpm;
  }
}

RunAttributeResult ErrorCalc::summarize() const {
  string precursorFailure, fragmentFailure;
  double precursorSigmaPpm, fragmentSigmaPpm, precursorPredictionPpm, fragmentPredictionTh;
  calcMassErrorDist(&precursorFailure, &fragmentFailure, &precursorSigmaPpm,
    &fragmentSigmaPpm, &precursorPredictionPpm, &fragmentPredictionTh);
  carp(CARP_INFO, "Precursor and fragment error summary:");
  vector<string> messages;
  char msg[64];
  if (precursorFailure.empty()) {
    carp(CARP_INFO, "precursor standard deviation: %f ppm", precursorSigmaPpm);
    sprintf(msg, "Precursor error: %.2f ppm", precursorPredictionPpm);
    messages.push_back(string(msg));
  } else {
    carp(CARP_ERROR, "Precursor error calculation failed: %s", precursorFailure.c_str());
  }
  if (fragmentFailure.empty()) {
    carp(CARP_INFO, "fragment standard deviation: %f ppm", fragmentSigmaPpm);
    sprintf(msg, "Fragment bin size: %.4f Th", fragmentPredictionTh);
    messages.push_back(string(msg));
  } else {
    carp(CARP_ERROR, "Fragment error calculation failed: %s", fragmentFailure.c_str());
  }
  carp(CARP_INFO, "");

  RunAttributeResult result;
  result.setValue(KEY_MESSAGES, StringUtils::Join(messages, '\n'));
  result.setValue(KEY_PRECURSOR_FAILURE, precursorFailure);
  result.setValue(KEY_FRAGMENT_FAILURE, fragmentFailure);
  result.setValue(KEY_PRECURSOR_SIGMA, StringUtils::ToString(precursorSigmaPpm));
  result.setValue(KEY_FRAGMENT_SIGMA, StringUtils::ToString(fragmentSigmaPpm));
  result.setValue(KEY_PRECURSOR_PREDICTION, StringUtils::ToString(precursorPredictionPpm));
  result.setValue(KEY_FRAGMENT_PREDICTION, StringUtils::ToString(fragmentPredictionTh));
  return result;
}

const string ErrorCalc::KEY_MESSAGES             = "search_param_messages";
const string ErrorCalc::KEY_PRECURSOR_FAILURE    = "precursor_failure";
const string ErrorCalc::KEY_FRAGMENT_FAILURE     = "fragment_failure";
const string ErrorCalc::KEY_PRECURSOR_SIGMA      = "precursor_sigma_ppm";
const string ErrorCalc::KEY_FRAGMENT_SIGMA       = "frag_sigma_ppm";
const string ErrorCalc::KEY_PRECURSOR_PREDICTION = "precursor_prediction_ppm";
const string ErrorCalc::KEY_FRAGMENT_PREDICTION  = "fragment_prediction_th";

PerChargeErrorCalc::PerChargeErrorCalc(int charge):
  charge_(charge), chargeForBinSize_(charge > 0 ? charge : 4), numTotalSpectra_(0),
  numPassingSpectra_(0), numSpectraSameBin_(0), numSpectraWithinPpm_(0),
  numSpectraWithinPpmAndScans_(0), numMultipleFragBins_(0), numSingleFragBins_(0) {
  if (!numeric_limits<double>::is_iec559) {
    carp(CARP_FATAL, "Something went wrong.");
  }
}

PerChargeErrorCalc::~PerChargeErrorCalc() {
}

void PerChargeErrorCalc::processSpectrum(
  const Spectrum* spectrum,
  const vector<double>& binnedSpectrum
) {
  ++numTotalSpectra_;

  if (spectrum->getNumPeaks() < Params::GetInt("pm-min-scan-frag-peaks")) {
    return;
  }

  bool chargeOk = false;
  const vector<SpectrumZState>& zStates = spectrum->getZStates();
  for (vector<SpectrumZState>::const_iterator i = zStates.begin(); i != zStates.end(); i++) {
    if (i->getCharge() == charge_) {
      chargeOk = true;
      break;
    }
  }
  if (!chargeOk) {
    return;
  }

  double precursorMz = getPrecursorMz(spectrum);
  if (std::isnan(precursorMz) ||
      !(Params::GetDouble("pm-min-precursor-mz") <= precursorMz &&
        precursorMz <= Params::GetDouble("pm-max-precursor-mz"))) {
    return;
  }

  ++numPassingSpectra_;
  // pull out the top fragments by intensity
  vector<Peak> peaks = spectrum->getPeaks();
  std::sort(peaks.begin(), peaks.end(), Peak::compareByIntensity);
  peaks.resize(Params::GetInt("pm-top-n-frag-peaks"));

  int precursorBinIndex = getBinIndexPrecursor(precursorMz);
  map< int, pair< const Spectrum*, vector<Peak> > >::const_iterator prevIter =
    spectra_.find(precursorBinIndex);
  if (prevIter != spectra_.end()) {
    // there was a previous spectrum in this bin; check to see if they're a pair
    const Spectrum* prev = prevIter->second.first;
    const double precursorMzPrev = getPrecursorMz(prev);
    const double precursorMzDiffPpm = (precursorMz - precursorMzPrev) * MILLION / precursorMz;
    ++numSpectraSameBin_;
    // check precursor
    if (abs(precursorMzDiffPpm) <= Params::GetDouble("pm-max-precursor-delta-ppm")) {
      // check scan count between the scans
      ++numSpectraWithinPpm_;
      if (abs(spectrum->getFirstScan() - prev->getFirstScan()) <= Params::GetInt("pm-max-scan-separation")) {
        // count the fragment peaks in common
        ++numSpectraWithinPpmAndScans_;
        vector< pair<Peak, Peak> > pairedFragments =
          pairFragments(prevIter->second.second, peaks);
        if (pairedFragments.size() >= Params::GetInt("pm-min-common-frag-peaks")) {
          // we've got a pair! record everything
          sort(pairedFragments.begin(), pairedFragments.end(), sortPairedFragments);
          vector< pair<Peak, Peak> >::const_iterator stop =
            pairedFragments.size() >= Params::GetInt("pm-pair-top-n-frag-peaks")
              ? pairedFragments.begin() + Params::GetInt("pm-pair-top-n-frag-peaks")
              : pairedFragments.end();
          for (vector< pair<Peak, Peak> >::const_iterator i = pairedFragments.begin();
              i != stop;
              i++) {
            pairedFragmentPeaks_.push_back(make_pair(i->first, i->second));
          }
          pairedPrecursorMzs_.push_back(make_pair(precursorMzPrev, precursorMz));
        }
      }
    }
  }
  // make the new spectrum its bin's representative
  spectra_[precursorBinIndex] = make_pair(spectrum, peaks);
}

void PerChargeErrorCalc::clearBins() {
  spectra_.clear();
}

void PerChargeErrorCalc::nextFile() {
  clearBins();
}

int PerChargeErrorCalc::getNumPassingSpectra() const{
  return numPassingSpectra_;
}

int PerChargeErrorCalc::getNumSpectraSameBin() const{
  return numSpectraSameBin_;
}

int PerChargeErrorCalc::getNumSpectraWithinPpm() const{
  return numSpectraWithinPpm_;
}

int PerChargeErrorCalc::getNumSpectraWithinPpmAndScans() const{
  return numSpectraWithinPpmAndScans_;
}

int PerChargeErrorCalc::getNumMultipleFragBins() const {
  return numMultipleFragBins_;
}

int PerChargeErrorCalc::getNumSingleFragBins() const {
  return numSingleFragBins_;
}

const vector< pair<Peak, Peak> >& PerChargeErrorCalc::getPairedFragmentPeaks() const {
  return pairedFragmentPeaks_;
}

const vector< pair<double, double> >& PerChargeErrorCalc::getPairedPrecursorMzs() const {
  return pairedPrecursorMzs_;
}

int PerChargeErrorCalc::getBinIndexPrecursor(double mz) {
  static double lowestPrecursorBinStartMz = Params::GetDouble("pm-min-precursor-mz") -
    fmod(Params::GetDouble("pm-min-precursor-mz"), AVERAGINE_PEAK_SEPARATION / chargeForBinSize_);
  return (int)((mz - lowestPrecursorBinStartMz) / (AVERAGINE_PEAK_SEPARATION / chargeForBinSize_));
}

int PerChargeErrorCalc::getBinIndexFragment(double mz) {
  static double lowestFragmentBinStartMz = Params::GetDouble("pm-min-frag-mz") -
    fmod(Params::GetDouble("pm-min-frag-mz"), AVERAGINE_PEAK_SEPARATION);
  return (int)((mz - lowestFragmentBinStartMz) / AVERAGINE_PEAK_SEPARATION);
}

double PerChargeErrorCalc::getPrecursorMz(const Spectrum* spectrum) const {
  const vector<SpectrumZState>& zStates = spectrum->getZStates();
  for (vector<SpectrumZState>::const_iterator i = zStates.begin(); i != zStates.end(); i++) {
    if (i->getCharge() == chargeForBinSize_) {
      return i->getMZ();
    }
  }
  return numeric_limits<double>::quiet_NaN();;
}

vector< pair<Peak, Peak> > PerChargeErrorCalc::pairFragments(
  const vector<Peak>& prev,
  const vector<Peak>& cur
) {
  map<int, Peak> mapPrev = binFragments(prev);
  map<int, Peak> mapCur = binFragments(cur);
  vector< pair<Peak, Peak> > pairs;
  for (map<int, Peak>::const_iterator i = mapPrev.begin(); i != mapPrev.end(); i++) {
    map<int, Peak>::const_iterator j = mapCur.find(i->first);
    if (j != mapCur.end()) {
      pairs.push_back(make_pair(i->second, j->second));
    }
  }
  return pairs;
}

map<int, Peak> PerChargeErrorCalc::binFragments(const vector<Peak>& peaks) {
  map<int, Peak> binFragmentMap;
  set<int> binsToRemove;
  for (vector<Peak>::const_iterator i = peaks.begin(); i != peaks.end(); i++) {
    FLOAT_T mz = i->getLocation();
    if (mz < Params::GetDouble("pm-min-frag-mz")) {
      continue;
    }
    int binIndex = getBinIndexFragment(mz);
    if (binFragmentMap.find(binIndex) != binFragmentMap.end()) {
      binsToRemove.insert(binIndex);
    } else {
      binFragmentMap[binIndex] = *i;
    }
  }
  for (set<int>::const_iterator i = binsToRemove.begin(); i != binsToRemove.end(); i++) {
    binFragmentMap.erase(*i);
  }
  numMultipleFragBins_ += binsToRemove.size();
  numSingleFragBins_ += binFragmentMap.size();
  return binFragmentMap;
}

bool PerChargeErrorCalc::sortPairedFragments(
  const pair<Peak, Peak>& x,
  const pair<Peak, Peak>& y
) {
  return min(x.first.getIntensity(), x.second.getIntensity()) <
         min(y.first.getIntensity(), y.second.getIntensity());
}

Model::Model(double nMean, double nStd, double nMinStd, double uStart, double uEnd):
  normal_(NormalDistribution(nMean, nStd, nMinStd)), uniform_(UniformDistribution(uStart, uEnd)) {
  weights_[0] = weights_[1] = log(0.5);
  summaries_[0] = summaries_[1] = 0;
}

Model::~Model() {
}

double Model::fit(const vector<double>& data) {
  const double stopThreshold = 0.1;
  const int maxIterations = 1e8;

  double initialLogProbSum = -numeric_limits<double>::infinity();
  double lastLogProbSum;
  int i = 0;
  double improvement = numeric_limits<double>::infinity();
  while (improvement > 0.1 && i < maxIterations + 1) {
    fromSummaries();
    double logProbSum = summarize(data);

    if (i++ == 0) {
      initialLogProbSum = logProbSum;
    } else {
      improvement = logProbSum - lastLogProbSum;
      carp(CARP_DETAILED_DEBUG, "Improvement: %f", improvement);
    }
    lastLogProbSum = logProbSum;
  }
  clearSummaries();
  return lastLogProbSum - initialLogProbSum;
}

double Model::getMu() const {
  return normal_.getMu();
}

double Model::getSigma() const {
  return normal_.getSigma();
}

double Model::summarize(const vector<double>& x) {
  vector<double> r;
  r.reserve(x.size() * 2);
  normal_.logProbability(x, &r);
  uniform_.logProbability(x, &r);

  double logProbSum = 0;
  for (size_t i = 0; i < x.size(); i++) {
    double total = -numeric_limits<double>::infinity();

    for (size_t j = 0; j < 2; j++) {
      r[j * x.size() + i] += weights_[j];
      total = pairLse(total, r[j * x.size() + i]);
    }

    for (size_t j = 0; j < 2; j++) {
      r[j * x.size() + i] = exp(r[j * x.size() + i] - total);
      summaries_[j] += r[j * x.size() + i];
    }

    logProbSum += total;
  }

  normal_.summarize(x, &r[0]);
  uniform_.summarize(x, &r[x.size()]);
  return logProbSum;
}

void Model::fromSummaries() {
  double sum = summaries_[0] + summaries_[1];
  if (sum == 0) {
    return;
  }

  summaries_[0] /= sum;
  summaries_[1] /= sum;

  normal_.fromSummaries();
  weights_[0] = log(summaries_[0]);
  summaries_[0] = 0;

  uniform_.fromSummaries();
  weights_[1] = log(summaries_[1]);
  summaries_[1] = 0;
}

void Model::clearSummaries() {
  summaries_[0] = summaries_[1] = 0;
  normal_.clearSummaries();
  uniform_.clearSummaries();
}

double Model::pairLse(double x, double y) {
  const double inf = numeric_limits<double>::infinity();
  if (x == inf || y == inf) {
    return inf;
  } else if (x == -inf) {
    return y;
  } else if (y == -inf) {
    return x;
  } else if (x > y) {
    return x + log(exp(y - x) + 1);
  }
  return y + log(exp(x - y) + 1);
}

Model::NormalDistribution::NormalDistribution(double mean, double std, double minStd):
  mu_(mean), sigma_(std), minStd_(minStd),
  logSigmaSqrt2Pi_(-log(std * SQRT_2_PI)), twoSigmaSquared_(2 * pow(std, 2)) {
  clearSummaries();
}

Model::NormalDistribution::~NormalDistribution() {
}

double Model::NormalDistribution::getMu() const {
  return mu_;
}

double Model::NormalDistribution::getSigma() const {
  return sigma_;
}

void Model::NormalDistribution::logProbability(const vector<double>& x, vector<double>* r) const {
  for (vector<double>::const_iterator i = x.begin(); i != x.end(); i++) {
    r->push_back(logSigmaSqrt2Pi_ - pow(*i - mu_, 2) / twoSigmaSquared_);
  }
}

void Model::NormalDistribution::summarize(const vector<double>& x, double* weights) {
  for (size_t i = 0; i < x.size(); i++) {
    summaries_[0] += weights[i];
    summaries_[1] += weights[i] * x[i];
    summaries_[2] += weights[i] * pow(x[i], 2);
  }
}

void Model::NormalDistribution::fromSummaries() {
  if (summaries_[0] == 0) {
    return;
  }
  mu_ = summaries_[1] / summaries_[0];
  sigma_ = sqrt(summaries_[2] / summaries_[0] - pow(summaries_[1], 2) / pow(summaries_[0], 2));
  if (sigma_ < minStd_) {
    sigma_ = minStd_;
  }
  clearSummaries();
  logSigmaSqrt2Pi_ = -log(sigma_ * SQRT_2_PI);
  twoSigmaSquared_ = 2 * pow(sigma_, 2);
}

void Model::NormalDistribution::clearSummaries() {
  summaries_[0] = summaries_[1] = summaries_[2] = 0;
}

Model::UniformDistribution::UniformDistribution(double start, double end):
  start_(start), end_(end), logP_(-log(end - start)) {
  clearSummaries();
}

Model::UniformDistribution::~UniformDistribution() {
}

void Model::UniformDistribution::logProbability(const vector<double>& x, vector<double>* r) const {
  for (vector<double>::const_iterator i = x.begin(); i != x.end(); i++) {
    r->push_back(start_ <= *i && *i <= end_ ? logP_ : -numeric_limits<double>::infinity());
  }
}

void Model::UniformDistribution::summarize(const vector<double>& x, double* weights) {
  for (size_t i = 0; i < x.size(); i++) {
    if (weights[i] <= 0) {
      continue;
    }
    double value = x[i];
    if (value < summaries_[0]) {
      summaries_[0] = value;
    }
    if (value > summaries_[1]) {
      summaries_[1] = value;
    }
  }
}

void Model::UniformDistribution::fromSummaries() {
  start_ = summaries_[0];
  end_ = summaries_[1];
  logP_ = -log(end_ - start_);
  clearSummaries();
}

void Model::UniformDistribution::clearSummaries() {
  summaries_[0] = numeric_limits<double>::infinity();
  summaries_[1] = -numeric_limits<double>::infinity();
}

const double BINNING_MIN_MZ = 50.5*AVERAGINE_PEAK_SEPARATION;
const double BINNING_MAX_MZ = 6000.5*AVERAGINE_PEAK_SEPARATION;

void estimateMuSigma(
  const vector<double>& data,
  double minSigma,
  double* muFit,
  double* sigmaFit
) {
  double dataMin = data[0];
  double dataMax = data[0];
  double muMixedDist = data[0];
  for (vector<double>::const_iterator i = data.begin() + 1; i != data.end(); i++) {
    muMixedDist += *i;
    if (*i < dataMin) {
      dataMin = *i;
    }
    if (*i > dataMax) {
      dataMax = *i;
    }
  }
  muMixedDist /= data.size();

  double sigmaMixedDist = 0;
  for (vector<double>::const_iterator i = data.begin(); i != data.end(); i++) {
    sigmaMixedDist += pow(*i - muMixedDist, 2);
  }
  sigmaMixedDist = sqrt(sigmaMixedDist / data.size());

  carp(CARP_DEBUG, "mixed distribution: min %f, max %f, mean %f, sd %f",
       dataMin, dataMax, muMixedDist, sigmaMixedDist);

  // model the observed distribution as a mixture of Gaussian and uniform
  Model model(muMixedDist, sigmaMixedDist, minSigma, dataMin, dataMax);
  // fit the mixture model with EM
  double improvement = model.fit(data);
  carp(CARP_DEBUG, "model improvement: %f", improvement);

  *muFit = model.getMu();
  *sigmaFit = model.getSigma();
  carp(CARP_DEBUG, "fit: mean=%f, sigma=%f", *muFit, *sigmaFit);
}

int calcBinIndexMassPrecursor(double mass) {
  return calcBinIndexMzFragment(mass);
}

int calcBinIndexMzFragment(double mz) {
  double relativeMz = mz - BINNING_MIN_MZ;
  int bin = (int)floor(relativeMz / AVERAGINE_PEAK_SEPARATION);
  return bin > 0 ? bin : 0;
}

double calcMH(double mz, int charge) {
  return (mz - MASS_H_AVERAGE) * charge + MASS_H_AVERAGE;
}

vector<double> binSpectrum(const Spectrum* spectrum) {
  //int numBins = (int)((BINNING_MAX_MZ - BINNING_MIN_MZ)/AVERAGINE_PEAK_SEPARATION) + 1;
  const int numBins = 5950;
  vector<Peak> peaks = spectrum->getPeaks();
  vector<double> matrix(numBins, 0);
  for (vector<Peak>::const_iterator i = peaks.begin(); i != peaks.end(); i++) {
    FLOAT_T mz = i->getLocation();
    FLOAT_T intensity = i->getIntensity();
    if (mz < BINNING_MIN_MZ || mz > BINNING_MAX_MZ) {
      continue;
    }
    int binIdx = calcBinIndexMzFragment(mz);
    if (binIdx < 0 || binIdx >= numBins) {
      continue;
    }
    if (intensity > matrix[binIdx]) {
      matrix[binIdx] = intensity;
    }
  }
  double sum = 0;
  for (vector<double>::const_iterator i = matrix.begin(); i != matrix.end(); i++) {
    sum += *i;
  }
  if (sum > 0) {
    for (vector<double>::iterator i = matrix.begin(); i != matrix.end(); i++) {
      *i /= sum;
    }
  }
  return matrix;
}

int processSpectra(const vector<string>& files, vector<RunAttributeDetector*> detectors) {
  int n = 0;
  for (vector<string>::const_iterator i = files.begin(); i != files.end(); i++) {
    carp(CARP_INFO, "param-medic processing input file %s...", i->c_str());
    if (i > files.begin()) {
      for (vector<RunAttributeDetector*>::const_iterator j = detectors.begin();
           j != detectors.end();
           j++) {
        (*j)->nextFile();
      }
    }
    SpectrumCollection* collection = SpectrumCollectionFactory::create(*i);
    collection->parse();
    for (SpectrumIterator j = collection->begin(); j != collection->end(); j++) {
      if ((*j)->getNumPeaks() < Params::GetInt("pm-min-scan-frag-peaks")) {
        continue;
      }
      // bin the spectrum peaks
      vector<double> binned = ParamMedic::binSpectrum(*j);
      // run each of the detectors on the binned peaks
      for (vector<RunAttributeDetector*>::const_iterator k = detectors.begin();
           k != detectors.end();
           k++) {
        (*k)->processSpectrum(*j, binned);
      }
      n++;
    }
    delete collection;
  }
  return n;
}

int processSpectra(const vector<string>& files, RunAttributeDetector* detector) {
  vector<RunAttributeDetector*> detectors(1, detector);
  return processSpectra(files, detectors);
}

}

