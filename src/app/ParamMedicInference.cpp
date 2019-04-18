#include "io/carp.h"
#include "util/AminoAcidUtil.h"
#include "util/mass.h"
#include "util/StringUtils.h"
#include "ParamMedicApplication.h"
#include "ParamMedicInference.h"

using namespace std;

namespace ParamMedic {

const double AVERAGINE_PEAK_SEPARATION = 1.0005079;

// Constants used by PhosphoLossProportionCalc
const double SEARCH_MOD_MASS_PHOSPHO = 79.966331;
const double DELTA_MASS_PHOSPHO_LOSS = 98.0;
const double PHOSPHO_ZSCORE_CUTOFF = 9.0;

PhosphoLossProportionCalc::PhosphoLossProportionCalc()
: sumProportionsInPhosphoLoss_(0), numSpectraUsed_(0) {
  sumsProportionsPerControlPeak_.insert(make_pair(-20, 0));
  sumsProportionsPerControlPeak_.insert(make_pair(-15, 0));
  sumsProportionsPerControlPeak_.insert(make_pair(-12, 0));
  sumsProportionsPerControlPeak_.insert(make_pair(-10, 0));
  sumsProportionsPerControlPeak_.insert(make_pair(10, 0));
  sumsProportionsPerControlPeak_.insert(make_pair(12, 0));
  sumsProportionsPerControlPeak_.insert(make_pair(15, 0));
  sumsProportionsPerControlPeak_.insert(make_pair(20, 0));
}

PhosphoLossProportionCalc::~PhosphoLossProportionCalc() {
}

void PhosphoLossProportionCalc::processSpectrum(
  const Crux::Spectrum* spectrum,
  const vector<double>& binnedSpectrum
) {
  const vector<SpectrumZState>& zStates = spectrum->getZStates();
  if (zStates.size() != 1) {
    return;
  }
  numSpectraUsed_++;
  double mz = spectrum->getPrecursorMz() - (DELTA_MASS_PHOSPHO_LOSS / zStates.front().getCharge());
  int bin = calcBinIndexMzFragment(mz);
  for (map<int, double>::iterator i = sumsProportionsPerControlPeak_.begin();
       i != sumsProportionsPerControlPeak_.end();
       i++) {
    i->second += binnedSpectrum[bin + i->first];
  }
  sumProportionsInPhosphoLoss_ += binnedSpectrum[bin];
}

RunAttributeResult PhosphoLossProportionCalc::summarize() const {
  RunAttributeResult result;
  if (numSpectraUsed_ == 0) {
    carp(CARP_WARNING, "No spectra usable for phosphorylation detection");
    result.setValue("phospho_present", RunAttributeResult::ERROR);
    result.setValue("phospho_statistic", RunAttributeResult::ERROR);
    return result;
  }
  vector<double> values;
  for (map<int, double>::const_iterator i = sumsProportionsPerControlPeak_.begin();
       i != sumsProportionsPerControlPeak_.end();
       i++) {
    values.push_back(i->second);
  }
  double mean = MathUtil::Mean(values);
  double sd = MathUtil::StdDev(values);
  carp(CARP_DEBUG, "Phospho control peaks (mean=%.3f):", mean);
  for (map<int, double>::const_iterator i = sumsProportionsPerControlPeak_.begin();
       i != sumsProportionsPerControlPeak_.end();
       i++) {
    carp(CARP_DEBUG, "    %d: %.4f", i->first, i->second);
  }
  double proportionToControl = sumProportionsInPhosphoLoss_ / mean;
  double zScoreToControl = (sumProportionsInPhosphoLoss_ - mean) / sd;
  carp(CARP_DEBUG, "Phospho-loss peak: %.3f", sumProportionsInPhosphoLoss_);
  carp(CARP_DEBUG, "Phospho: ratio phospho-loss to control peaks: %.5f (z=%.3f)",
       proportionToControl, zScoreToControl);

  // summarize results
  bool phospho = false;
  if (zScoreToControl > PHOSPHO_ZSCORE_CUTOFF) {
    carp(CARP_INFO, "Phosphorylation: detected");
    result.addMod(Modification::LOCATION_CTERM, SEARCH_MOD_MASS_PHOSPHO, true);
    phospho = true;
  } else {
    carp(CARP_INFO, "Phosphorylation: not detected");
  }
  result.setValue("phospho_present", phospho ? "T" : "F");
  result.setValue("phospho_statistic", StringUtils::ToString(zScoreToControl));
  return result;
}

// Constants used by Tmt6vs10Detector
const int TMT610_MIN_PEAKS_TMT_PRESENT_FOR_DECISION = 2;
const double TMT610_MIN_TMT10_PROPORTION_FOR_DECISION = 0.2;
const double TMT10_MASS_DIFF_FROM_6 = 0.00632;
const double TMT10_PADDING_EACH_SIDE = TMT10_MASS_DIFF_FROM_6 * 5;
const double TMT610_WINDOW_WIDTH = TMT10_PADDING_EACH_SIDE * 2 + TMT10_MASS_DIFF_FROM_6;
const double TMT610_PEAK_WIDTH_FOR_DETECT = 0.003;
const int TMT610_SMALLEST_NOMINAL = 127;
const double TMT610_SMALLEST_BIN_CENTER = TMT610_SMALLEST_NOMINAL * AVERAGINE_PEAK_SEPARATION;
const double TMT610_SMALLEST_BIN_MINMASS = TMT610_SMALLEST_BIN_CENTER - AVERAGINE_PEAK_SEPARATION / 2;
const int TMT610_MINPEAKS_FOR_DETECT = 30;

Tmt6vs10Detector::Tmt6vs10Detector() {
  nominalMassTmt6Mass_[127] = 127.12476;
  nominalMassTmt6Mass_[128] = 128.13443;
  nominalMassTmt6Mass_[129] = 129.13147;
  nominalMassTmt6Mass_[130] = 130.14114;
  nominalMassTmt10Mass_[127] = 127.13108;
  nominalMassTmt10Mass_[128] = 128.12811;
  nominalMassTmt10Mass_[129] = 129.13779;
  nominalMassTmt10Mass_[130] = 130.13482;
  for (map<int, double>::const_iterator i = nominalMassTmt6Mass_.begin();
       i != nominalMassTmt6Mass_.end();
       i++) {
    nominalMassMinBinMass_[i->first] =
      min(nominalMassTmt10Mass_[i->first], i->second ) - TMT10_PADDING_EACH_SIDE;
  }
}

Tmt6vs10Detector::~Tmt6vs10Detector() {
}

void Tmt6vs10Detector::processSpectrum(
  const Crux::Spectrum* spectrum,
  const vector<double>& binnedSpectrum
) {
  const vector<SpectrumZState>& zStates = spectrum->getZStates();
  if (zStates.size() != 1) {
    return;
  }

  // find the smallest peak in this scan that might be in one of the double-reagent bins
  vector<Peak> peaks = spectrum->getPeaks();
  std::sort(peaks.begin(), peaks.end(), Peak::compareByMZ);
  size_t curIdx = 0;
  for ( ; curIdx < peaks.size(); curIdx++) {
    if (peaks[curIdx].getLocation() >= TMT610_SMALLEST_BIN_MINMASS) {
      break;
    }
  }
  const double tmt610BiggestBinMaxMass = nominalMassMinBinMass_[130] + TMT610_WINDOW_WIDTH;
  for ( ; curIdx < peaks.size(); curIdx++) {
    double curMz = peaks[curIdx].getLocation();
    if (curMz >= tmt610BiggestBinMaxMass) {
      break;
    }

    // int() does floor. Add 1 because 0-index vs. 1-index
    int curBin = int((curMz - AVERAGINE_PEAK_SEPARATION / 2) / AVERAGINE_PEAK_SEPARATION) + 1;
    if (nominalMassMinBinMass_[curBin] < curMz &&
        curMz < nominalMassMinBinMass_[curBin] + TMT610_WINDOW_WIDTH) {
      map< int, vector<double> >::iterator i = nominalMassAllPeaks_.find(curBin);
      if (i == nominalMassAllPeaks_.end()) {
        nominalMassAllPeaks_[curBin] = vector<double>(1, curMz);
      } else  {
        i->second.push_back(curMz);
      }
    }
  }
}

RunAttributeResult Tmt6vs10Detector::summarize() const {
  int nPeaksWithEnoughTmt10 = 0;
  for (map< int, vector<double> >::const_iterator i = nominalMassAllPeaks_.begin();
       i != nominalMassAllPeaks_.end();
       i++) {
    const vector<double>& curPeaks = i->second;
    carp(CARP_DEBUG, "nominal mass %d: %d", i->first, curPeaks.size());
    if (curPeaks.size() >= TMT610_MINPEAKS_FOR_DETECT) {
      double smallestMzThisBin = i->first * AVERAGINE_PEAK_SEPARATION - AVERAGINE_PEAK_SEPARATION / 2;
      carp(CARP_DEBUG, "    %f-%f", curPeaks.front(), curPeaks.back());
      // find the proportion of peaks falling in the expected range for each TMT ion at this nominal mass
      map<int, double>::const_iterator lookup;
      lookup = nominalMassTmt6Mass_.find(i->first);
      double peak6Mz = (lookup != nominalMassTmt6Mass_.end()) ?
        lookup->second : numeric_limits<double>::quiet_NaN();
      lookup = nominalMassTmt10Mass_.find(i->first);
      double peak10Mz = (lookup != nominalMassTmt10Mass_.end()) ?
        lookup->second : numeric_limits<double>::quiet_NaN();
      double minMassPeak6 = peak6Mz - TMT610_PEAK_WIDTH_FOR_DETECT / 2;
      double maxMassPeak6 = peak6Mz + TMT610_PEAK_WIDTH_FOR_DETECT / 2;
      double minMassPeak10 = peak10Mz - TMT610_PEAK_WIDTH_FOR_DETECT / 2;
      double maxMassPeak10 = peak10Mz + TMT610_PEAK_WIDTH_FOR_DETECT / 2;
      int nPeaksNearPeak6 = 0, nPeaksNearPeak10 = 0;
      for (vector<double>::const_iterator j = curPeaks.begin(); j != curPeaks.end(); j++) {
        if (minMassPeak6 <= *j && *j <= maxMassPeak6) {
          nPeaksNearPeak6++;
        } else if (minMassPeak10 <= *j && *j <= maxMassPeak10) {
          nPeaksNearPeak10++;
        }
      }
      double proportionNearPeak6 = (double)nPeaksNearPeak6 / curPeaks.size();
      double proportionNearPeak10 = (double)nPeaksNearPeak10 / curPeaks.size();
      carp(CARP_DEBUG, "Near 6: %f. Near 10: %f", proportionNearPeak6, proportionNearPeak10);
      if (proportionNearPeak10 > TMT610_MIN_TMT10_PROPORTION_FOR_DECISION) {
        carp(CARP_DEBUG, "Peak %d has TMT10 signal", i->first);
        nPeaksWithEnoughTmt10++;
      }
    } else {
      carp(CARP_DEBUG, "Too few, failing peak.");
    }
  }
  RunAttributeResult result;
  result.setValue("n_peaks_with_enough_TMT10_signal",
                  StringUtils::ToString(nPeaksWithEnoughTmt10));
  if (nPeaksWithEnoughTmt10 >= TMT610_MIN_PEAKS_TMT_PRESENT_FOR_DECISION) {
    carp(CARP_DEBUG, "Detected TMT 10-plex! %d peaks had signal.",
         nPeaksWithEnoughTmt10);
    result.setValue("TMT10", "T");
    return result;
  }
  carp(CARP_DEBUG, "Did NOT detect TMT 10-plex. Only %d peaks had signal.",
       nPeaksWithEnoughTmt10);
  result.setValue("TMT10", "F");
  return result;
}

// Constants used by ReporterIonProportionCalc
const double SEARCH_MOD_MASS_ITRAQ_4PLEX = 144.10253;
const double SEARCH_MOD_MASS_ITRAQ_8PLEX = 304.2022;
const double SEARCH_MOD_MASS_TMT_2PLEX = 225.155833;
const double SEARCH_MOD_MASS_TMT_6PLEX = 229.162932;
const double TMT_REPORTER_ION_ZSCORE_CUTOFF = 3.5;

ReporterIonProportionCalc::ReporterIonProportionCalc()
: foundMs3Scans_(false) {
  tmt6PlexOnlyReporterIonMzs_.push_back(128.1344);
  tmt6PlexOnlyReporterIonMzs_.push_back(129.1378);
  tmt6PlexOnlyReporterIonMzs_.push_back(130.1411);
  tmt6PlexOnlyReporterIonMzs_.push_back(131.1382);

  REPORTER_ION_TYPE types[] = { TMT_2PLEX, TMT_6PLEX, ITRAQ_4PLEX, ITRAQ_8PLEX, CONTROL };
  for (size_t i = 0; i < sizeof(types) / sizeof(REPORTER_ION_TYPE); i++) {
    REPORTER_ION_TYPE t = types[i];
    reporterIonTypeBins_[t] = vector<double>();
    reporterTypeBinSumProportion_[t] = map<int, double>();
    vector<double> mzs;
    switch (t) {
    case TMT_2PLEX:
      mzs.push_back(126.1277);
      mzs.push_back(127.1311);
      break;
    case TMT_6PLEX:
      mzs = tmt6PlexOnlyReporterIonMzs_;
      break;
    case ITRAQ_4PLEX:
      mzs.push_back(114.0);
      mzs.push_back(115.0);
      mzs.push_back(116.0);
      mzs.push_back(117.0);
      break;
    case ITRAQ_8PLEX:
      mzs.push_back(117.0);
      mzs.push_back(118.0);
      mzs.push_back(119.0);
      mzs.push_back(121.0);
      break;
    case CONTROL:
      mzs.push_back(111.0);
      mzs.push_back(112.0);
      mzs.push_back(120.0);
      mzs.push_back(122.0);
      mzs.push_back(123.0);
      mzs.push_back(124.0);
      mzs.push_back(125.0);
      mzs.push_back(133.0);
      mzs.push_back(134.0);
      break;
    }
    for (vector<double>::const_iterator j = mzs.begin(); j != mzs.end(); j++) {
      int bin = calcBinIndexMzFragment(*j);
      reporterIonTypeBins_[t].push_back(bin);
      reporterTypeBinSumProportion_[t][bin] = 0.0;
    }
  }

  reporterIonTStatThresholds_[TMT_2PLEX] = 2.0;
  reporterIonTStatThresholds_[TMT_6PLEX] = 2.0;
  reporterIonTStatThresholds_[ITRAQ_4PLEX] = 1.5;
  reporterIonTStatThresholds_[ITRAQ_8PLEX] = 1.5;

  carp(CARP_DEBUG, "Reporter ion type count (including control): %d",
       reporterIonTypeBins_.size());
}

ReporterIonProportionCalc::~ReporterIonProportionCalc() {
}

void ReporterIonProportionCalc::processSpectrum(
  const Crux::Spectrum* spectrum,
  const vector<double>& binnedSpectrum
) {
  const vector<SpectrumZState>& zStates = spectrum->getZStates();
  if (zStates.size() != 1) {
    return;
  /*} else if ( ? ) {
    foundMs3Scans_ = true;*/
  }
  for (map< REPORTER_ION_TYPE, vector<double> >::const_iterator i = reporterIonTypeBins_.begin();
       i != reporterIonTypeBins_.end();
       i++) {
    map<int, double>& m = reporterTypeBinSumProportion_[i->first];
    for (map<int, double>::const_iterator j = m.begin(); j != m.end(); j++) {
      m[j->first] += binnedSpectrum[j->first];
    }
  }
  tmt10Detector_.processSpectrum(spectrum, binnedSpectrum);
}

RunAttributeResult ReporterIonProportionCalc::summarize() const {
  // summarize the control bins. These sums are used in the t-statistic calculation
  const map<int, double>& mControl = reporterTypeBinSumProportion_.find(CONTROL)->second;
  vector<double> controlBinSums;
  carp(CARP_DEBUG, "Reporter ion control bin sums:");
  for (map<int, double>::const_iterator i = mControl.begin(); i != mControl.end(); i++) {
    // adding 1 to bin number to convert from zero-based index
    carp(CARP_DEBUG, "  %d: %.2f", i->first + 1, i->second);
    controlBinSums.push_back(i->second);
  }
  double controlBinMean = MathUtil::Mean(controlBinSums);
  double controlBinSd = MathUtil::StdDev(controlBinSums);
  carp(CARP_DEBUG, "Reporter ion control bin mean: %.2f", controlBinMean);

  // check each reporter ion type, determine which are significantly elevated
  set<REPORTER_ION_TYPE> significant;
  map<REPORTER_ION_TYPE, double> reporterTypeTStatistics;

  vector<double> tmt6ReporterZScores;

  for (map< REPORTER_ION_TYPE, map<int, double> >::const_iterator i = reporterTypeBinSumProportion_.begin();
       i != reporterTypeBinSumProportion_.end();
       i++) {
    if (i->first == CONTROL) {
      continue;
    }
    // reporter bin sums are used in t-statistic collection
    const map<int, double>& m = i->second;
    vector<double> sums, ionZScores;
    for (map<int, double>::const_iterator j = m.begin(); j != m.end(); j++) {
      double sum = j->second;
      sums.push_back(j->second);
      // z-score is only calculated for debug output purposes. Not used in significance determination
      double zScore = (sum - controlBinMean) / controlBinSd;
      // adding 1 to bin number to convert from zero-based index
      carp(CARP_DEBUG, "    %d, mz=%d: ratio=%.2f, zscore=%.2f",
           i->first, j->first + 1, sum / controlBinMean, zScore);
      ionZScores.push_back(zScore);
    }
    // reporter bin mean is only for debug output
    double mean = MathUtil::Mean(sums);
    carp(CARP_DEBUG, "%d, ion zscores: %s", i->first, StringUtils::Join(ionZScores, '\t').c_str());
    if (i->first == TMT_6PLEX) {
      tmt6ReporterZScores = ionZScores;
    }
    carp(CARP_DEBUG, "%d bin mean: %.2f", i->first, mean);
    // calculate t-statistic of the reporter bin intensity sums vs. the control bin intensity sums
    double tStatistic = tTestInd(sums, controlBinSums, false);
    reporterTypeTStatistics[i->first] = tStatistic;
    // ratio is only calculated for debug output
    double ratio = mean / controlBinMean;
    carp(CARP_DEBUG, "%d, overall: reporter/control mean ratio: %.4f. t-statistic: %.4f",
         i->first, ratio, tStatistic);

    // check the t-statistic against the appropriate threshold and conditionally declare significance
    map<REPORTER_ION_TYPE, double>::const_iterator lookup =
      reporterIonTStatThresholds_.find(i->first);
    if (tStatistic > lookup->second) {
      significant.insert(i->first);
    }
  }

  // Now, all the sets or reporter ions are analyzed and tested.
  // Create the properly-formatted result that summarizes everything.
  RunAttributeResult result;

  // handle iTRAQ
  bool itraq8 = significant.find(ITRAQ_8PLEX) != significant.end();
  bool itraq4 = significant.find(ITRAQ_4PLEX) != significant.end();
  if (itraq8) {
    carp(CARP_INFO, "iTRAQ: 8-plex reporter ions detected");
    result.addMod("K", SEARCH_MOD_MASS_ITRAQ_8PLEX, true);
    result.addMod(Modification::LOCATION_NTERM, SEARCH_MOD_MASS_ITRAQ_8PLEX, true);
    if (!itraq4) {
      carp(CARP_WARNING, "    No iTRAQ 4-plex reporters detected, only 8-plex.");
    }
    itraq4 = false;
  } else if (itraq4) {
    carp(CARP_INFO, "iTRAQ: 4-plex reporter ions detected");
    // 8plex mass same as 4plex, more or less
    result.addMod("K", SEARCH_MOD_MASS_ITRAQ_4PLEX, true);
    result.addMod(Modification::LOCATION_NTERM, SEARCH_MOD_MASS_ITRAQ_4PLEX, true);
    itraq8 = false;
  } else {
    carp(CARP_INFO, "iTRAQ: no reporter ions detected");
  }

  // handle TMT

  bool tmt2 = significant.find(TMT_2PLEX) != significant.end();
  bool tmt6 = significant.find(TMT_6PLEX) != significant.end();

  // special handling for TMT 6 vs. 2: if at least one 6plex ion is present, say we've got 6plex
  vector<double> tmt6IonsPresent;
  if (tmt2 && !tmt6) {
    for (size_t i = 0; i < tmt6ReporterZScores.size(); i++) {
      if (tmt6ReporterZScores[i] > TMT_REPORTER_ION_ZSCORE_CUTOFF) {
        tmt6IonsPresent.push_back(tmt6PlexOnlyReporterIonMzs_[i]);
      }
    }
    if (!tmt6IonsPresent.empty()) {
      // declaring TMT6 present
      carp(CARP_DEBUG, "TMT 6plex ions present: %d. m/z values: %s",
           tmt6IonsPresent.size(), StringUtils::Join(tmt6IonsPresent, ',').c_str());
      significant.insert(TMT_6PLEX);
      tmt6 = true;
    }
  }

  // first check for TMT10, if checking is justified
  int nTmt10PeaksDetected = 0;
  carp(CARP_DEBUG, "Found MS3 scans? %d", foundMs3Scans_);
  bool tmt10 = false;
  if (tmt6 || foundMs3Scans_) {
    // Either we have TMT6/10, or there are MS3 scans. Either way, we might have TMT10.
    // Let's ask tmt6vs10detector
    RunAttributeResult tmtDetectorResult = tmt10Detector_.summarize();
    tmt10 = tmtDetectorResult.getValue("TMT10") == "T";
    nTmt10PeaksDetected = StringUtils::FromString<int>(
      tmtDetectorResult.getValue("n_peaks_with_enough_TMT10_signal"));
    if (tmt10) {
      carp(CARP_INFO, "TMT10 is present, %d TMT10 peaks detected.", nTmt10PeaksDetected);
      // declare TMT6 and TMT2 to be absent
      //tmt6 = false;
    } else {
      carp(CARP_INFO, "TMT10 is not present, %d TMT10 peaks detected.", nTmt10PeaksDetected);
    }
  }
  if (tmt6) {
    carp(CARP_INFO, "TMT: 6-plex reporter ions detected");
    result.addMod("K", SEARCH_MOD_MASS_TMT_6PLEX, true);
    result.addMod(Modification::LOCATION_NTERM, SEARCH_MOD_MASS_TMT_6PLEX, true);
    if (!tmt2) {
      carp(CARP_WARNING, "    No TMT 2-plex reporters detected, only 6-plex");
    }
    // special warning for TMT 6 vs. 2
    if (tmt6IonsPresent.size() < tmt6PlexOnlyReporterIonMzs_.size()) {
      carp(CARP_WARNING, "    TMT 6-plex detected, but not all ions present. Ions present:");
      for (vector<double>::const_iterator i = tmt6IonsPresent.begin(); i != tmt6IonsPresent.end(); i++) {
        carp(CARP_WARNING, "        %.4f", *i);
      }
    }
    tmt2 = tmt10 = false;
  } else if (tmt10) {
    result.addMod("K", SEARCH_MOD_MASS_TMT_6PLEX, true);
    result.addMod(Modification::LOCATION_NTERM, SEARCH_MOD_MASS_TMT_6PLEX, true);
    if (!tmt2) {
      carp(CARP_WARNING, "    No TMT 2-plex reporters detected, only 10-plex");
    }
    tmt2 = tmt6 = false;
  } else if (tmt2) {
    carp(CARP_INFO, "TMT: 2-plex reporter ions detected");
    result.addMod("K", SEARCH_MOD_MASS_TMT_2PLEX, true);
    result.addMod(Modification::LOCATION_NTERM, SEARCH_MOD_MASS_TMT_2PLEX, true);
    tmt6 = tmt10 = false;
  } else {
    carp(CARP_INFO, "TMT: no reporter ions detected");
  }

  // declare label to be present or not, and report the appropriate statistic.
  result.setValue("iTRAQ_8plex_present", itraq8 ? "T" : "F");
  result.setValue("iTRAQ_8plex_statistic",
                  StringUtils::ToString(reporterTypeTStatistics[ITRAQ_8PLEX]));
  result.setValue("iTRAQ_4plex_present", itraq4 ? "T" : "F");
  result.setValue("iTRAQ_4plex_statistic",
                  StringUtils::ToString(reporterTypeTStatistics[ITRAQ_4PLEX]));
  result.setValue("TMT_6plex_present", tmt6 ? "T" : "F");
  result.setValue("TMT_6plex_statistic",
                  StringUtils::ToString(reporterTypeTStatistics[TMT_6PLEX]));
  result.setValue("TMT_10plex_present", tmt10 ? "T" : "F");
  result.setValue("TMT_10plex_statistic",
                  StringUtils::ToString(nTmt10PeaksDetected));
  result.setValue("TMT_2plex_present", tmt2 ? "T" : "F");
  result.setValue("TMT_2plex_statistic",
                  StringUtils::ToString(reporterTypeTStatistics[TMT_2PLEX]));
  return result;
}

// Constants used by SilacDetector
const double SILAC_ZSCORE_CUTOFF = 4.0;
const int MAX_SCAN_SEPARATION = 50;

SilacDetector::SilacDetector() {
  modBinDistances_.push_back(4);
  modBinDistances_.push_back(6);
  modBinDistances_.push_back(8);
  modBinDistances_.push_back(10);
  controlBinDistances_.push_back(11);
  controlBinDistances_.push_back(14);
  controlBinDistances_.push_back(15);
  controlBinDistances_.push_back(21);
  controlBinDistances_.push_back(23);
  controlBinDistances_.push_back(27);
  modKExactMasses_[4] = 4.025107;
  modKExactMasses_[6] = 6.020129;
  modKExactMasses_[8] = 8.014199;
  modKExactMasses_[10] = 10.008269;
  modRExactMasses_[6] = 6.020129;
  modRExactMasses_[10] = 10.008269;
}

SilacDetector::~SilacDetector() {
}

void SilacDetector::processSpectrum(
  const Crux::Spectrum* spectrum,
  const vector<double>& binnedSpectrum
) {
  const vector<SpectrumZState>& zStates = spectrum->getZStates();
  if (zStates.size() != 1) {
    return;
  }
  scans_.push_back(spectrum->getFirstScan());
  double mass = zStates.front().getSinglyChargedMass();
  bins_.push_back(calcBinIndexMassPrecursor(mass));
}

RunAttributeResult SilacDetector::summarize() const {
  // construct the full set of all separations to summarize, both control and SILAC-label
  set<int> separations;
  for (vector<int>::const_iterator i = modBinDistances_.begin(); i != modBinDistances_.end(); i++) {
    separations.insert(*i);
  }
  for (vector<int>::const_iterator i = controlBinDistances_.begin();
       i != controlBinDistances_.end();
       i++) {
    separations.insert(*i);
  }
  // paranoia
  if (separations.size() < modBinDistances_.size() + controlBinDistances_.size()) {
    carp(CARP_WARNING, "A specified separation is also a control separation! Specified: %s",
         StringUtils::Join(modBinDistances_, ',').c_str());
  }

  // initialize a map from separation distances to counts of pairs with that separation
  map<int, int> countsWithSeparations;
  for (set<int>::const_iterator i = separations.begin(); i != separations.end(); i++) {
    countsWithSeparations[*i] = 0;
  }

  // keep track of the scan window defined by the minimum and maximum scan index to consider
  size_t minIdx = 0, maxIdx = 0;
  for (size_t i = 0; i < scans_.size(); i++) {
    // determine the minimum and maximum scan number currently in range
    int scan = scans_[i];
    int minScan = scan - MAX_SCAN_SEPARATION;
    int maxScan = scan + MAX_SCAN_SEPARATION;
    while (scans_[minIdx] < minScan) {
      minIdx++;
    }
    while (scans_[maxIdx] < maxScan && maxIdx < scans_.size() - 1) {
      maxIdx++;
    }
    // within the scan window, increment the separations that we care about with any that involve this scan
    for (size_t j = minIdx; j < maxIdx; j++) {
      int separation = abs(bins_[i] - bins_[j]);
      if (separations.find(separation) != separations.end()) {
        countsWithSeparations[separation]++;
      }
    }
  }

  // summarize the control separations
  double meanControlCount = 0.0;
  vector<int> controlCounts;
  carp(CARP_DEBUG, "SILAC: Control separation counts:");
  for (vector<int>::const_iterator i = controlBinDistances_.begin();
       i != controlBinDistances_.end();
       i++) {
    int count = countsWithSeparations[*i];
    meanControlCount += count;
    controlCounts.push_back(count);
    carp(CARP_DEBUG, "  %d: %d", *i, count);
  }
  meanControlCount /= controlBinDistances_.size();
  carp(CARP_DEBUG, "SILAC: Mean control separation count: %.5f", meanControlCount);
  if (meanControlCount <= 0) {
    carp(CARP_WARNING, "SILAC: No counts for any control separation pairs! Cannot estimate "
         "prevalance of SILAC separations.");
    // make a dummy result with no significant inferences
    RunAttributeResult result;
    for (vector<int>::const_iterator i = modBinDistances_.begin();
         i != modBinDistances_.end();
         i++) {
      result.setValue("SILAC_" + StringUtils::ToString(*i) + "Da_present",
                      RunAttributeResult::ERROR);
      result.setValue("SILAC_" + StringUtils::ToString(*i) + "Da_statistic",
                      RunAttributeResult::ERROR);
    }
    return result;
  }
  carp(CARP_DEBUG, "SILAC: Counts for each separation:");
  for (vector<int>::const_iterator i = modBinDistances_.begin();
       i != modBinDistances_.end();
       i++) {
    int count = countsWithSeparations[*i];
    carp(CARP_DEBUG, "  %d: %d (proportion=%.5f)", *i, count, count / meanControlCount);
  }
  double controlSd = MathUtil::StdDev(controlCounts);

  // determine any separations with a significantly elevate number of representatives
  RunAttributeResult result;
  set<int> significant;
  carp(CARP_DEBUG, "SILAC: Ratios of mass separations to control separations:");
  for (vector<int>::const_iterator i = modBinDistances_.begin(); i != modBinDistances_.end(); i++) {
    // z-score is checked against a cutoff to determine significance
    int count = countsWithSeparations[*i];
    double zScoreToControl = (count - meanControlCount) / controlSd;
    if (zScoreToControl > SILAC_ZSCORE_CUTOFF) {
      significant.insert(*i);
      carp(CARP_INFO, "SILAC: %dDa separation detected.", *i);
      // paranoia
      map<int, double>::const_iterator findK = modKExactMasses_.find(*i);
      map<int, double>::const_iterator findR = modRExactMasses_.find(*i);
      bool inK = findK != modKExactMasses_.end();
      bool inR = findR != modRExactMasses_.end();
      if (!inK && !inR) {
        throw runtime_error("Unknown SILAC separation " + StringUtils::ToString(*i));
      }
      // find out the exact appropriate mass for the search
      if (inK) {
        result.addMod("K", findK->second, true);
      }
      if (inR) {
        result.addMod("R", findR->second, true);
      }
      result.setValue("SILAC_" + StringUtils::ToString(*i) + "Da_present", "T");
    } else {
      result.setValue("SILAC_" + StringUtils::ToString(*i) + "Da_present", "F");
    }
    result.setValue("SILAC_" + StringUtils::ToString(*i) + "Da_statistic",
                    StringUtils::ToString(zScoreToControl));
    // show some details for debug output
    carp(CARP_DEBUG, "SILAC:     %dDa: %.5f (z=%.3f)",
         *i, count / meanControlCount, zScoreToControl);
  }
  if (significant.empty()) {
    carp(CARP_INFO, "SILAC: no labeling detected");
  }
  return result;
}

// Constants used by EnzymeDetector
const double MOD_OXIDATION_MASS_DIFF = 15.994915;
const double CTERM_MASS_UNMOD = 2 * MASS_H_AVERAGE + MOD_OXIDATION_MASS_DIFF;
const double PROPORTION_ONE_THIRD_PRECURSOR_Y1_THRESHOLD = 0.4;
const int MIN_SPECTRA_FOR_CONFIDENCE = 2000;
const double MIN_TRYPSIN_ZSCORE_THRESHOLD = 2.0;
const double MIN_TRYPSIN_ZSCORE_LIKELY_THRESHOLD = 1.5 ;
const double MIN_ARGC_LYSC_ZSCORE_THRESHOLD = 5.0;
const double MIN_PEPSIN_TSTAT_THRESHOLD = 5.0;
const double MIN_CHYMOTRYPSIN_NO_PEPSIN_TSTAT_THRESHOLD = 5.0;
const double MIN_BIN_PROPORTION_FOR_COUNT = 0.005;

EnzymeDetector::EnzymeDetector(const vector<Modification>& mods):
  numSpectra_(0), sumProportionMs2SignalOneThirdPrecursor_(0.0), cTermMass_(CTERM_MASS_UNMOD) {
  for (char c = 'A'; c <= 'Z'; c++) {
    double mass = AminoAcidUtil::GetMass(c, true);
    if (!std::isnan(mass)) {
      aaMasses_[c] = mass;
    }
  }

  set<string> modLocations;
  for (vector<Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
    carp(CARP_INFO, "Accounting for modification: %s", i->str().c_str());
    string location = i->getLocation();
    if (modLocations.find(location) != modLocations.end()) {
      carp(CARP_WARNING, "Multiple modifications on %s", i->getLocation().c_str());
    }
    if (location == Modification::LOCATION_CTERM) {
      // TODO: what if it's variable
      cTermMass_ += i->getMassDiff();
    } else if (location == Modification::LOCATION_NTERM) {
    } else if (location.length() == 1) {
      // TODO: what if it's variable
      aaMasses_[location[0]] += i->getMassDiff();
    }
  }

  for (map<char, double>::const_iterator i = aaMasses_.begin(); i != aaMasses_.end(); i++) {
    // y1 m/z is the c-terminus msas plus the AA mass plus a hydrogen
    const char aa = i->first;
    double mz = i->second + cTermMass_ + MASS_H_AVERAGE;
    aaY1Charge1Mz_[aa] = mz;
    aaY1Charge1Bin_[aa] = calcBinIndexMzFragment(mz);
    sumProportionsY1_[aa] = 0.0;
    sumProportionsBNMinus1_[aa] = 0.0;
  }
  carp(CARP_DEBUG, "Y1 ion mzs:");
  vector<string> y1Strings;
  for (map<char, double>::const_iterator i = aaY1Charge1Mz_.begin();
       i != aaY1Charge1Mz_.end();
       i++) {
    y1Strings.push_back(StringUtils::ToString(i->first) + ": " + StringUtils::ToString(i->second));
  }
  carp(CARP_DEBUG, "%s", StringUtils::Join(y1Strings, ',').c_str());

  controlAas_.insert('A');
  controlAas_.insert('C');
  controlAas_.insert('D');
  controlAas_.insert('E');
  controlAas_.insert('G');
  controlAas_.insert('H');
  controlAas_.insert('I');
  controlAas_.insert('K');
  controlAas_.insert('M');
  controlAas_.insert('N');
  controlAas_.insert('P');
  controlAas_.insert('R');
  controlAas_.insert('S');
  controlAas_.insert('V');
}

EnzymeDetector::~EnzymeDetector() {
}

void EnzymeDetector::processSpectrum(
  const Crux::Spectrum* spectrum,
  const vector<double>& binnedSpectrum
) {
  const vector<SpectrumZState>& zStates = spectrum->getZStates();
  if (zStates.size() != 1) {
    return;
  }
  numSpectra_++;

  // track how much signal is contained below 1/3 of the precursor m/z.
  // I'll use this to decide which ions to use to determine enzyme
  const double precursor = spectrum->getPrecursorMz();
  vector<Peak> peaks = spectrum->getPeaks();
  std::sort(peaks.begin(), peaks.end(), Peak::compareByMZ);
  double proportionOneThirdPrecursor = 0.0;
  for (vector<Peak>::const_iterator i = peaks.begin(); i != peaks.end(); i++) {
    proportionOneThirdPrecursor += i->getIntensity();
    if (i->getLocation() > precursor) {
      break;
    }
  }
  sumProportionMs2SignalOneThirdPrecursor_ +=
    proportionOneThirdPrecursor / spectrum->getTotalEnergy();

  int charge = zStates.front().getCharge();
  for (map<char, double>::const_iterator i = aaMasses_.begin(); i != aaMasses_.end(); i++) {
    int bNMinus1 = calcBinIndexMzFragment(
      (precursor - MASS_H_AVERAGE) * charge - i->second - cTermMass_ + MASS_H_AVERAGE);
    // if the appropriate bin for this AA is over the threshold, count it
    const char aa = i->first;
    double binY1 = binnedSpectrum[aaY1Charge1Bin_[aa]];
    if (binY1 > MIN_BIN_PROPORTION_FOR_COUNT) {
      sumProportionsY1_[aa]++;
    }
    if (bNMinus1 < binnedSpectrum.size()) {
      // if the appropriate bin for this AA is over the threshold, count it
      if (binnedSpectrum[bNMinus1] > MIN_BIN_PROPORTION_FOR_COUNT) {
        sumProportionsY1_[aa]++;
      }
    }
  }
}

RunAttributeResult EnzymeDetector::summarize() const {
  RunAttributeResult result;
  if (numSpectra_ < MIN_SPECTRA_FOR_CONFIDENCE) {
    carp(CARP_WARNING, "only %d spectra were analyzed. Enzyme determination suspect.",
         numSpectra_);
  }
  double proportionOneThirdPrecursor =
    sumProportionMs2SignalOneThirdPrecursor_ / numSpectra_;
  carp(CARP_DEBUG, "Signal proportion under 1/3 precursor mz: %f", proportionOneThirdPrecursor);
  const map<char, double>* aaProportionSumsPtr;
  if (proportionOneThirdPrecursor > PROPORTION_ONE_THIRD_PRECURSOR_Y1_THRESHOLD) {
    aaProportionSumsPtr = &sumProportionsY1_;
    carp(CARP_DEBUG, "Using y1 ions for enzyme determination.");
  } else {
    aaProportionSumsPtr = &sumProportionsBNMinus1_;
    carp(CARP_DEBUG, "Using b(n-1) ions for enzyme determination.");
  }
  const map<char, double>& aaProportionSums = *aaProportionSumsPtr;

  // Trypsin, ArgC and LysC.
  // Trypsin is just ArgC + LysC. So test R and K individually.
  // If the lower of the two is significant, Trypsin.
  // If not, then if R or K significant, then that one.
  set<char> controlAas = controlAas_;
  controlAas.erase('K');
  controlAas.erase('R');
  vector<double> controlValues;
  for (set<char>::const_iterator i = controlAas.begin(); i != controlAas.end(); i++) {
    map<char, double>::const_iterator j = aaProportionSums.find(*i);
    if (j != aaProportionSums.end()) {
      controlValues.push_back(j->second);
    }
  }
  double controlMean = MathUtil::Mean(controlValues);
  double controlSd = MathUtil::StdDev(controlValues);

  map<char, double> aaZScores;
  // debug report on each AA. Save up the messages. If enzyme = unknown, we'll print them
  vector<string> perAaMessages;
  for (map<char, double>::const_iterator i = aaProportionSums.begin();
       i != aaProportionSums.end();
       i++) {
    const char aa = i->first;
    const double x = i->second;
    const double zScore = (x - controlMean) / controlSd;
    aaZScores[aa] = zScore;
    char msg[64];
    sprintf(msg, "%c: spectra: %.4f. z-score = %.4f", aa, x / numSpectra_, zScore);
    carp(CARP_DEBUG, "%s", msg);
    perAaMessages.push_back(string(msg));
  }

  double zR = aaZScores['R'], zK = aaZScores['K'];
  carp(CARP_DEBUG, "  ArgC z-score: %f", zR);
  carp(CARP_DEBUG, "  LysC z-score: %f", zK);
  double trypsinMinZScore = min(zR, zK);
  carp(CARP_DEBUG, "  min(ArgC, LysC): %f", trypsinMinZScore);

  if (trypsinMinZScore > MIN_TRYPSIN_ZSCORE_THRESHOLD) {
    result.setValue("enzyme", "trypsin");
    return result;
  }

  // OK, we don't think it's trypsin. That's the main finding. But can we get more specific?
  if (zK > MIN_ARGC_LYSC_ZSCORE_THRESHOLD) {
    result.setValue("enzyme", "Lys-C");
    return result;
  } else if (zR > MIN_ARGC_LYSC_ZSCORE_THRESHOLD) {
    result.setValue("enzyme", "Arg-C");
    return result;
  }

  // Try Pepsin. Don't test against Pepsin or Chymotrypsin controls
  controlAas = controlAas_;
  controlAas.erase('F');
  controlAas.erase('L');
  controlAas.erase('W');
  controlAas.erase('Y');
  vector<double> testValues;
  map<char, double>::const_iterator lookup = aaProportionSums.find('F');
  if (lookup != aaProportionSums.end()) {
    testValues.push_back(lookup->second);
  }
  lookup = aaProportionSums.find('L');
  if (lookup != aaProportionSums.end()) {
    testValues.push_back(lookup->second);
  }
  controlValues.clear();
  for (set<char>::const_iterator i = controlAas.begin(); i != controlAas.end(); i++) {
    map<char, double>::const_iterator j = aaProportionSums.find(*i);
    if (j != aaProportionSums.end()) {
      controlValues.push_back(j->second);
    }
  }
  double pepsinTStatistic = tTestInd(testValues, controlValues);
  carp(CARP_DEBUG, "  Pepsin t-statistic: %f", pepsinTStatistic);
  // Try Chymotrypsin
  testValues.clear();
  lookup = aaProportionSums.find('W');
  if (lookup != aaProportionSums.end()) {
    testValues.push_back(lookup->second);
  }
  lookup = aaProportionSums.find('Y');
  if (lookup != aaProportionSums.end()) {
    testValues.push_back(lookup->second);
  }
  double chymotrypsinTStatistic = tTestInd(testValues, controlValues);
  carp(CARP_DEBUG, "  Chymotrypsin (but not Pepsin) t-statistic: %f", chymotrypsinTStatistic);

  // Test Pepsin. If passes, test separately for chymotrypsin. Return one or the other.
  if (pepsinTStatistic > MIN_PEPSIN_TSTAT_THRESHOLD) {
    result.setValue("enzyme", chymotrypsinTStatistic > MIN_CHYMOTRYPSIN_NO_PEPSIN_TSTAT_THRESHOLD ?
                    "Chymotrypsin" : "Pepsin");
    return result;
  }

  // Didn't detect anything else. If we can say it's "likely" trypsin, do so
  if (trypsinMinZScore > MIN_TRYPSIN_ZSCORE_LIKELY_THRESHOLD) {
    carp(CARP_INFO, "Trypsin test passes bare minimum threshold, and no other enzyme dtected.");
    result.setValue("enzyme", "Likely Trypsin");
    return result;
  }

  // We got nothin'. Give the user the information so they can decide.
  carp(CARP_INFO, "Unable to determine enzyme. Summary of each amino acid:");
  for (vector<string>::const_iterator i = perAaMessages.begin(); i != perAaMessages.end(); i++) {
    carp(CARP_INFO, "%s", i->c_str());
  }
  result.setValue("enzyme", "Unknown");
  return result;
}

}

