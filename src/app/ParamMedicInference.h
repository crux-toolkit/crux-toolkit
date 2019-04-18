#ifndef PARAM_MEDIC_INFERENCE_H
#define PARAM_MEDIC_INFERENCE_H

#include "ParamMedicApplication.h"
#include "util/MathUtil.h"
#include <map>

namespace ParamMedic {

class PhosphoLossProportionCalc : public RunAttributeDetector {
 public:
  PhosphoLossProportionCalc();
  ~PhosphoLossProportionCalc();
  void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum);
  RunAttributeResult summarize() const;
 private:
  double sumProportionsInPhosphoLoss_;
  int numSpectraUsed_;
  std::map<int, double> sumsProportionsPerControlPeak_;
};

class Tmt6vs10Detector : public RunAttributeDetector {
 public:
  Tmt6vs10Detector();
  ~Tmt6vs10Detector();
  void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum);
  RunAttributeResult summarize() const;
 private:
  std::map< int, std::vector<double> > nominalMassAllPeaks_;
  std::map<int, double> nominalMassTmt6Mass_; // map from nominal mass to precise TMT6 mass
  std::map<int, double> nominalMassTmt10Mass_; // map from nominal mass to precise TMT10 mass
  std::map<int, double> nominalMassMinBinMass_;
};

class ReporterIonProportionCalc : public RunAttributeDetector {
 public:
  ReporterIonProportionCalc();
  ~ReporterIonProportionCalc();
  void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum);
  RunAttributeResult summarize() const;
 private:
  enum REPORTER_ION_TYPE { TMT_2PLEX, TMT_6PLEX, ITRAQ_4PLEX, ITRAQ_8PLEX, CONTROL };
  std::vector<double> tmt6PlexOnlyReporterIonMzs_;
  std::map< REPORTER_ION_TYPE, vector<double> > reporterIonTypeBins_;
  std::map< REPORTER_ION_TYPE, std::map<int, double> > reporterTypeBinSumProportion_;
  std::map<REPORTER_ION_TYPE, double> reporterIonTStatThresholds_;
  Tmt6vs10Detector tmt10Detector_;
  bool foundMs3Scans_;
};

class SilacDetector : public RunAttributeDetector {
 public:
  SilacDetector();
  ~SilacDetector();
  void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum);
  RunAttributeResult summarize() const;
 private:
  std::vector<int> modBinDistances_;
  std::vector<int> controlBinDistances_;
  std::map<int, double> modKExactMasses_, modRExactMasses_;
  std::vector<int> scans_;
  std::vector<int> bins_;
};

class EnzymeDetector : public RunAttributeDetector {
 public:
  EnzymeDetector(const vector<Modification>& mods);
  ~EnzymeDetector();
  void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum);
  RunAttributeResult summarize() const;
 private:
  double tTest() const;

  int numSpectra_;
  double sumProportionMs2SignalOneThirdPrecursor_;
  double cTermMass_;
  std::map<char, double> sumProportionsY1_;
  std::map<char, double> sumProportionsBNMinus1_;
  std::map<char, double> aaY1Charge1Mz_;
  std::map<char, int> aaY1Charge1Bin_;
  std::map<char, double> aaMasses_;
  std::set<char> controlAas_;
};

template<typename T>
static double tTestInd(const T& x, const T& y, bool equalVar = true) {
  return equalVar
    ? (MathUtil::Mean(x) - MathUtil::Mean(y)) / (
      sqrt(((x.size() - 1) * MathUtil::Variance(x, false) + (y.size() - 1) * MathUtil::Variance(y, false)) /
        (x.size() + y.size() - 2)) *
      sqrt((double)1 / x.size() + (double)1 / y.size()))

    : (MathUtil::Mean(x) - MathUtil::Mean(y)) /
      sqrt((MathUtil::Variance(x, false) / x.size()) + (MathUtil::Variance(y, false) / y.size()));
}

}

#endif

