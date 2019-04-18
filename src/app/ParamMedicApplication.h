#ifndef PARAM_MEDIC_APPLICATION_H
#define PARAM_MEDIC_APPLICATION_H

#include "CruxApplication.h"
#include "model/Spectrum.h"

namespace ParamMedic {
  class RunAttributeResult;
  class Modification;
}

class ParamMedicApplication : public CruxApplication {
 public:
  ParamMedicApplication();
  virtual ~ParamMedicApplication();
  virtual int main(int argc, char** argv);
  static void processFiles(
    const std::vector<std::string>& files,
    bool errorCalcEnabled,
    bool modsEnabled,
    ParamMedic::RunAttributeResult* errorCalcResult,
    std::vector<ParamMedic::RunAttributeResult>* modsResult
  );
  virtual std::string getName() const;
  virtual std::string getDescription() const;
  virtual bool hidden() const;
  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;
  virtual bool needsOutputDirectory() const;
};

namespace ParamMedic {

class Modification {
 public:
  Modification();
  Modification(const std::string& location, double massDiff, bool isVariable);
  Modification(const Modification& other);
  virtual ~Modification();

  friend void swap(Modification& x, Modification& y);
  Modification& operator=(Modification rhs);

  std::string getLocation() const;
  double getMassDiff() const;
  bool getVariable() const;
  std::string str() const;

  static std::vector<Modification> GetFromResults(const std::vector<RunAttributeResult>& results);

  static const std::string LOCATION_NTERM;
  static const std::string LOCATION_CTERM;
 private:
  std::string location_;
  double massDiff_;
  bool variable_;
};

class RunAttributeResult {
 public:
  RunAttributeResult();
  RunAttributeResult(const RunAttributeResult& other);
  virtual ~RunAttributeResult();

  friend void swap(RunAttributeResult& x, RunAttributeResult& y);
  RunAttributeResult& operator=(RunAttributeResult rhs);

  std::string getValue(const std::string& name) const;
  void setValue(const std::string& name, const std::string& value);
  void addMod(const std::string& location, double massDiff, bool isVariable);
  const std::vector<Modification>& getMods() const;

  static std::string getFirstValue(const std::vector<RunAttributeResult>& results, const std::string& name);

  static const std::string ERROR;
 private:
  // modifications recommended as a result of analysis
  std::vector<Modification> modifications_;
  // name-value pairs to output that summarize the analysis
  std::map<std::string, std::string> nameValuePairs_;
};

class RunAttributeDetector {
 public:
  virtual void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum) = 0;
  virtual void nextFile() {}
  virtual RunAttributeResult summarize() const { return RunAttributeResult(); }
};

class PerChargeErrorCalc;

class ErrorCalc : public RunAttributeDetector {
 public:
  ErrorCalc();
  ~ErrorCalc();

  void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum);

  void nextFile();

  void calcMassErrorDist(
    std::string* precursorFailure,
    std::string* fragmentFailure,
    double* precursorSigmaPpm,
    double* fragmentSigmaPpm,
    double* precursorPredictionPpm,
    double* fragmentPredictionTh
  ) const;

  RunAttributeResult summarize() const;

  static const std::string KEY_MESSAGES;
  static const std::string KEY_PRECURSOR_FAILURE;
  static const std::string KEY_FRAGMENT_FAILURE;
  static const std::string KEY_PRECURSOR_SIGMA;
  static const std::string KEY_FRAGMENT_SIGMA;
  static const std::string KEY_PRECURSOR_PREDICTION;
  static const std::string KEY_FRAGMENT_PREDICTION;

 private:
  std::map<int, PerChargeErrorCalc*> calcs_;
  int numTotalSpectra_;
};

class PerChargeErrorCalc : public RunAttributeDetector {
 public:
  PerChargeErrorCalc(int charge);
  ~PerChargeErrorCalc();

  void processSpectrum(
    const Crux::Spectrum* spectrum,
    const std::vector<double>& binnedSpectrum);
  void clearBins();
  void nextFile();

  int getNumPassingSpectra() const;
  int getNumSpectraSameBin() const;
  int getNumSpectraWithinPpm() const;
  int getNumSpectraWithinPpmAndScans() const;
  int getNumMultipleFragBins() const;
  int getNumSingleFragBins() const;
  const std::vector< std::pair<Peak, Peak> >& getPairedFragmentPeaks() const;
  const std::vector< std::pair<double, double> >& getPairedPrecursorMzs() const;

 protected:
  int getBinIndexPrecursor(double mz);
  int getBinIndexFragment(double mz);
  double getPrecursorMz(const Crux::Spectrum* spectrum) const;

  // given two spectra, pair up their fragments that are in the same bin
  std::vector< std::pair<Peak, Peak> > pairFragments(
    const std::vector<Peak>& prev,
    const std::vector<Peak>& cur
  );

  // keep only one fragment per bin; if another fragment wants to be in the bin,
  // toss them both out - this reduces ambiguity
  std::map<int, Peak> binFragments(const std::vector<Peak>& peaks);

  static bool sortPairedFragments(
    const std::pair<Peak, Peak>& x,
    const std::pair<Peak, Peak>& y
  );

  int charge_;
  int chargeForBinSize_;
  // count the spectra that go by
  int numTotalSpectra_;
  int numPassingSpectra_;
  int numSpectraSameBin_;
  int numSpectraWithinPpm_;
  int numSpectraWithinPpmAndScans_;
  // number and position of bins
  int numMultipleFragBins_;
  int numSingleFragBins_;
  // map from bin index to current spectrum
  std::map< int, std::pair< const Crux::Spectrum*, std::vector<Peak> > > spectra_;
  // the paired peak values that we'll use to estimate mass error
  std::vector< std::pair<Peak, Peak> > pairedFragmentPeaks_;
  std::vector< std::pair<double, double> > pairedPrecursorMzs_;
};

class Model {
 public:
  Model(double nMean, double nStd, double nMinStd, double uStart, double uEnd);
  virtual ~Model();

  // fit the model to new data using EM.
  // this method fits the components of the model to new data using the EM method.
  // it will iterate until either max iterations has been reached, or the stop
  // threshold has been passed.
  double fit(const std::vector<double>& data);

  double getMu() const;
  double getSigma() const;

  // summarize a batch of data and store sufficient statistics.
  // this will run the expectation step of EM and store sufficient statistics in
  // the appropriate distribution objects.
  // the summarization can be thought of as a chunk of the E step, and the
  // fromSummaries method as the M step.
  double summarize(const std::vector<double>& x);

  // fit the model to the collected sufficient statistics.
  // fit the parameters of the model to the sufficient statistics gathered during
  // the summarize calls.
  // this should return an exact update.
  void fromSummaries();

  // clear the summary statistics stored in the object
  void clearSummaries();

  // perform log-sum-exp on a pair of numbers in a log space.
  // this is calculated as z = log( e^x + e^y ).
  // however, this causes underflow sometimes when x or y are too negative.
  // a simplification of this is thus z = x + log( e^(y-x) + 1 ), where x is the
  // greater number.
  // if either of the inputs are infinity, return infinity, and if either of the
  // inputs are negative infinity, then simply return the other input.
  static double pairLse(double x, double y);
 protected:
  class NormalDistribution {
   public:
    NormalDistribution(double mean, double std, double minStd);
    virtual ~NormalDistribution();
    double getMu() const;
    double getSigma() const;
    void logProbability(const std::vector<double>& x, std::vector<double>* r) const;
    void summarize(const std::vector<double>& x, double* r);
    void fromSummaries();
    void clearSummaries();
   protected:
    double mu_;
    double sigma_;
    double minStd_;
    double logSigmaSqrt2Pi_;
    double twoSigmaSquared_;
    double summaries_[3];
  };

  class UniformDistribution {
   public:
    UniformDistribution(double start, double end);
    virtual ~UniformDistribution();
    void logProbability(const std::vector<double>& x, std::vector<double>* r) const;
    void summarize(const std::vector<double>& x, double* r);
    void summarize();
    void fromSummaries();
    void clearSummaries();
   protected:
    double start_;
    double end_;
    double logP_;
    double summaries_[2];
  };

  NormalDistribution normal_;
  UniformDistribution uniform_;
  double weights_[2];
  double summaries_[2];
};

// estimate mu and sigma of the mixed distribution, as initial estimate for Gaussian
static void estimateMuSigma(
  const std::vector<double>& data,
  double minSigma,
  double* muFit,
  double* sigmaFit
);

int calcBinIndexMassPrecursor(double mass);
int calcBinIndexMzFragment(double mz);
double calcMH(double mz, int charge);
std::vector<double> binSpectrum(const Crux::Spectrum* spectrum);
int processSpectra(
  const std::vector<std::string>& files,
  std::vector<RunAttributeDetector*> detectors);
int processSpectra(
  const std::vector<std::string>& files,
  RunAttributeDetector* detector);
}

#endif

