#ifndef PARAMMEDIC_H
#define PARAMMEDIC_H

#include "CruxApplication.h"
#include "Spectrum.h"

class ParamMedicApplication : public CruxApplication {
 public:
  ParamMedicApplication();
  virtual ~ParamMedicApplication();
  virtual int main(int argc, char** argv);
  virtual std::string getName() const;
  virtual std::string getDescription() const;
  virtual bool hidden() const;
  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;
  virtual bool needsOutputDirectory() const;
};

class ParamMedicErrorCalculator {
 public:
  ParamMedicErrorCalculator();
  virtual ~ParamMedicErrorCalculator();

  void processFiles(const std::vector<std::string>& files);
  void processSpectrum(Crux::Spectrum* spectrum);
  void clearBins();

  // this is to be run after all spectra have been processed;
  // fits the mixed model to the mixed distributions of m/z differences
  void calcMassErrorDist(
    std::string* precursorFailure,
    std::string* fragmentFailure,
    double* precursorSigmaPpm,
    double* fragmentSigmaPpm,
    double* precursorPredictionPpm,
    double* fragmentPredictionTh
  );

  // estimate mu and sigma of the mixed distribution, as initial estimate for Gaussian
  static void estimateMuSigma(
    const std::vector<double>& data,
    double minSigma,
    double* muFit,
    double* sigmaFit
  );
 protected:
  int getBinIndexPrecursor(double mz) const;
  int getBinIndexFragment(double mz) const;
  double getPrecursorMz(const Crux::Spectrum* spectrum) const;

  // given two spectra, pair up their fragments that are in the same bin
  std::vector< std::pair<const Peak*, const Peak*> > pairFragments(
    const Crux::Spectrum* prev,
    const Crux::Spectrum* cur
  );

  // keep only one fragment per bin; if another fragment wants to be in the bin,
  // toss them both out - this reduces ambiguity
  std::map<int, const Peak*> binFragments(const Crux::Spectrum* spectrum);

  static bool sortPairedFragments(
    const std::pair<const Peak*, const Peak*> x,
    const std::pair<const Peak*, const Peak*> y
  );

  // count the spectra that go by
  int numTotalSpectra_;
  int numPassingSpectra_;
  int numSpectraSameBin_;
  int numSpectraWithinPpm_;
  int numSpectraWithinPpmAndScans_;
  // number and position of bins
  double lowestPrecursorBinStartMz_;
  double lowestFragmentBinStartMz_;
  int numPrecursorBins_;
  int numFragmentBins_;
  int numMultipleFragBins_;
  int numSingleFragBins_;
  // map from bin index to current spectrum
  std::map<int, Crux::Spectrum*> spectra_;
  // the paired peak values that we'll use to estimate mass error
  std::vector< std::pair<const Peak*, const Peak*> > pairedFragmentPeaks_;
  std::vector< std::pair<double, double> > pairedPrecursorMzs_;
};

class ParamMedicModel {
 public:
  ParamMedicModel(double nMean, double nStd, double nMinStd, double uStart, double uEnd);
  virtual ~ParamMedicModel();

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

#endif

