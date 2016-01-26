#ifndef XHHC_SCORE_PEPTIDE_SPECTRUM_H
#define XHHC_SCORE_PEPTIDE_SPECTRUM_H

#include "app/CruxApplication.h"
#include "model/Spectrum.h"
#include "LinkedIonSeries.h"

class XLinkScoreSpectrum : public CruxApplication {
 public:
  XLinkScoreSpectrum();
  ~XLinkScoreSpectrum();

  virtual int main(int argc, char** argv);
  virtual void processParams();
  virtual std::string getName() const;
  virtual std::string getDescription() const;
  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

 private:
  double get_concat_score(char* peptideA, char* peptideB, int link_site,
                          int charge, Crux::Spectrum* spectrum);
  FLOAT_T* get_observed_raw(Crux::Spectrum* spectrum, int charge);
  void print_spectrum(Crux::Spectrum* spectrum, LinkedIonSeries& ion_series);
};

#endif

