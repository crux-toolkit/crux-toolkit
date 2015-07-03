/**
 * \file PinWriter.h
 * \brief Writes search results in the .pin format.
 */
#ifndef PINWRITER_H
#define PINWRITER_H

#include <set>
#include <string>
#include <vector>
#include "objects.h"
#include "model/Spectrum.h"
#include "util/mass.h"
#include "model/Match.h"
#include "model/MatchCollection.h"
#include "model/SpectrumZState.h"
#include "model/Spectrum.h"
#include "model/Peptide.h"
#include <limits>

class PinWriter {
 public:
  PinWriter();
  ~PinWriter();
  PinWriter(const char* output_file);
 
  void write(
    MatchCollection* target_collection,
    const std::vector<MatchCollection*>& decoys,
    int top_rank
  );
  
  void printHeader();

  void closeFile();
  void openFile(
    const std::string& filename, 
    const std::string& output_directory,
    bool overwrite
  );

  void setEnabledStatus(const std::string& name, bool enabled);

 protected:
  std::vector< std::pair<std::string, bool> > features_;
  std::vector<std::string> enabledFeatures_;
  std::ofstream* out_;
  ENZYME_T enzyme_; 
  int precision_;
  int mass_precision_;

  void printPSM(Crux::Match* match);

  std::string getPeptide(Crux::Peptide* peptide);
  bool isInfinite(FLOAT_T x);
  std::string getId(Crux::Match* match, int scan_number); 

  struct IsFeature : public std::unary_function<const std::pair<std::string, bool>&, bool> {
    IsFeature(const std::string& name): search_(name) {}
    bool operator() (const std::pair<std::string, bool>& check) {
      return search_ == check.first;
    }
    std::string search_;
  };

  struct FeatureCopy {
    FeatureCopy(std::vector<std::string>* target): target_(target) {}
    void operator() (const std::pair<std::string, bool>& feature) {
      if (feature.second) target_->push_back(feature.first);
    }
    std::vector<std::string>* target_;
  };
};

#endif // PINWRITER_H

