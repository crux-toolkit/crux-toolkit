// Tahmina Baker
//
// This file contains classes for reporting search results. 
//
// The Reporter class is a generic base class for displaying peptide-spectrum
// matches found during a search. It is intended to be used as follows:
//    1. Call ReportSpectrum for each unique spectrum being searched
//    2. Call ReportMatch for each peptide that matches the current spectrum
//    3. Call WriteReport when all the matches for the current spectrum have
//       been reported.
//
// The TextReporter is a very basic Reporter that outputs results to stdout.
//
// The PBReporter class outputs results to the Results protocol buffer (see 
// Results.proto), which can later be parsed by a postprocessor.


#ifndef REPORT_H
#define REPORT_H

#include <iostream>
#include <vector>
#include "raw_proteins.pb.h"
#include "spectrum_collection.h"
#include "results.pb.h"
#include "header.pb.h"
#include "peptide.h"
#include "abspath.h"
#include "records.h"

using namespace std;
class Spectrum;

#define CHECK(x) GOOGLE_CHECK((x))

// Base class to display peptide-spectrum matches.
class Reporter {
 public:
  virtual void ReportSpectrum(const Spectrum* spectrum, int charge, 
                              int spectrum_index, pb::Stats* stats) = 0;
  virtual void ReportMatch(int score, const Peptide& peptide) = 0;
  virtual void WriteReport() = 0;
};

// Very basic output
class TextReporter : public Reporter {
 public:
  virtual void ReportSpectrum(const Spectrum* spectrum, int charge,
                              int spectrum_index, pb::Stats* stats);
  virtual void ReportMatch(int score, const Peptide& peptide); 
  virtual void WriteReport();

 private:
  const Spectrum *spectrum_;
  int charge_;
};

// Reporter that writes the results to a Results protocol buffer file that 
// can later be used by a postprocessor.
class PBReporter : public Reporter {
 public:
  PBReporter(const string& filename, pb::Header &results_header);

  // Reporter methods
  virtual void ReportSpectrum(const Spectrum* spectrum, int charge,
                              int spectrum_index, pb::Stats* stats);
  virtual void ReportMatch(int score, const Peptide& peptide); 
  virtual void WriteReport();
   
 private:
  const Spectrum *spectrum_;
  HeadedRecordWriter results_writer_;
  pb::Results pb_results_;
};

#endif // REPORT_H
