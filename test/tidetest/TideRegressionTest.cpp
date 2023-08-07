#include "TideRegressionTest.h"
#include "TideRegressionTestDriver.h"

#include <cstring>
#include <iostream>
#include <libgen.h>

using namespace std;

void initSettings() {
  addTest("standard", TideRegressionSettings(
    fasta_, "trypsin", "full-digest", 0, 6, 50, 200.0, 7200.0, false,
    "C+57.02146", spectrumRecords_, 3.0, 5));
  addTest("enzymatic digestion", TideRegressionSettings(
    fasta_, "chymotrypsin", "partial-digest", 2, 6, 50, 200.0, 7200.0, false,
    "C+57.02146", spectrumRecords_, 3.0, 5));
  addTest("peptide length", TideRegressionSettings(
    fasta_, "trypsin", "full-digest", 0, 10, 12, 200.0, 7200.0, false,
    "C+57.02146", spectrumRecords_, 3.0, 5));
  addTest("peptide mass", TideRegressionSettings(
    fasta_, "trypsin", "full-digest", 0, 6, 50, 1200.0, 1300.0, false,
    "C+57.02146", spectrumRecords_, 3.0, 5));
  addTest("monoisotopic precursor", TideRegressionSettings(
    fasta_, "trypsin", "full-digest", 0, 6, 50, 200.0, 7200.0, true,
    "C+57.02146", spectrumRecords_, 3.0, 5));
  addTest("variable mods", TideRegressionSettings(
    fasta_, "trypsin", "full-digest", 0, 6, 50, 200.0, 7200.0, false,
    "C+57.02146,2M+15.9949", spectrumRecords_, 3.0, 5));
}

void addTest(const string& name, const TideRegressionSettings& settings) {
  settings_[name] = settings;
}

int main(int argc, char** argv) {

  if (argc < 2 || argc > 3) {
    printUsage(argv[0]);
    return 1;
  } else if (argc == 2 && string(argv[1]) == "clean") {
    TideRegressionTestDriver::cleanTestFiles();
    return 0;
  }

  setPaths(argv[1]);
  initSettings();

  if (argc == 2) {
    for (map<string, TideRegressionSettings>::const_iterator i = settings_.begin();
         i != settings_.end();
         ++i) {
      runTest(i);
    }
  } else {
    map<string, TideRegressionSettings>::const_iterator i = settings_.find(argv[2]);
    if (i == settings_.end()) {
      cout << "Test \"" << argv[2] << "\" not found. Available tests are:" << endl;
      listTests();
      return 1;
    }
    runTest(i);
  }

  return 0;
}

void setPaths(const char* cruxPath) {
  char cruxDirChars[strlen(cruxPath) + 1];
  strcpy(cruxDirChars, cruxPath);
  string cruxDir = dirname(cruxDirChars);
  tideIndex_ = cruxDir + "/tide/tide-index";
  tideSearch_ = cruxDir + "/tide/tide-search";
  crux_ = cruxPath;
  //fasta_ = cruxDir + "/test/smoke-tests/small-yeast.fasta";
  //spectrumRecords_ = cruxDir + "/test/smoke-tests/demo.spectrumrecords";
  fasta_ = cruxDir + "/tidetest/worm.fasta";
  spectrumRecords_ = cruxDir + "/tidetest/worm-06-10000.spectrumrecords";
}

void runTest(map<string, TideRegressionSettings>::const_iterator i) {
  cout << "Running test \"" << i->first << "\"..." << flush;
  TideRegressionTestDriver driver(i->second);
  bool result = driver.runTest(tideIndex_, tideSearch_, crux_);
  cout << (result ? "passed" : "failed (" + driver.getError() + ")")
       << endl << flush;
}

void listTests() {
  for (map<string, TideRegressionSettings>::const_iterator i = settings_.begin();
       i != settings_.end();
       ++i) {
    cout << '\t' << i->first << endl;
  }
}

void printUsage(const char* argv0) {
  char arg[strlen(argv0) + 1];
  strcpy(arg, argv0);
  cout << "usage: " << basename(arg) << " "
          "(clean|<crux path> [<test name>])" << endl;
}

