#include "TideRegressionTestDriver.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>

using namespace std;

const string TideRegressionTestDriver::TIDEPEPIX = "tidereg-pepix";
const string TideRegressionTestDriver::TIDEPROTIX = "tidereg-protix";
const string TideRegressionTestDriver::TIDEAUXLOCS = "tidereg-auxlocs";
const string TideRegressionTestDriver::TIDEOUT = "tidereg-tidesearch";
const string TideRegressionTestDriver::CRUXIDX = "tidereg-cruxidx";
const string TideRegressionTestDriver::CRUXOUT = "tidereg-cruxsearch";

TideRegressionTestDriver::TideRegressionTestDriver(
  const TideRegressionSettings& settings):
  settings_(settings) {
}

TideRegressionTestDriver::TideRegressionTestDriver() {
}

TideRegressionTestDriver::~TideRegressionTestDriver() {
}

bool TideRegressionTestDriver::runTest(
  const string& tideIndexPath,
  const string& tideSearchPath,
  const string& cruxPath
) {
  error_.clear();

  if (!isRegularFile(tideIndexPath)) {
    error_ = "tide-index path \"" + tideIndexPath + "\" is invalid";
    return false;
  } else if (!isRegularFile(tideSearchPath)) {
    error_ = "tide-index path \"" + tideSearchPath + "\" is invalid";
    return false;
  } else if (!isRegularFile(cruxPath)) {
    error_ = "crux path \"" + cruxPath + "\" is invalid";
    return false;
  }

  stringstream cmdBuilder;
  string cmd;

  // clean up files from previous testing
  cmd = "rm -f " + TIDEPEPIX + "* " + TIDEPROTIX + " " + TIDEAUXLOCS;
  system(cmd.c_str());

  // tide-index
  cmdBuilder << tideIndexPath << " "
             "--fasta=" << settings_.fasta << " "
             "--enzyme=" << settings_.enzyme << " "
             "--digestion=" << settings_.digestion << " "
             "--max_missed_cleavages=" << settings_.missedCleavages << " "
             "--max_length=" << settings_.maxLength << " "
             "--max_mass=" << settings_.maxMass << " "
             "--min_length=" << settings_.minLength << " "
             "--min_mass=" << settings_.minMass << " "
             "--monoisotopic_precursor=" << (settings_.monoisotopicPrecursor ?
                                             'T' : 'F') << " "
             "--mods_spec=" << settings_.modsSpec << " "
             "--peptides=" << TIDEPEPIX << " "
             "--proteins=" << TIDEPROTIX << " "
             "--aux_locations=" << TIDEAUXLOCS << " "
             "&>/dev/null";
  cmd = cmdBuilder.str();
  system(cmd.c_str());

  // check that we got the expected output files
  if (!isRegularFile(TIDEPEPIX) || !isRegularFile(TIDEPROTIX) ||
      !isRegularFile(TIDEAUXLOCS)) {
    error_ = "tide-index failed to run";
    return false;
  }

  // clean up files from previous testing
  cmd = "rm -f " + TIDEOUT;
  system(cmd.c_str());

  // tide-search
  cmdBuilder.str("");
  cmdBuilder << tideSearchPath << " "
             "--mass_window=" << settings_.massWindow << " "
             "--top_matches=" << settings_.topMatch << " "
             "--spectra=" << settings_.spectrumRecords << " "
             "--peptides=" << TIDEPEPIX << " "
             "--proteins=" << TIDEPROTIX << " "
             "> " << TIDEOUT;
  cmd = cmdBuilder.str();
  system(cmd.c_str());

  // check that we got the expected output files
  if (!isRegularFile(TIDEOUT)) {
    error_ = "tide-search failed to run";
    return false;
  }

  // clean up files from previous testing
  cmd = "rm -rf " + CRUXIDX;
  system(cmd.c_str());

  // crux tide-index
  cmdBuilder.str("");
  cmdBuilder << cruxPath << " tide-index "
             "--enzyme " << settings_.enzyme << " "
             "--digestion " << settings_.digestion << " "
             "--missed-cleavages " << settings_.missedCleavages << " "
             "--max-length " << settings_.maxLength << " "
             "--min-length " << settings_.minLength << " "
             "--max-mass " << settings_.maxMass << " "
             "--min-mass " << settings_.minMass << " "
             "--monoisotopic-precursor " << (settings_.monoisotopicPrecursor ?
                                             'T' : 'F') << " "
             "--use-flanking-peaks T "
             "--mods-spec " << settings_.modsSpec << " "
             "--decoy-format none "
             "--overwrite T "
             "--output-dir " << CRUXOUT << " " <<
             settings_.fasta << ' ' << CRUXIDX << " "
             "&>/dev/null";
  cmd = cmdBuilder.str();
  system(cmd.c_str());

  // check that we got the expected output files
  if (!isDir(CRUXIDX)) {
    error_ = "crux tide-index failed to run";
    return false;
  }

  // clean up files from previous testing
  cmd = "rm -rf " + CRUXOUT;
  system(cmd.c_str());

  // crux tide-search
  cmdBuilder.str("");
  cmdBuilder << cruxPath << " tide-search "
             "--precursor-window " << settings_.massWindow << " "
             "--precursor-window-type mass "
             "--top-match " << settings_.topMatch << " "
             "--overwrite T "
             "--output-dir " << CRUXOUT << " " <<
             settings_.spectrumRecords << ' ' << CRUXIDX << " "
             "&>/dev/null";
  cmd = cmdBuilder.str();
  system(cmd.c_str());

  // check that we got the expected output files
  if (!isDir(CRUXOUT)) {
    error_ = "crux tide-search failed to run";
    return false;
  }

  return compareResults();
}

void TideRegressionTestDriver::cleanTestFiles() {
  string cmd;
  // clean tide-index files
  cmd = "rm -f " + TIDEPEPIX + "* " + TIDEPROTIX + " " + TIDEAUXLOCS;
  system(cmd.c_str());
  // clean tide-search file
  cmd = "rm -f " + TIDEOUT;
  system(cmd.c_str());
  // clean crux tide-index files
  cmd = "rm -rf " + CRUXIDX;
  system(cmd.c_str());
  // clean crux tide-search file
  cmd = "rm -rf " + CRUXOUT;
  system(cmd.c_str());
}

bool TideRegressionTestDriver::isRegularFile(const string& path) {
  struct stat sb;
  if (stat(path.c_str(), &sb) == -1) {
    return false;
  }
  return S_ISREG(sb.st_mode);
}

bool TideRegressionTestDriver::isDir(const string& path) {
  struct stat sb;
  if (stat(path.c_str(), &sb) == -1) {
    return false;
  }
  return S_ISDIR(sb.st_mode);
}

string TideRegressionTestDriver::getError() {
  return error_;
}

void TideRegressionTestDriver::loadSettings(
  const TideRegressionSettings& settings
) {
  settings_ = settings;
}

void TideRegressionTestDriver::setFasta(const string& fasta) {
  settings_.fasta = fasta;
}

void TideRegressionTestDriver::setEnzyme(const string& enzyme) {
  settings_.enzyme = enzyme;
}

void TideRegressionTestDriver::setDigestion(const string& digestion) {
  settings_.digestion = digestion;
}

void TideRegressionTestDriver::setMissedCleavages(int missedCleavages) {
  settings_.missedCleavages = missedCleavages;
}

void TideRegressionTestDriver::setMinLength(int minLength) {
  settings_.minLength = minLength;
}

void TideRegressionTestDriver::setMaxLength(int maxLength) {
  settings_.maxLength = maxLength;
}

void TideRegressionTestDriver::setMinMass(double minMass) {
  settings_.minMass = minMass;
}

void TideRegressionTestDriver::setMaxMass(double maxMass) {
  settings_.maxMass = maxMass;
}

void TideRegressionTestDriver::setMonoisotopicPrecursor(bool monoisotopicPrecursor) {
  settings_.monoisotopicPrecursor = monoisotopicPrecursor;
}

void TideRegressionTestDriver::setModsSpec(const string& modsSpec) {
  settings_.modsSpec = modsSpec;
}

void TideRegressionTestDriver::setSpectrumRecords(const string& spectrumRecords) {
  settings_.spectrumRecords = spectrumRecords;
}

void TideRegressionTestDriver::setMassWindow(double massWindow) {
  settings_.massWindow = massWindow;
}

void TideRegressionTestDriver::setTopMatch(int topMatch) {
  settings_.topMatch = topMatch;
}

vector<string> TideRegressionTestDriver::splitTab(const string& row) {
  vector<string> tokens;

  string::const_iterator from = row.begin();
  for (string::const_iterator i = from; i != row.end(); ++i) {
    if (*i == '\t') {
      tokens.push_back(string(from, i));
      from = i + 1;
    }
  }
  tokens.push_back(string(from, row.end()));
  return tokens;
}

bool TideRegressionTestDriver::compareResults() {
  ifstream f1(string(CRUXOUT + "/tide-search.target.txt").c_str());
  ifstream f2(TIDEOUT.c_str());

  if (!f1.good() || !f2.good()) {
    error_ = "couldn't open search results file";
    return false;
  }

  string s1, s2;
  vector<string> v1, v2;
  map<string, int> colsCrux, colsTide;

  // read header line from crux output and set up crux column map
  getline(f1, s1);
  v1 = splitTab(s1);
  for (size_t i = 0; i < v1.size(); ++i) {
    colsCrux[v1[i]] = i;
  }

  // set up tide column map
  colsTide["scan"] = 0;
  colsTide["spectrum precursor m/z"] = 1;
  colsTide["charge"] = 2;
  colsTide["xcorr score"] = 3;
  colsTide["sequence"] = 4;

  // set up comparison map
  map<string, COMPARE_TYPE> comparisons;
  comparisons["scan"] = EXACT;
  comparisons["spectrum precursor m/z"] = NUMERIC;
  comparisons["charge"] = EXACT;
  comparisons["xcorr score"] = NUMERIC;
  comparisons["sequence"] = PEPTIDE;

  getline(f1, s1);
  getline(f2, s2);
  size_t psmNum = 1;
  while (f1.good() && f2.good()) {
    v1 = splitTab(s1);
    v2 = splitTab(s2);
    if (v1.size() != colsCrux.size() || v2.size() != colsTide.size()) {
      // column header count did not match number of values
      error_ = "column count mismatch";
      return false;
    }
    // compare columns
    for (map<string, COMPARE_TYPE>::const_iterator i = comparisons.begin();
         i != comparisons.end();
         ++i) {
      map<string, int>::const_iterator cruxLookup = colsCrux.find(i->first);
      if (cruxLookup == colsCrux.end()) {
        // column not found in crux
        error_ = "crux results file missing column \"" + i->first + "\"";
        return false;
      }
      map<string, int>::const_iterator tideLookup = colsTide.find(i->first);
      if (tideLookup == colsCrux.end()) {
        // column not found in tide
        error_ = "tide results file missing column \"" + i->first + "\"";
        return false;
      }
      if (!valueCompare(v1[cruxLookup->second], v2[tideLookup->second], i->second)) {
        // values did not match
        stringstream errorBuilder;
        errorBuilder << "on psm " << psmNum << ", "
                     "\"" + v1[cruxLookup->second] + "\" != "
                     "\"" + v2[tideLookup->second] + "\"";
        error_ = errorBuilder.str();
        return false;
      }
    }
    getline(f1, s1);
    getline(f2, s2);
    ++psmNum;
  }

  if (f1.good() || f2.good()) {
    // one of the files still had more rows
    error_ = "psm count differs";
    return false;
  }

  return true;
}

bool TideRegressionTestDriver::valueCompare(
  const string& v1,
  const string& v2,
  COMPARE_TYPE type
) {
  switch (type) {
  case EXACT:
    return v1 == v2;
  case NUMERIC:
    return numCompare(v1, v2);
  case PEPTIDE:
    return peptideCompare(v1, v2);
  default:
    return false;
  }
}

bool TideRegressionTestDriver::peptideCompare(const string& v1, const string& v2) {
  for (string::const_iterator i = v1.begin(), j = v2.begin();
       i != v1.end() && j != v2.end();
       ++i, ++j) {
    if (*i == '[' && *j == '[') {
      // found mod, compare the value between brackets
      string::const_iterator k = ++i;
      while (*(++k) != ']');
      string mod1(i, k);
      i = k;
      k = ++j;
      while (*(++k) != ']');
      string mod2(j, k);
      j = k;
      if (!numCompare(mod1, mod2)) {
        // delta did not match
        return false;
      }
    } else if (*i != *j) {
      return false;
    }
  }

  return true;
}

bool TideRegressionTestDriver::numCompare(const string& v1, const string& v2) {
  if (v1 == v2) {
    return true;
  }

  // didn't compare equal, so check if decimal places are different
  size_t decimal = v1.find('.', 0);
  size_t precision1 = (decimal == string::npos) ? 0 : v1.length() - decimal - 1;
  decimal = v2.find('.', 0);
  size_t precision2 = (decimal == string::npos) ? 0 : v2.length() - decimal - 1;

  size_t minPrecision = min(precision1, precision2);

  string newV1 = makePrecise(v1, minPrecision);
  string newV2 = makePrecise(v2, minPrecision);

  if (newV1 == newV2) {
    return true;
  } else if (minPrecision > 0) {
    // check if +1 or -1 on most precise digit makes them same
    double tiny = pow(0.1, minPrecision);
    double m1 = atof(newV1.c_str());
    double m2 = atof(newV2.c_str());
    if (preciseCompare(m1 + tiny, m2, minPrecision) ||
        preciseCompare(m1 - tiny, m2, minPrecision)) {
      return true;
    }
  }

  return false;
}

string TideRegressionTestDriver::makePrecise(const string& v, size_t decimalPlaces) {
  stringstream ss;
  ss << fixed << setprecision(decimalPlaces) << atof(v.c_str());
  return ss.str();
}

bool TideRegressionTestDriver::preciseCompare(
  double x,
  double y,
  size_t decimalPlaces
) {
  stringstream ss;
  ss << fixed << setprecision(decimalPlaces) << x;
  string s = ss.str();
  ss.str("");
  ss << y;
  return s == ss.str();
}

