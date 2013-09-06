#ifndef TIDEREGRESSIONTESTDRIVER_H
#define TIDEREGRESSIONTESTDRIVER_H

#include <string>
#include <vector>

struct TideRegressionSettings {
public:
  TideRegressionSettings():
    fasta(""), enzyme("trypsin"), digestion("full-digest"), missedCleavages(0),
    minLength(6), maxLength(50), minMass(200.0), maxMass(7200.0),
    monoisotopicPrecursor(false), modsSpec("C+57.02146"), spectrumRecords(""), 
    massWindow(3.0), topMatch(5)
    {}
  TideRegressionSettings(const std::string& vFasta, const std::string& vEnzyme,
    const std::string& vDigestion, int vMissedCleavages, int vMinLength,
    int vMaxLength, double vMinMass, double vMaxMass, bool vMonoisotopicPrecursor,
    const std::string& vModsSpec, const std::string& vSpectrumRecords,
    double vMassWindow, int vTopMatch):
    fasta(vFasta), enzyme(vEnzyme), digestion(vDigestion),
    missedCleavages(vMissedCleavages), minLength(vMinLength), maxLength(vMaxLength),
    minMass(vMinMass), maxMass(vMaxMass),
    monoisotopicPrecursor(vMonoisotopicPrecursor), modsSpec(vModsSpec),
    spectrumRecords(vSpectrumRecords), massWindow(vMassWindow), topMatch(vTopMatch)
    {}
  ~TideRegressionSettings() {}
  // tide-index settings
  std::string fasta;
  std::string enzyme;
  std::string digestion;
  int missedCleavages;
  int minLength, maxLength;
  double minMass, maxMass;
  bool monoisotopicPrecursor;
  std::string modsSpec;
  // tide-search settings
  std::string spectrumRecords;
  double massWindow;
  int topMatch;
};

class TideRegressionTestDriver {
public:
  TideRegressionTestDriver(const TideRegressionSettings& settings);
  TideRegressionTestDriver();
  ~TideRegressionTestDriver();

  bool runTest(
    const std::string& tideIndexPath,
    const std::string& tideSearchPath,
    const std::string& cruxPath
  );
  static void cleanTestFiles();

  static bool isRegularFile(const std::string& path);
  static bool isDir(const std::string& path);

  std::string getError();

  void loadSettings(const TideRegressionSettings& settings);
  void setFasta(const std::string& fasta);
  void setEnzyme(const std::string& enzyme);
  void setDigestion(const std::string& digestion);
  void setMissedCleavages(int missedCleavages);
  void setMinLength(int minLength);
  void setMaxLength(int maxLength);
  void setMinMass(double minMass);
  void setMaxMass(double maxMass);
  void setMonoisotopicPrecursor(bool monoisotopicPrecursor);
  void setModsSpec(const std::string& modsSpec);
  void setSpectrumRecords(const std::string& spectrumRecords);
  void setMassWindow(double massWindow);
  void setTopMatch(int topMatch);

private:
  const static std::string TIDEPEPIX;
  const static std::string TIDEPROTIX;
  const static std::string TIDEAUXLOCS;
  const static std::string TIDEOUT;
  const static std::string CRUXIDX;
  const static std::string CRUXOUT;
  enum COMPARE_TYPE { EXACT, NUMERIC, PEPTIDE };

  TideRegressionSettings settings_;
  std::string error_;

  static std::vector<std::string> splitTab(const std::string& row);
  bool compareResults();
  static bool valueCompare(
    const std::string& v1,
    const std::string& v2,
    COMPARE_TYPE type
  );
  static bool peptideCompare(const std::string& v1, const std::string& v2);
  static bool numCompare(const std::string& v1, const std::string& v2);
  static std::string makePrecise(const std::string& v, size_t decimalPlaces);
  static bool preciseCompare(double x, double y, size_t decimalPlaces);
};

#endif

