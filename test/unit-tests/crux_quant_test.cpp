#include <gtest/gtest.h>
#include "crux-quant/Utils.h"


std::string spectrum_file = TEST_DATA_PATH "/test-data/test_filtered.mzML";
std::string psm_file = TEST_DATA_PATH "/test-data/tide-search.txt";

typedef pwiz::msdata::SpectrumListPtr SpectrumListPtr;

TEST(CruxQuant, TestLoadSpectra) {
  SpectrumListPtr spectra_ms1 = CruxQuant::loadSpectra(spectrum_file, 1);
  EXPECT_EQ(spectra_ms1->size(), 13);
 
}

// TODO rewrite this test, make it more robust
TEST(CruxQuant, TestIndexes) {
  SpectrumListPtr spectra_ms1 = CruxQuant::loadSpectra(spectrum_file, 1);
  CruxQuant::IndexedSpectralResults indexResults = CruxQuant::indexedMassSpectralPeaks(spectra_ms1, spectrum_file);
  std::map<int, std::map<int, CruxQuant::IndexedMassSpectralPeak>> indexes = indexResults._indexedPeaks;

  for(const auto& pair: indexes ){
      int index = pair.first;
      const std::map<int, CruxQuant::IndexedMassSpectralPeak>& peaksList = pair.second;
  
      // Loop through the list of IndexedMassSpectralPeak for this index
      for (const auto& peak : peaksList) {
          // Access the class variables of CruxQuant::IndexedMassSpectralPeak for each peak
          FLOAT_T mz = peak.second.mz;
          int roundedMz = static_cast<int>(std::round(mz * CruxQuant::BINS_PER_DALTON));
          EXPECT_EQ(index, roundedMz);
      }
  }
 
}

TEST(CruxQuant, TestCreateIdentifications){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  SpectrumListPtr spectra_ms2 = CruxQuant::loadSpectra(spectrum_file, 2);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader, spectrum_file, spectra_ms2);
  EXPECT_EQ(8700, allIdentifications.size());
}

TEST(CruxQuant, TestcalCulateTheoreticalIsotopeDistributions){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  SpectrumListPtr spectra_ms2 = CruxQuant::loadSpectra(spectrum_file, 2);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader, spectrum_file, spectra_ms2);
  unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxQuant::calculateTheoreticalIsotopeDistributions(allIdentifications);
  EXPECT_EQ(67, modifiedSequenceToIsotopicDistribution.size());
  auto it = modifiedSequenceToIsotopicDistribution.find(allIdentifications.at(0).sequence);
  ASSERT_NE(it, modifiedSequenceToIsotopicDistribution.end());
}

TEST(CruxQuant, TestSetPeakFindingMass){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  SpectrumListPtr spectra_ms2 = CruxQuant::loadSpectra(spectrum_file, 2);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader, spectrum_file, spectra_ms2);
  unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxQuant::calculateTheoreticalIsotopeDistributions(allIdentifications);
  CruxQuant::setPeakFindingMass(allIdentifications, modifiedSequenceToIsotopicDistribution);
  double actual = std::round(244.62680618106413 * 100)/100;
  double expected = std::round(allIdentifications.at(0).peakfindingMass * 100)/100;
  EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(CruxQuant, TestCreateChargeStates){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  SpectrumListPtr spectra_ms2 = CruxQuant::loadSpectra(spectrum_file, 2);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader, spectrum_file, spectra_ms2);
  vector<double> chargeStates = CruxQuant::createChargeStates(allIdentifications);
  EXPECT_EQ(66, chargeStates.size()) << chargeStates.size();
}

// TODO Write the actual tests
TEST(CruxQuant, TestQuantifyMs2IdentifiedPeptides){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  SpectrumListPtr spectra_ms2 = CruxQuant::loadSpectra(spectrum_file, 2);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader, spectrum_file, spectra_ms2);
  
  // CruxQuant::quantifyMs2IdentifiedPeptides(spectrum_file, allIdentifications);
  EXPECT_EQ(1, 1);
}

// TODO Write the actual tests
TEST(CruxQuant, TestGetIndexedPeak){
  EXPECT_EQ(1, 1);
}
