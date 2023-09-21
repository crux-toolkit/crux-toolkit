#include <gtest/gtest.h>
#include "crux-quant/Utils.h"


std::string spectrum_file = TEST_DATA_PATH "/test-data/test_filtered.mzML";
std::string psm_file = TEST_DATA_PATH "/test-data/tide-search.txt";


TEST(CruxQuant, TestLoadSpectra) {
  Crux::SpectrumCollection* spectra_ms1 = CruxQuant::loadSpectra(spectrum_file, 1);
  EXPECT_EQ(spectra_ms1->getNumSpectra(), 13);
 
}

TEST(CruxQuant, TestIndexes) {
  Crux::SpectrumCollection* spectra_ms1 = CruxQuant::loadSpectra(spectrum_file, 1);
  std::unordered_map<int, std::vector<CruxQuant::IndexedMassSpectralPeak>> indexes =  CruxQuant::indexedMassSpectralPeaks(spectra_ms1);

  for(const auto& pair: indexes ){
      int index = pair.first;
      const std::vector<CruxQuant::IndexedMassSpectralPeak>& peaksList = pair.second;
  
      // Loop through the list of IndexedMassSpectralPeak for this index
      for (const auto& peak : peaksList) {
          // Access the class variables of CruxQuant::IndexedMassSpectralPeak for each peak
          FLOAT_T mz = peak.mz;
          int roundedMz = static_cast<int>(std::round(mz * CruxQuant::BINS_PER_DALTON));
          EXPECT_EQ(index, roundedMz);
      }
  }
 
}

TEST(CruxQuant, TestCreateIdentifications){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader);
  EXPECT_EQ(100, allIdentifications.size());
}

TEST(CruxQuant, TestcalCulateTheoreticalIsotopeDistributions){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader);
  unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxQuant::calculateTheoreticalIsotopeDistributions(allIdentifications);
  EXPECT_EQ(67, modifiedSequenceToIsotopicDistribution.size());
  auto it = modifiedSequenceToIsotopicDistribution.find(allIdentifications.at(0).Sequence);
  ASSERT_NE(it, modifiedSequenceToIsotopicDistribution.end());
}

TEST(CruxQuant, TestSetPeakFindingMass){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader);
  unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxQuant::calculateTheoreticalIsotopeDistributions(allIdentifications);
  CruxQuant::SetPeakFindingMass(allIdentifications, modifiedSequenceToIsotopicDistribution);
  double actual = std::round(244.62680618106413 * 100)/100;
  double expected = std::round(allIdentifications.at(0).PeakfindingMass * 100)/100;
  EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(CruxQuant, TestCreateChargeStates){
  MatchFileReader *matchFileReader = new MatchFileReader(psm_file);
  vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader);
  vector<double> chargeStates = CruxQuant::createChargeStates(allIdentifications);
  EXPECT_EQ(66, chargeStates.size()) << chargeStates.size();
}