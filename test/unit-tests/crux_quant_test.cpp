#include <gtest/gtest.h>
#include "CruxQuantApplication.h" // The import of this header has side effects


std::string spectrum_file = TEST_DATA_PATH "/test-data/sliced-mzml.mzML";


TEST(TestIndexes, BasicAssertions) {
  // Expect two strings not to be equal.
  Crux::SpectrumCollection* spectra_ms1 = CruxQuantApplication::loadSpectra(spectrum_file, 1);
  std::unordered_map<int, std::list<CruxQuant::IndexedMassSpectralPeak>> indexes =  CruxQuantApplication::IndexedMassSpectralPeaks(spectra_ms1);

  for(const auto& pair: indexes ){
      int index = pair.first;
      const std::list<CruxQuant::IndexedMassSpectralPeak>& peaksList = pair.second;
  
      // Loop through the list of IndexedMassSpectralPeak for this index
      for (const auto& peak : peaksList) {
          // Access the class variables of CruxQuant::IndexedMassSpectralPeak for each peak
          FLOAT_T mz = peak.mz;
          int roundedMz = static_cast<int>(std::round(mz * CruxQuantApplication::BinsPerDalton));
          EXPECT_EQ(index, roundedMz);
      }
  }
 
}


TEST(LoadSpectra, BasicAssertions) {
  // Expect two strings not to be equal.
  Crux::SpectrumCollection* spectra_ms1 = CruxQuantApplication::loadSpectra(spectrum_file, 1);
  EXPECT_EQ(spectra_ms1->getNumSpectra(), 10);
 
}
