#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include "IndexedMassSpectralPeak.h"
#include "Ms1ScanInfo.h"

using std::string;
using std::unordered_map;
using std::vector;

namespace CruxLFQ {

class LFQMetaData {
   private:
    vector<vector<IndexedMassSpectralPeak>*>* indexedPeaks;
    unordered_map<string, vector<Ms1ScanInfo>>* ms1Scans;

    LFQMetaData() {}                                      // Private constructor to prevent instantiation
    LFQMetaData(const LFQMetaData&) = delete;             // Disable copy constructor
    LFQMetaData& operator=(const LFQMetaData&) = delete;  // Disable assignment operator

   public:
    static LFQMetaData& getInstance() {
        static LFQMetaData instance;
        return instance;
    }

    void setIndexedPeaks(vector<vector<IndexedMassSpectralPeak>*>* peaks) {
        indexedPeaks = peaks;
    }

    void setMs1Scans(unordered_map<string, vector<Ms1ScanInfo>>* scans) {
        ms1Scans = scans;
    }

    vector<vector<IndexedMassSpectralPeak>*>* getIndexedPeaks() {
        return indexedPeaks;
    }

    unordered_map<string, vector<Ms1ScanInfo>>* getMs1Scans() {
        return ms1Scans;
    }
};
}  // namespace CruxLFQ