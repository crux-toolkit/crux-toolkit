#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include "IndexedMassSpectralPeak.h"

using std::string;
using std::unordered_map;
using std::vector;

namespace CruxLFQ {
class LFQMetaData {
   private:
    std::vector<std::vector<CruxLFQ::IndexedMassSpectralPeak>*>* indexedPeaks;

    LFQMetaData() {}                                      // Private constructor to prevent instantiation
    LFQMetaData(const LFQMetaData&) = delete;             // Disable copy constructor
    LFQMetaData& operator=(const LFQMetaData&) = delete;  // Disable assignment operator

   public:
    static LFQMetaData& getInstance() {
        static LFQMetaData instance;
        return instance;
    }

    void setIndexedPeaks(std::vector<std::vector<CruxLFQ::IndexedMassSpectralPeak>*>* peaks) {
        indexedPeaks = peaks;
    }

    vector<vector<IndexedMassSpectralPeak>*>* getIndexedPeaks() {
        return indexedPeaks;
    }
};
}  // namespace CruxLFQ