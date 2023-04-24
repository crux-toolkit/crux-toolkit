#include <vector>

#include "IndexedMassSpectralPeak.h"

using std::vector;

namespace CruxQuant{
    class PeakIndexingEngine{
    private:
        vector<IndexedMassSpectralPeak> indexedPeaks;
        const int BinsPerDalton = 100;
    public:
        PeakIndexingEngine(vector<IndexedMassSpectralPeak> indexedPeaks);
        void IndexMassSpectralPeaks();
        void cleaIndex();
        vector<IndexedMassSpectralPeak> GetIndexedPeaks(double theorMass, int zeroBasedScanIndex, double tolerance, int chargeState);
        int  BinarySearchForIndexedPeak(vector<IndexedMassSpectralPeak> indexedPeaks, int zeroBasedScanIndex);
    };
}