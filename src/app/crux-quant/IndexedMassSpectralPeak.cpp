#include "IndexedMassSpectralPeak.h"
#include <sstream>

namespace CruxQuant {

    IndexedMassSpectralPeak::IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime) :
        mz(mz), zeroBasedMs1ScanIndex(zeroBasedMs1ScanIndex), retentionTime(retentionTime), intensity(intensity) {}

   
    std::string IndexedMassSpectralPeak::ToString() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << "MZ Value: " << mz << "; " << "MS1 Scan Index: "<< zeroBasedMs1ScanIndex << "; " << "Intensity: " << intensity << "; " << "Rentention Time: " << retentionTime ;
        return oss.str();
    }

} // namespace CruxQuant
