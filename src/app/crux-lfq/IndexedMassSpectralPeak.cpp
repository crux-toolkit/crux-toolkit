#include "IndexedMassSpectralPeak.h"

#include <sstream>

namespace CruxLFQ {

IndexedMassSpectralPeak::IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime) : mz(mz), zeroBasedMs1ScanIndex(zeroBasedMs1ScanIndex), retentionTime(retentionTime), intensity(intensity) {}

std::string IndexedMassSpectralPeak::ToString() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << "MZ Value: " << mz << "; "
        << "MS1 Scan Index: " << zeroBasedMs1ScanIndex << "; "
        << "Intensity: " << intensity << "; "
        << "Rentention Time: " << retentionTime;
    return oss.str();
}

bool IndexedMassSpectralPeak::operator==(const IndexedMassSpectralPeak& other) const {
    return &other != nullptr &&
           this->mz == other.mz &&
           this->zeroBasedMs1ScanIndex == other.zeroBasedMs1ScanIndex;
}

bool IndexedMassSpectralPeak::operator!=(const IndexedMassSpectralPeak& other) const {
    return !(*this == other);
}

int IndexedMassSpectralPeak::GetHashCode() const {
    return std::hash<double>()(this->mz);
}

IndexedMassSpectralPeak& IndexedMassSpectralPeak::operator=(const IndexedMassSpectralPeak& other) {
    if (this != &other) {
        mz = other.mz;
        intensity = other.intensity;
        zeroBasedMs1ScanIndex = other.zeroBasedMs1ScanIndex;
        retentionTime = other.retentionTime;
    }
    return *this;
}

}  // namespace CruxLFQ
