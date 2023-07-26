#include "IndexedMassSpectralPeak.h"
#include <sstream>

namespace CruxQuant {

    IndexedMassSpectralPeak::IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime) :
        mz(mz), zeroBasedMs1ScanIndex(zeroBasedMs1ScanIndex), retentionTime(retentionTime), intensity(intensity) {}

    bool IndexedMassSpectralPeak::operator==(const IndexedMassSpectralPeak& otherPeak) const {
        return otherPeak.mz == mz && otherPeak.zeroBasedMs1ScanIndex == zeroBasedMs1ScanIndex;
    }

    bool IndexedMassSpectralPeak::operator!=(const IndexedMassSpectralPeak& otherPeak) const {
        return !(*this == otherPeak);
    }

    size_t IndexedMassSpectralPeak::GetHashCode() const {
        size_t hash1 = std::hash<double>()(mz);
        size_t hash2 = std::hash<int>()(zeroBasedMs1ScanIndex);
        return hash1 ^ (hash2 << 1);
    }

    std::string IndexedMassSpectralPeak::ToString() const {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << mz << "; " << zeroBasedMs1ScanIndex;
        return oss.str();
    }

} // namespace CruxQuant

namespace std {
    size_t hash<CruxQuant::IndexedMassSpectralPeak>::operator()(const CruxQuant::IndexedMassSpectralPeak& peak) const {
        return peak.GetHashCode();
    }
}
