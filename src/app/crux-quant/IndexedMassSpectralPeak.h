#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <functional>

namespace CruxQuant {

    class IndexedMassSpectralPeak {
    public:
        const int zeroBasedMs1ScanIndex;
        const double mz;
        const double retentionTime;
        const double intensity;

        IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime);

        bool operator==(const IndexedMassSpectralPeak& otherPeak) const;
        bool operator!=(const IndexedMassSpectralPeak& otherPeak) const;
        size_t GetHashCode() const;
        std::string ToString() const;
    };

} // namespace CruxQuant

namespace std {
    template<>
    struct hash<CruxQuant::IndexedMassSpectralPeak> {
        size_t operator()(const CruxQuant::IndexedMassSpectralPeak& peak) const;
    };
}
