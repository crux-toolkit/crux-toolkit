#pragma once

#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

namespace CruxLFQ {

class IndexedMassSpectralPeak {
   public:
    double mz = 0.0;
    double intensity = 0.0;
    int zeroBasedMs1ScanIndex = 0;
    double retentionTime = 0.0;

    // Default constructor
    IndexedMassSpectralPeak() = default;
    IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime);

    bool operator==(const IndexedMassSpectralPeak& other) const;
    bool operator!=(const IndexedMassSpectralPeak& other) const;
    IndexedMassSpectralPeak& operator=(const IndexedMassSpectralPeak& other);
    int GetHashCode() const;
    std::string ToString() const;
};

}  // namespace CruxLFQ

// Declare the explicit specialization of std::hash before use
namespace std {
template <>
struct hash<CruxLFQ::IndexedMassSpectralPeak> {
    size_t operator()(const CruxLFQ::IndexedMassSpectralPeak& peak) const {
        return peak.GetHashCode();
    }
};
}  // namespace std