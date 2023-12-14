#pragma once

#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

namespace CruxLFQ {

class IndexedMassSpectralPeak {
   public:
    const double mz;
    const double intensity;
    const int zeroBasedMs1ScanIndex;
    const double retentionTime;

    // Default constructor
    IndexedMassSpectralPeak()
        : mz(0.0),                   // Initialize mz to a suitable default
          intensity(0.0),            // Initialize intensity to a suitable default
          zeroBasedMs1ScanIndex(0),  // Initialize zeroBasedMs1ScanIndex to a suitable default
          retentionTime(0.0)         // Initialize retentionTime to a suitable default
    {}

    IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime);

    bool operator==(const IndexedMassSpectralPeak& other) const;
    bool operator!=(const IndexedMassSpectralPeak& other) const;
    IndexedMassSpectralPeak& operator=(const IndexedMassSpectralPeak& other);
    int GetHashCode() const;
    std::string ToString() const;
};

}  // namespace CruxQuant


// Declare the explicit specialization of std::hash before use
namespace std {
    template <>
    struct hash<CruxLFQ::IndexedMassSpectralPeak> {
        size_t operator()(const CruxLFQ::IndexedMassSpectralPeak& peak) const {
            return peak.GetHashCode();
        }
    };
}