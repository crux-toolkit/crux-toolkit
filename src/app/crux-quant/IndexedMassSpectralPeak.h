#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <functional>

namespace CruxQuant {

    class IndexedMassSpectralPeak {
    public:
        const double mz;
        const double intensity;
        const int zeroBasedMs1ScanIndex;
        const double retentionTime;
        
        IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime);

        std::string ToString() const;
    };

} // namespace CruxQuant
