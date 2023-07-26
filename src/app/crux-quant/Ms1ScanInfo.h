#pragma once

#include <string>

namespace CruxQuant{
    class Ms1ScanInfo{
    public:
        const int oneBasedScanNumber;
        const int zeroBasedMs1ScanIndex;
        const double retentionTime;

        Ms1ScanInfo(int oneBasedScanNumber, int zeroBasedMs1ScanIndex, double retentionTime);

        std::string ToString() const;
    };
} // namespace CruxQuant
