#pragma once

#include "Ms1ScanInfo.h"

namespace CruxLFQ {
struct Ms1ScanInfo {
    int oneBasedScanNumber;
    int zeroBasedMs1ScanIndex;
    double retentionTime;

    Ms1ScanInfo(int oneBasedScan, int zeroBasedIndex, double rtime)
        : oneBasedScanNumber(oneBasedScan), zeroBasedMs1ScanIndex(zeroBasedIndex), retentionTime(rtime) {}
};
}  // namespace CruxLFQ