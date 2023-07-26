#include "Ms1ScanInfo.h"
#include <sstream>

namespace CruxQuant{
    Ms1ScanInfo::Ms1ScanInfo(int oneBasedScanNumber, int zeroBasedMs1ScanIndex, double retentionTime)
        : oneBasedScanNumber(oneBasedScanNumber),
          zeroBasedMs1ScanIndex(zeroBasedMs1ScanIndex),
          retentionTime(retentionTime){}

    std::string Ms1ScanInfo::ToString() const{
        std::stringstream ss;
        ss << zeroBasedMs1ScanIndex << "; " << oneBasedScanNumber << "; " << retentionTime;
        return ss.str();
    }
} // namespace CruxQuant
