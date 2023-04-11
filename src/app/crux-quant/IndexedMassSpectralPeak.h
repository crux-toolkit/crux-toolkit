#pragma once

#include <string>

using std::string;


namespace CruxQuant{
    class IndexedMassSpectralPeak{
    private:
        double mz;
        double intensity;
        int zeroBasedMs1ScanIndex;
        double retentionTime;
    public:
        IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime);
        ~IndexedMassSpectralPeak(void);

        double getMz() const;
        double getIntensity() const;
        int getZeroBasedMs1ScanIndex() const;
        double getRetentionTime() const;

        void setMz(double mz);
        void setIntensity(double intensity);
        void setZeroBasedMs1ScanIndex(int zeroBasedMs1ScanIndex);
        void setRetentionTime(double retentionTime);
        
        bool operator==(const IndexedMassSpectralPeak& other) const;
        bool operator!=(const IndexedMassSpectralPeak& other) const;
        int GetHashCode();
        string ToString() const;
    };
}