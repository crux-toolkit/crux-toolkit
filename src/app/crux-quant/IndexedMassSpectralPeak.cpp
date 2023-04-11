#include "IndexedMassSpectralPeak.h"

#include <string>

using std::string;

using namespace CruxQuant;

IndexedMassSpectralPeak::IndexedMassSpectralPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime){
    this->mz = mz;
    this->intensity = intensity;
    this->zeroBasedMs1ScanIndex = zeroBasedMs1ScanIndex;
    this->retentionTime = retentionTime;
}

IndexedMassSpectralPeak::~IndexedMassSpectralPeak(void){}

double IndexedMassSpectralPeak::getMz() const{
    return this->mz;
}

double IndexedMassSpectralPeak::getIntensity() const{
    return this->intensity;
}

int IndexedMassSpectralPeak::getZeroBasedMs1ScanIndex() const{
    return this->zeroBasedMs1ScanIndex;
}

double IndexedMassSpectralPeak::getRetentionTime() const{
    return this->retentionTime;
}

void IndexedMassSpectralPeak::setMz(double mz){
    this->mz = mz;
}

void IndexedMassSpectralPeak::setIntensity(double intensity){
    this->intensity = intensity;
}

void IndexedMassSpectralPeak::setZeroBasedMs1ScanIndex(int zeroBasedMs1ScanIndex){
    this->zeroBasedMs1ScanIndex = zeroBasedMs1ScanIndex;
}

void IndexedMassSpectralPeak::setRetentionTime(double retentionTime){
    this->retentionTime = retentionTime;
}

bool IndexedMassSpectralPeak::operator==(const IndexedMassSpectralPeak& other) const{
    
    return &other != nullptr && 
        this->mz == other.mz &&
        this->zeroBasedMs1ScanIndex == other.zeroBasedMs1ScanIndex;
}

bool IndexedMassSpectralPeak::operator!=(const IndexedMassSpectralPeak& other) const{
    return !(*this == other);
}

int IndexedMassSpectralPeak::GetHashCode(){
    return std::hash<double>()(this->mz);
}

string IndexedMassSpectralPeak::ToString() const{
    return "mz: " + std::to_string(this->mz) + ", intensity: " + std::to_string(this->intensity) + ", zeroBasedMs1ScanIndex: " + std::to_string(this->zeroBasedMs1ScanIndex) + ", retentionTime: " + std::to_string(this->retentionTime);
}
