#include "SpectrumConverter.h"
using namespace std;
SpectrumConverter::SpectrumConverter(const std::string &filename) : filename_(filename) {
}
std::string SpectrumConverter::GetFilename() const {
    return filename_;
}