#ifndef SPECTRUM_CONVERTER_H
#define SPECTRUM_CONVERTER_H
#include "spectrum.pb.h"
using namespace std;
class SpectrumConverter {
public:
    SpectrumConverter(const std::string &filename);
    virtual bool parse(std::vector<pb::Spectrum> &all_spectra, int ms_level=2, bool dia_mode = false);
    std::string GetFilename() const;
private:
    std::string filename_;
};
#endif