#ifndef SPECTRUM_MSTOOLKIT_CONVERTER_H
#define SPECTRUM_MSTOOLKIT_CONVERTER_H
#include "SpectrumConverter.h"
using namespace std;
class MSToolkitSpectrumConverter : public SpectrumConverter {
public:
    MSToolkitSpectrumConverter(const std::string &filename);
    bool parse(vector<pb::Spectrum> &all_spectra, int ms_level = 2, bool dia_mode = false) override;
};
#endif