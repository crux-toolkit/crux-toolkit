#ifndef _S2N_H
#define _S2N_H

#include "MSToolkit/Spectrum.h"
#include "Smooth.h"
#include <algorithm>

Spectrum signalToNoise(Spectrum& s, int start, int stop, float sig, float* cutoff, bool skipZero=true, bool subtract=false);
Spectrum SNPeaks(Spectrum& s, int start, int stop, float SN, float FWHM, float max, float base);
Spectrum SNSubtracted(Spectrum& s, int start, int stop, float SN, float FWHM, float max, float base);
float findSNCutoff(Spectrum& s, int start, int stop, float sig, bool skipZero=true);
float findSNCutoff2(Spectrum& s, int start, int stop, float sig, double& max, bool skipZero=true);

#endif
