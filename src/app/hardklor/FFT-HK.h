#ifndef _FFT_HK_H
#define _FFT_HK_H

#include "MSReader.h"
#include "FFT.h"
#include <cmath>
#include <vector>

void FFTCharge(double *f, MSToolkit::Spectrum& s, unsigned int start, unsigned int stop,
							 unsigned int lowCharge, unsigned int highCharge, double interval, bool bSpline=false);
double GetIntensity(MSToolkit::Spectrum& s, unsigned int start, unsigned int stop, double mz);
void Patterson(double *f, MSToolkit::Spectrum& s, unsigned int start, unsigned int stop,
							 unsigned int lowCharge, unsigned int highCharge, double interval);
void SenkoCharge(std::vector<int> *charges, MSToolkit::Spectrum& s, unsigned int start, unsigned int stop, 
								 unsigned int lowCharge, unsigned int highCharge, double interval, char method);

#endif
