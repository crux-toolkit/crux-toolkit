#ifndef _FFT_HK_H
#define _FFT_HK_H

#include "MSReader.h"
#include "FFT.h"
#include <cmath>
#include <vector>

using namespace std;
using namespace MSToolkit;

void FFTCharge(double *f, Spectrum& s, unsigned int start, unsigned int stop,
							 unsigned int lowCharge, unsigned int highCharge, double interval, bool bSpline=false);
double GetIntensity(Spectrum& s, unsigned int start, unsigned int stop, double mz);
void Patterson(double *f, Spectrum& s, unsigned int start, unsigned int stop,
	 unsigned int lowCharge, unsigned int highCharge/*, double interval*/);
void SenkoCharge(vector<int> *charges, Spectrum& s, unsigned int start, unsigned int stop, 
								 unsigned int lowCharge, unsigned int highCharge, double interval, char method);

#endif
