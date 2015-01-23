#ifndef _FFT_H
#define _FFT_H

#include <cmath>

#ifndef PI
#define PI      3.14159265358979323846
#endif

namespace Hardklor {

typedef struct complex{
	double real;
	double imag;
} complex;

};

void BitReverse(Hardklor::complex* data, int size);
void FFT(Hardklor::complex* data, int size, bool forward);
void FFTreal(Hardklor::complex* data, int size);

#endif
