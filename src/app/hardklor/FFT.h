#ifndef _FFT_H
#define _FFT_H

#include <cmath>

#ifndef PI
#define PI      3.14159265358979323846
#endif

typedef struct hardklor_complex{
	double real;
	double imag;
} hardklor_complex;


void BitReverse(hardklor_complex* data, int size);
void FFT(hardklor_complex* data, int size, bool forward);
void FFTreal(hardklor_complex* data, int size);

#endif
