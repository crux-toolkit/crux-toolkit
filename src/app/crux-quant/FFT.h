#ifndef _FFT_H
#define _FFT_H

#include <cmath>

#ifndef PI
#define PI      3.14159265358979323846
#endif

typedef struct mercury_complex{
	double real;
	double imag;
} mercury_complex;


void BitReverse(mercury_complex* data, int size);
void FFT(mercury_complex* data, int size, bool forward);
void FFTreal(mercury_complex* data, int size);

#endif
