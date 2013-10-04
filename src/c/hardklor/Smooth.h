#ifndef _SMOOTH_H
#define _SMOOTH_H

#include "MSToolkit/Spectrum.h"

using namespace std;
using namespace MSToolkit;

double SG_GenFact(int, int);
double SG_GramPoly(int, int, int, int);
void SG_Smooth(Spectrum&, int, int);
void SG_SmoothD(float*, int, int, int);
double SG_Weight(int, int, int, int, int);

#endif
