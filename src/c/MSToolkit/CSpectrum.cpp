/*************************************************************************//**
 * \file CSpectrum.cpp
 * $Revision: 1.0 $
 * \brief: C Wrapper for MSToolkit Spectrum Object.
 ****************************************************************************/
#include "CSpectrum.h"
#include <iostream>

using namespace std;
using namespace MSToolkit;


/**
 * \returns A new Spectrum Object.
 */
MST_SPECTRUM_T* newMST_Spectrum() {
  return new Spectrum();
}

/**
 * frees A Spectrum object
 */
void freeMST_Spectrum(MST_SPECTRUM_T* spectrum) {
  delete spectrum;
}

/**
 * get the precursor mz of the spectrum object.
 */
double MST_Spectrum_getMZ(MST_SPECTRUM_T* spectrum) {
  return spectrum->getMZ();
}

/**
 * get the scan number of the spectrum object.
 */
int MST_Spectrum_getScanNumber(MST_SPECTRUM_T* spectrum) {
  return spectrum->getScanNumber();
}

/**
 * get the charge from the Z line
 */
int MST_Spectrum_atZ(MST_SPECTRUM_T* spectrum, int z_idx) {
  return spectrum->atZ(z_idx).z;
}

/**
 * get the number of possible charges from the z line
 */
int MST_Spectrum_sizeZ(MST_SPECTRUM_T* spectrum) {
  return spectrum->sizeZ();
}

/**
 * get the size / or number of peaks in the spectrum
 */
int MST_Spectrum_size(MST_SPECTRUM_T* spectrum) {
  return spectrum->size();
}

/**
 * get the intensity of the peak
 */
double MST_Spectrum_intensity(MST_SPECTRUM_T* spectrum, int peak_idx) {
  return spectrum->at(peak_idx).intensity;
}

/**
 * get the mz of the peak
 */
double MST_Spectrum_mz(MST_SPECTRUM_T* spectrum, int peak_idx) {
  return spectrum->at(peak_idx).mz;
}

