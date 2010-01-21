/**
 * \file CSpectrum.h
 * \brief C Wrapper for MSToolkit Specturm Object.
 */

/* AUTHOR: Sean Joseph McIlwain
 * PROJECT: shared
 * COPYRIGHT: 2009 - University of Washington
 * VERSION: $Revision: 1.0 $
 ********************************************************************/

#ifndef _CSPECTRUM_H_
#define _CSPECTRUM_H_

#ifdef __cplusplus
#include "Spectrum.h"
#endif

//Make Spectrum visible to C as an opaque struct.
#ifdef __cplusplus
typedef MSToolkit::Spectrum MST_SPECTRUM_T;
extern "C" {

#else
typedef struct Spectrum_ MST_SPECTRUM_T;
#endif

  /**
   * \returns A new Spectrum Object.
   */
  MST_SPECTRUM_T* newMST_Spectrum();

  /**
   * frees A Spectrum object
   */
  void freeMST_Spectrum(MST_SPECTRUM_T* spectrum);

  /**
   * get the precursor mz of the spectrum object.
   */
  double MST_Spectrum_getMZ(MST_SPECTRUM_T* spectrum);

  /**
   * get the scan number of the spectrum object.
   */
  int MST_Spectrum_getScanNumber(MST_SPECTRUM_T* spectrum);

  /**
   * get the charge from the Z line
   */
  int MST_Spectrum_atZ(MST_SPECTRUM_T* spectrum, int z_idx);

  /**
   * get the number of possible charges from the z line
   */
  int MST_Spectrum_sizeZ(MST_SPECTRUM_T* spectrum);

  /**
   * get the size / or number of peaks in the spectrum
   */
  int MST_Spectrum_size(MST_SPECTRUM_T* spectrum);

  /**
   * get the intensity of the peak
   */
  double MST_Spectrum_intensity(MST_SPECTRUM_T* spectrum, int peak_idx);

  /**
   * get the mz of the peak
   */
  double MST_Spectrum_mz(MST_SPECTRUM_T* spectrum, int peak_idx);

#ifdef __cplusplus
}
#endif

#endif
