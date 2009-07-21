/**
 * \file CMSReader.h
 * \brief C Wrapper for MS Toolkit MSReader Object.
 */

/* AUTHOR: Sean Joseph McIlwain
 * PROJECT: shared
 * COPYRIGHT: 2009 - University of Washington
 * VERSION: $Revision: 1.0 $
 ********************************************************************/
#include "CSpectrum.h"

#ifdef __cplusplus
#include "MSReader.h"
using namespace MSToolkit;
#endif

//Make MSReader visible to C as an opaque struct.
#ifdef __cplusplus
extern "C" {
typedef MSReader MST_MSREADER_T;
#else
typedef struct MSReader_ MST_MSREADER_T;
#endif

  /**
   * \returns A new MSReader Object.
   */
  MST_MSREADER_T* newMST_MSReader();

  /**
   * frees A MSReader object
   */
  void freeMST_MSReader(MST_MSREADER_T* reader);

  /**
   * calls readFile on the MSReader Object.
   */
  void MST_MSReader_readFile(MST_MSREADER_T *reader , char* filename, MST_SPECTRUM_T* spectrum);

  /**
   * calls readFile on the MSReader Object with a scan number.
   */
  void MST_MSReader_readFileScan(MST_MSREADER_T* reader, char* filename, MST_SPECTRUM_T* spectrum, int scNum);





#ifdef __cplusplus
}
#endif
