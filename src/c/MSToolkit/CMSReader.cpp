/*************************************************************************//**
 * \file CMSReader.cpp
 * $Revision: 1.0 $
 * \brief: C Wrapper for MSToolkit MSReader Object.
 ****************************************************************************/

#include "CMSReader.h"



/**
 * \returns A new MSReader Object.
 */
MST_MSREADER_T* newMST_MSReader() {
  return new MST_MSREADER_T();
}

/**
 * frees A MSReader object
 */
void freeMST_MSReader(MST_MSREADER_T* reader) {
  delete reader;
}

/**
 * calls readFile on the MSReader Object.
 */
void MST_MSReader_readFile(MST_MSREADER_T *reader, 
			  char* filename, 
			  MST_SPECTRUM_T* spectrum) {
  reader->readFile(filename, *spectrum);
}

/**
 * calls readFile on the MSReader Object with a scan number.
 */
void MST_MSReader_readFileScan(MST_MSREADER_T* reader, 
			      char* filename, 
			      Spectrum* spectrum, 
			      int scNum) {
  reader->readFile(filename, *spectrum, scNum);
}
