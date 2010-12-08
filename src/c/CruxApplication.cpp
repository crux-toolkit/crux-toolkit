/**
 * \file CruxApplication.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Abstract Object for a CruxApplication
 *****************************************************************************/
#include "CruxApplication.h"

#include "carp.h"
#include "parameter.h"

using namespace std;


/**
 * Frees an allocated CruxApplication
 */
CruxApplication::~CruxApplication() {
}

/**
 * \returns the file stem of the application, default blank.
 */
string CruxApplication::getFileStem() {
  return "";
}
