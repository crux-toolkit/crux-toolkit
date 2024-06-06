#include "LibraWrapper.hpp"
using namespace mzParser;

/**
 * constructor for main.cpp where LibraConditionHandler has been filled prior
 * to constructor.  
 * @param pCond pointer to LibraConditionHandler instance
 * @param mzXMLFileName name of mzXML file
 */
LibraWrapper::LibraWrapper(LibraConditionHandler* pCond, char* mzXMLFileName, RAMPFILE* FI) {
  pFI = FI;
  mzXMLFile = mzXMLFileName;
  pLibraConditionHandler = pCond;
}


/**
 * constructor for case when LibraConditionHandler has not been filled ahead of time.
 * @param conditionFileName name condition xml file
 * @param mzXMLFileName name of mzXML file
 */
LibraWrapper::LibraWrapper(char* conditionFileName, char* mzXMLFileName, RAMPFILE* FI) {
  mzXMLFile = mzXMLFileName;
  pFI = FI;
  pLibraConditionHandler = new LibraConditionHandler();
  pLibraConditionHandler->setFileName( conditionFileName );
  pLibraConditionHandler->readFile();
}


LibraWrapper::~LibraWrapper() {
}


/**
 * get a pointer to a LibraSummary instance
 * @return pointer to a LibraSummary object
 */
LibraSummary* LibraWrapper::getLibraSummary() {
  LibraSummary* pLibraSummary = new LibraSummary(pLibraConditionHandler);
  return pLibraSummary;
}


/**
 * get a pointer to a LibraResult instance
 * @return pointer to a LibraResult object
 */
LibraResult* LibraWrapper::getLibraResult(int scanNum) {
  LibraResult* pLibraResult = new LibraResult(pLibraConditionHandler, mzXMLFile, pFI, scanNum);
  return pLibraResult;
}

