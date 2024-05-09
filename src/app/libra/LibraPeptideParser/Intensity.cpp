#include "Intensity.hpp"
#include <cmath>

using namespace mzParser;

/**
 * Constructs an Intensity object with existing LibraConditionHandler
 * object and start and end scan numbers.  It turns off Quantitation's
 * print to outfile and print to stdout.
 * @param pCond is a pointer to a LibraConditionHandler object
 * @param startScanNum is first scan number to sum over
 * @param endScanNum is last scan number to sum over
 * @param mzXMLFile is name of mzXML file
 */
Intensity::Intensity(LibraConditionHandler* pCond, char* mzXMLFile, RAMPFILE* FI, int scan) {
  if (dbug) cout << "Intensity::Intensity(LibraConditionHandler*,char*,int)" << endl;

  scanNum = scan;
  mzXMLFileName = mzXMLFile;
  pLibraConditionHandler = pCond;
  pQuantitation = new Quantitation(pLibraConditionHandler, mzXMLFile, FI);

  bool b = false;
  pQuantitation->setWriteToOutFile(b);
  pQuantitation->setPrintToStdOut(b);
  pQuantitation->calculateIntensities(scan, scan);
}


Intensity::Intensity() {
  pQuantitation = NULL;
}

Intensity::~Intensity() {
  delete pQuantitation;
}


int Intensity::getScanNumber() {
  return scanNum;
}


/**
 * get vector of target masses using key scanNumber
 * @param scan number
 * @return target masses
 */
vector<double> Intensity::getTargetMass() {
  return pQuantitation->getTargetMasses(scanNum);
}


int Intensity::getNumberOfChannels() {
  return pQuantitation->getNumberOfChannels(scanNum);
}


/**
 * get vector of absolute intensities using key scanNumber
 * @return absolute intensities
 */
vector<double> Intensity::getAbsoluteIntensities() {
  return pQuantitation->getAbsoluteIntensities(scanNum);
}


vector<double> Intensity::getNormalizedIntensities() {
  return pQuantitation->getNormalizedIntensities(scanNum);
}
