#include "LibraSummary.hpp"

/**
* constructor
* @param pointer to condition instance
*/
LibraSummary::LibraSummary( LibraConditionHandler* pCond) {
  pLibraConditionHandler = pCond;
}

LibraSummary::~LibraSummary() {
}


/**
* @return mass tolerance
*/
double LibraSummary::getMassTolerance() {
  return pLibraConditionHandler->getTolerance();
}


/**
* @return 0 for no centroiding, 1 for average, 2 for weighted average
*/
int LibraSummary::getCentroidingPref() {
  return pLibraConditionHandler->getCentroidingPref();
}


/**
* @return number of iterations used in centroiding
*/
int LibraSummary::getNumCentroidingIterations() {
    return pLibraConditionHandler->getNumCentroidingIterations();
}


/**
*  0:      Normalize against most intense,
* -2:      Normalize against TIC,
* default: Normalize on one isotope
* @return 0, -2, or default (?)
*/
int LibraSummary::getNormalization() {
  return pLibraConditionHandler->getNormalPosition();
}


int LibraSummary::getOutputType() {
  return pLibraConditionHandler->getOutputPrefs();
}


vector<double> LibraSummary::getFragmentMassValues(int targetID) {
  int scanNum = targetID;
  vvf masses = pLibraConditionHandler->getMassIsotopes();
  vector<double> tmp;

  for (int i = 0; i < (int) masses.size(); i++) {
    tmp.push_back( masses[scanNum][i] );
  }

  return tmp;
}
