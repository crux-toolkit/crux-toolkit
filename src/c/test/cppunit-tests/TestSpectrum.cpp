#include <cppunit/config/SourcePrefix.h>
#include "TestSpectrum.h"
#include "Spectrum.h"
#include "parameter.h" 

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( TestSpectrum );

void TestSpectrum::setUp(){
  initialize_parameters();  
  // initialize variables to use for testing
  default_s = new Spectrum();
}

void TestSpectrum::tearDown(){
  delete default_s;
}

void TestSpectrum::defaultGetters(){
  // available tests include 
  CPPUNIT_ASSERT(default_s->getFirstScan() == 0);
  CPPUNIT_ASSERT(default_s->getPrecursorMz() == 0);
  CPPUNIT_ASSERT(default_s->getNumPossibleZ() == 0);

  const vector<int>& charges = default_s->getPossibleZ(); 
  CPPUNIT_ASSERT(charges.size() == 0);
  vector <int> charges_to_search = default_s->getChargesToSearch();
  CPPUNIT_ASSERT(charges_to_search.size() == 0);

  CPPUNIT_ASSERT(default_s->getMinPeakMz() == 0);
  CPPUNIT_ASSERT(default_s->getMaxPeakMz() == 0);
  CPPUNIT_ASSERT(default_s->getNumPeaks() == 0);
  CPPUNIT_ASSERT(default_s->getTotalEnergy() == 0);
  CPPUNIT_ASSERT(default_s->getMaxPeakIntensity() == -1);
  CPPUNIT_ASSERT(default_s->getMass(1) == 0);
  CPPUNIT_ASSERT(default_s->getNeutralMass(1) < 0);
  CPPUNIT_ASSERT(default_s->getSinglyChargedMass(1) == 0);

}










