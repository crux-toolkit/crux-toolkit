#include <cppunit/config/SourcePrefix.h>
#include "TestSpectrum.h"
#include "Spectrum.h"
#include "SpectrumZState.h"
#include "parameter.h" 

using namespace std;
using namespace Crux;

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
  CPPUNIT_ASSERT(default_s->getNumZStates() == 0);

  const vector<SpectrumZState>& zstates = default_s->getZStates(); 
  CPPUNIT_ASSERT(zstates.size() == 0);
  vector <SpectrumZState> zstates_to_search = default_s->getZStatesToSearch();
  CPPUNIT_ASSERT(zstates_to_search.size() == 0);

  CPPUNIT_ASSERT(default_s->getMinPeakMz() == 0);
  CPPUNIT_ASSERT(default_s->getMaxPeakMz() == 0);
  CPPUNIT_ASSERT(default_s->getNumPeaks() == 0);
  CPPUNIT_ASSERT(default_s->getTotalEnergy() == 0);
  CPPUNIT_ASSERT(default_s->getMaxPeakIntensity() == -1);
  //CPPUNIT_ASSERT(default_s->getMass(1) == 0);
  //CPPUNIT_ASSERT(default_s->getNeutralMass(1) < 0);
  //CPPUNIT_ASSERT(default_s->getSinglyChargedMass(1) == 0);

}










