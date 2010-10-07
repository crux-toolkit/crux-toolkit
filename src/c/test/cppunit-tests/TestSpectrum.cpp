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
  CPPUNIT_ASSERT(default_s->get_first_scan() == 0);
  CPPUNIT_ASSERT(default_s->get_precursor_mz() == 0);
  CPPUNIT_ASSERT(default_s->get_num_possible_z() == 0);

  const vector<int>& charges = default_s->get_possible_z(); 
  CPPUNIT_ASSERT(charges.size() == 0);
  vector <int> charges_to_search = default_s->get_charges_to_search();
  CPPUNIT_ASSERT(charges_to_search.size() == 0);

  CPPUNIT_ASSERT(default_s->get_min_peak_mz() == 0);
  CPPUNIT_ASSERT(default_s->get_max_peak_mz() == 0);
  CPPUNIT_ASSERT(default_s->get_num_peaks() == 0);
  CPPUNIT_ASSERT(default_s->get_total_energy() == 0);
  CPPUNIT_ASSERT(default_s->get_max_peak_intensity() == -1);
  CPPUNIT_ASSERT(default_s->get_mass(1) == 0);
  CPPUNIT_ASSERT(default_s->get_neutral_mass(1) < 0);
  CPPUNIT_ASSERT(default_s->get_singly_charged_mass(1) == 0);

}










