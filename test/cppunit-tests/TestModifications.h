#ifndef CPP_UNIT_TESTMODIFICATIONS_H
#define CPP_UNIT_TESTMODIFICATIONS_H

#include <cppunit/extensions/HelperMacros.h>
#include "crux-utils.h"
#include "modifications.h"

/*
 * Test of the methods in modifications.cpp.  
 */

class TestModifications : public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE( TestModifications );
  CPPUNIT_TEST( convertToModAaSeq );
  CPPUNIT_TEST_SUITE_END();
  
 protected:
  const char* seq1;
  const char* seq2;
  MODIFIED_AA_T* modSeq1;
  MODIFIED_AA_T* modSeq2;
  int len2;
  AA_MOD_T *amod1, *amod2, *amod3;
  AA_MOD_T* amod_list[3];


 public:
  void setUp();
  void tearDown();

 protected:
  void convertToModAaSeq();
};

#endif //CPP_UNIT_TESTMODIFICATIONS_H
