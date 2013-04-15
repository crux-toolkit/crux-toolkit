#ifndef CPP_UNIT_PROTEIN_H
#define CPP_UNIT_PROTEIN_H

#include <cppunit/extensions/HelperMacros.h>
#include "Protein.h"

class TestProtein : public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE( TestProtein );
  CPPUNIT_TEST( testPeptideShuffle );
  CPPUNIT_TEST_SUITE_END();
  
 protected:
  // variables to use in testing
  const char* tryptic_seq;
  Crux::Protein* prot1;

 public:
  void setUp();
  void tearDown();

 protected:
  void testPeptideShuffle();
};

#endif //CPP_UNIT_PROTEIN_H
