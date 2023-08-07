#ifndef CPP_UNIT_TESTMATCHFILEREADER_H
#define CPP_UNIT_TESTMATCHFILEREADER_H

#include <cppunit/extensions/HelperMacros.h>
#include "MatchFileReader.h"

class TestMatchFileReader : public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE( TestMatchFileReader );
  CPPUNIT_TEST( getColumns );
  CPPUNIT_TEST_SUITE_END();
  
 protected:
  // variables to use in testing
  MatchFileReader defaultReader;
  MatchFileReader* tinyReader;//("sample-files/tiny-tab-file.txt");

 public:
  void setUp();
  void tearDown();

 protected:
  void getColumns();
};

#endif //CPP_UNIT_TESTMATCHFILEREADER_H
