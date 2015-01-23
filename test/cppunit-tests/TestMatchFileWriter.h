#ifndef CPP_UNIT_MatchFileWriter_H
#define CPP_UNIT_MatchFileWriter_H

#include <cppunit/extensions/HelperMacros.h>
#include "MatchFileWriter.h"

class TestMatchFileWriter : public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE( TestMatchFileWriter );
  CPPUNIT_TEST( setSingleColumns );
  CPPUNIT_TEST( setColumns );
  CPPUNIT_TEST( setColumnsBothMethods );
  CPPUNIT_TEST( setValues );
  CPPUNIT_TEST( precision );
  CPPUNIT_TEST( setColumnsCommand );
  CPPUNIT_TEST_SUITE_END();
  
 protected:
  MatchFileWriter* defaultWriterPtr;
  CruxApplication* defaultApplicationPtr;
  const char* filename;

 public:
  void setUp();
  void tearDown();

 protected:
  void setSingleColumns();
  void setColumns();
  void setColumnsBothMethods();
  void setValues();
  void precision();
  void setColumnsCommand();
};

#endif //CPP_UNIT_MatchFileWriter_H
