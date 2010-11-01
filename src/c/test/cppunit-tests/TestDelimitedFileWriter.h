#ifndef CPP_UNIT_TESTDELIMITEDFILEWRITER_H
#define CPP_UNIT_TESTDELIMITEDFILEWRITER_H

#include <cppunit/extensions/HelperMacros.h>
#include <cstdio>
#include "DelimitedFileWriter.h"

class TestDelimitedFileWriter : public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE( TestDelimitedFileWriter );
  CPPUNIT_TEST( setColumns );
  CPPUNIT_TEST( setIntValues );
  CPPUNIT_TEST( setStringValues );
  CPPUNIT_TEST( setFloatValues );
  CPPUNIT_TEST_SUITE_END();
  
 protected:
  // variables to use in testing
  DelimitedFileWriter defaultWriter;
  DelimitedFileWriter* defaultWriterPtr;
  const char* filename;
  std::vector<std::string> colNames;
  std::string header;


 public:
  void setUp();
  void tearDown();

 protected:
  void setColumns();
  void setIntValues();
  void setStringValues();
  void setFloatValues();
};

#endif //CPP_UNIT_TESTDELIMITEDFILEWRITER_H
