#ifndef CPP_UNIT_TESTXML_H
#define CPP_UNIT_TESTXML_H

#include <cppunit/extensions/HelperMacros.h>
#include <map>
#include "crux-utils.h"
#include "Match.h"
#include "Peptide.h"
#include "modifications.h"


/*
 * Test of the methods in match.cpp
 */


class TestXml: public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE( TestXml );
  CPPUNIT_TEST( getNumInternalCleavageNone );
  CPPUNIT_TEST( getNumInternalCleavageTwo );
  CPPUNIT_TEST( getNumTerminalCleavageTwoDash );
  CPPUNIT_TEST( getNumTerminalCleavageOnePrev );
  CPPUNIT_TEST( getNumTerminalCleavageOneNext );
  CPPUNIT_TEST( findVariableModificationsNone );
  CPPUNIT_TEST( findVariableModificationsThree );
  CPPUNIT_TEST( findStaticModificationsNoneFromVariable );
  CPPUNIT_TEST( findStaticModificationsNone );
  CPPUNIT_TEST( findStaticModificationsThree );
  CPPUNIT_TEST_SUITE_END();

 protected:
  std::map<int, double> var_mods;
  std::map<int, double> static_mods;
  double mass_v, mass_p;
  MASS_TYPE_T isotopic_type;
  std::string ord_pep_seq;

 public:
  void setUp();
  void tearDown();

 protected:
  void getNumInternalCleavageNone();
  void getNumInternalCleavageTwo();
  void getNumTerminalCleavageTwoDash();
  void getNumTerminalCleavageOnePrev();
  void getNumTerminalCleavageOneNext();
  void findVariableModificationsNone();
  void findVariableModificationsThree();
  void findStaticModificationsThree();
  void findStaticModificationsNone();
  void findStaticModificationsNoneFromVariable();
  
};

#endif //CPP_UNIT_TESTXML_H
