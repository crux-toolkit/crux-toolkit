#include <cppunit/config/SourcePrefix.h>
#include "TestMatchFileReader.h"
#include "parameter.h" 

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( TestMatchFileReader );

void TestMatchFileReader::setUp(){
  // initialize_parameters();  // accessing any parameter values requires this

  // initialize variables to use for testing
  tinyReader = new MatchFileReader("sample-files/tiny-tab-file.txt");
}

void TestMatchFileReader::tearDown(){
  // delete anything you allocated
  delete tinyReader;
}

void TestMatchFileReader::getColumns(){
  // try getting columns from a MFR with no file
  vector<bool> columns_present;
  defaultReader.getMatchColumnsPresent(columns_present);
  CPPUNIT_ASSERT(columns_present.empty());

  // open file and get columns
  defaultReader.loadData("sample-files/tiny-tab-file.txt");
  // after reading the file, it should look like this
  vector<bool> correct_values(NUMBER_MATCH_COLUMNS, false);
  for(int i = 0; i < 5; i++){
    correct_values.at(i) = true;
  }

  defaultReader.getMatchColumnsPresent(columns_present);
  CPPUNIT_ASSERT(columns_present == correct_values );

  // test a reader created with a file
  tinyReader->getMatchColumnsPresent(columns_present);
  CPPUNIT_ASSERT(columns_present == correct_values );

  /* in case you want to test them individually
  // the first 5 should be true and remaining false
  for(int i = 0; i < 5; i++){
    CPPUNIT_ASSERT(columns_present.at(i) == true );
  }
  for(int i = 5; i < NUMBER_MATCH_COLUMNS; i++){
    CPPUNIT_ASSERT(columns_present.at(i) == false );
  }
  */
}










