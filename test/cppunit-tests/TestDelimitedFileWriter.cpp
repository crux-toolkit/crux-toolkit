#include <cppunit/config/SourcePrefix.h>
#include "TestDelimitedFileWriter.h"
#include "parameter.h" 

// for enum types
#include "MatchColumns.h"

CPPUNIT_TEST_SUITE_REGISTRATION( TestDelimitedFileWriter );

using namespace std;

void TestDelimitedFileWriter::setUp(){
  initialize_parameters();  // accessing any parameter values requires this

  // initialize variables to use for testing
  filename = "tiny-delim.txt";
  remove(filename);
  defaultWriterPtr = new DelimitedFileWriter(filename);

  // a default header
  colNames.push_back("column#1");
  colNames.push_back("column#2");
  colNames.push_back("column#3");
  header = "column#1	column#2	column#3";

}

void TestDelimitedFileWriter::tearDown(){
  // delete anything you allocated
  delete defaultWriterPtr;
}

void TestDelimitedFileWriter::setColumns(){
  defaultWriterPtr->setColumnNames(colNames);
  defaultWriterPtr->writeHeader();

  defaultWriterPtr->setColumnName("addedcolumn#4", 3); // add column to the end
  defaultWriterPtr->setColumnName("mycol#7", 6);// add past end
  defaultWriterPtr->setColumnName("changedcolumn#2", 1);// change a col
  defaultWriterPtr->writeHeader();

  delete defaultWriterPtr; // to close the file
  defaultWriterPtr = NULL;

  fstream file(filename);
  string line;
  getline(file, line);
  CPPUNIT_ASSERT(line == header);
  
  getline(file, line);
  CPPUNIT_ASSERT(line == 
                 "column#1	changedcolumn#2	column#3	addedcolumn#4	column_5	column_6	mycol#7");

}

// test setting fields in a row and writing them to file
// writing should eliminate the previous row
// with header, row length should be same length as header
void TestDelimitedFileWriter::setIntValues(){
  // set values for file with no header
  defaultWriterPtr->setColumnCurrentRow(0, 777, 0);// col, value, precision
  defaultWriterPtr->setColumnCurrentRow(1, 888, 0);
  defaultWriterPtr->writeRow();

  // set values out of order and with gaps
  defaultWriterPtr->setColumnCurrentRow(0, 777, 0);// col, value, precision
  defaultWriterPtr->setColumnCurrentRow(3, 888, 0);
  defaultWriterPtr->setColumnCurrentRow(1, 999, 0);
  defaultWriterPtr->writeRow();

  // write a header
  defaultWriterPtr->setColumnNames(colNames);
  defaultWriterPtr->writeHeader();
  // set just one of the three cols
  defaultWriterPtr->setColumnCurrentRow(1, 999, 0);
  defaultWriterPtr->writeRow();

  // now look at file
  delete defaultWriterPtr; // to close the file
  defaultWriterPtr = NULL;

  fstream file(filename);
  string line;
  getline(file, line);
  CPPUNIT_ASSERT(line == "777	888");

  getline(file, line);
  CPPUNIT_ASSERT(line == "777	999		888");

  getline(file, line);// header
  getline(file, line);
  CPPUNIT_ASSERT(line == "	999	");

  file.close();
}

// try writing char* and std::string values
// also test that closing a file writes the last row
void TestDelimitedFileWriter::setStringValues(){
  defaultWriterPtr->setColumnCurrentRow(0, "const char*");
  string a = "string";
  defaultWriterPtr->setColumnCurrentRow(1, a);
  defaultWriterPtr->writeRow();

  delete defaultWriterPtr;
  defaultWriterPtr = NULL;

  fstream file(filename);
  string line;
  getline(file, line);
  CPPUNIT_ASSERT(line == "const char*	string");

  file.close();
}

// try writing a float and a double
// use the default precision and a specific precision
void TestDelimitedFileWriter::setFloatValues(){
  float a = 77.77;
  double b = 88.88;
  defaultWriterPtr->setColumnCurrentRow(0, a); // default precision
  defaultWriterPtr->setColumnCurrentRow(1, a, 1); // low precision
  defaultWriterPtr->setColumnCurrentRow(2, a, 3); // high precision
  defaultWriterPtr->setColumnCurrentRow(3, b); // default precision
  defaultWriterPtr->setColumnCurrentRow(4, b, 1); // low precision
  defaultWriterPtr->setColumnCurrentRow(5, b, 3); // high precision
  defaultWriterPtr->writeRow();

  delete defaultWriterPtr;
  defaultWriterPtr = NULL;

  fstream file(filename);
  string line;
  getline(file, line);
  CPPUNIT_ASSERT(line == 
           "77.76999664	77.8	77.770	88.88000000	88.9	88.880");

  file.close();
}




