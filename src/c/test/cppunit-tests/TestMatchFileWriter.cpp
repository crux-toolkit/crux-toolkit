#include <cppunit/config/SourcePrefix.h>
#include "TestMatchFileWriter.h"
#include "parameter.h" 
#include "MatchSearch.h"
#include "SequestSearch.h"
//#include "Percolator.h"
#include "QRanker.h"
#include "ComputeQValues.h"

#include <fstream>
#include <iostream>

// from parameter.cpp
bool reset_parameter(const char* name, const char* value);

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION( TestMatchFileWriter );

void TestMatchFileWriter::setUp(){
  initialize_parameters();  // accessing any parameter values requires this
  reset_parameter("overwrite", "TRUE");

  // initialize variables to use for testing
  filename = "tiny-match.txt";
  defaultWriterPtr = new MatchFileWriter(filename);
}

void TestMatchFileWriter::tearDown(){
  delete defaultWriterPtr;
}

// todo add a bool skipHeader = false arg
void testFileLine(const char* fileName, const char* correctLine){
  fstream file(fileName);
  string line;
  getline(file, line);
  //cerr << "line is '" << line << "'" << endl;
  //cerr << "should be '" << correctLine << "'" << endl;
  CPPUNIT_ASSERT(line == correctLine);
  file.close();
}

// try setting column names, by match_column_t and printing header
void TestMatchFileWriter::setSingleColumns(){
  // order printed shouldn't depend on order set
  defaultWriterPtr->addColumnName(CHARGE_COL);
  defaultWriterPtr->addColumnName(SCAN_COL);
  defaultWriterPtr->writeHeader();

  delete defaultWriterPtr;
  defaultWriterPtr = NULL;

  testFileLine(filename, "scan	charge");
}


// try setting column names by bool vector and printing header
void TestMatchFileWriter::setColumns(){
  vector<bool> useThese(5, false);
  useThese[1] = true; // charge
  useThese[0] = true; // scan
  useThese[4] = true; // peptide mass

  defaultWriterPtr->addColumnNames(useThese);
  defaultWriterPtr->writeHeader();

  delete defaultWriterPtr;
  defaultWriterPtr = NULL;

  testFileLine(filename, "scan	charge	peptide mass");
}

// try setting column names by both methods
void TestMatchFileWriter::setColumnsBothMethods(){
  vector<bool> useThese(5, false);
  useThese[1] = true; // charge
  useThese[0] = true; // scan
  useThese[4] = true; // peptide mass

  defaultWriterPtr->addColumnName(DELTA_CN_COL); // this should not be first
  defaultWriterPtr->addColumnNames(useThese);
  defaultWriterPtr->addColumnName(SCAN_COL); // this shouldn't change anything
  // this should come before peptide mass even though it's being set after
  defaultWriterPtr->addColumnName(SPECTRUM_NEUTRAL_MASS_COL);
  defaultWriterPtr->writeHeader();

  delete defaultWriterPtr;
  defaultWriterPtr = NULL;

  testFileLine(filename, "scan	charge	spectrum neutral mass	peptide mass	delta_cn");
}

// set some values by MATCH_COLUMN
void TestMatchFileWriter::setValues(){
  // define which columns to use
  defaultWriterPtr->addColumnName(CHARGE_COL);
  defaultWriterPtr->addColumnName(SCAN_COL);
  defaultWriterPtr->writeHeader();

  // set some values
  defaultWriterPtr->setColumnCurrentRow(CHARGE_COL, 3);
  defaultWriterPtr->setColumnCurrentRow(SCAN_COL, 10);
  defaultWriterPtr->setColumnCurrentRow(DELTA_CN_COL, 10);//shouldn't be printed

  defaultWriterPtr->writeRow();

  delete defaultWriterPtr;
  defaultWriterPtr = NULL;

  fstream file(filename);
  string line;
  getline(file, line);// header
  getline(file, line);
  CPPUNIT_ASSERT(line == "10	3");
}

// write a file with all values and confirm that precision is correct
void TestMatchFileWriter::precision(){
  double value = 777.123456789;

  for(int i=0; i < NUMBER_MATCH_COLUMNS; i++){
    defaultWriterPtr->addColumnName((MATCH_COLUMNS_T)i);
  }
  defaultWriterPtr->writeHeader();

  // set all values as the same
  for(int i=0; i < NUMBER_MATCH_COLUMNS; i++){
    defaultWriterPtr->setColumnCurrentRow((MATCH_COLUMNS_T)i, value);
  }
  defaultWriterPtr->writeRow();

  delete defaultWriterPtr;
  defaultWriterPtr = NULL;

  fstream file(filename);
  string line;
  getline(file, line);// header
  //cerr<< line << endl;
  getline(file, line);
  //cerr << line << endl;
  CPPUNIT_ASSERT(line ==
                 "777	777	777.1235	777.1235	777.1235	777.12346	777.12346	777	777.12346	777	777.12346	777.12346	777.12346	777.12346	777.12346	777.12346	777	777.12346	777.12346	777.12346	777.12346	777.12346	777	777	777	777	777	777	777	777	777.123	777.123	777.123	777.123	777	777.12346	777.12346	777.12346	777.12346	777	777");


}

// try setting the columns for each command type
void TestMatchFileWriter::setColumnsCommand(){

  // search-for-matches target
  defaultApplicationPtr = new MatchSearch();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false);
  defaultWriterPtr->writeHeader();
  delete defaultApplicationPtr;
  delete defaultWriterPtr;
  
  defaultApplicationPtr = NULL;
  defaultWriterPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	matches/spectrum	sequence	cleavage type	protein id	flanking aa");

  // search-for-matches decoy
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new MatchSearch();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, true);
  defaultWriterPtr->writeHeader();
  delete defaultApplicationPtr;
  delete defaultWriterPtr;
  defaultApplicationPtr = NULL;
  defaultWriterPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	matches/spectrum	sequence	cleavage type	protein id	flanking aa	unshuffled sequence	decoy matches/spectrum");

  // search-for-matches target, p-values
  reset_parameter("compute-p-values", "TRUE");
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new MatchSearch();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	p-value	matches/spectrum	sequence	cleavage type	protein id	flanking aa	eta	beta	shift	corr");

  // search-for-matches target, sp
  reset_parameter("compute-p-values", "FALSE");// undo previous
  reset_parameter("compute-sp", "TRUE");
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new MatchSearch();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	sp score	sp rank	xcorr score	xcorr rank	b/y ions matched	b/y ions total	matches/spectrum	sequence	cleavage type	protein id	flanking aa");

  // sequest-search
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new SequestSearch();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	sp score	sp rank	xcorr score	xcorr rank	b/y ions matched	b/y ions total	matches/spectrum	sequence	cleavage type	protein id	flanking aa");

  // sequest-search no sp
  reset_parameter("max-rank-preliminary", "0");
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new SequestSearch();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	matches/spectrum	sequence	cleavage type	protein id	flanking aa");

  // for post-search commands, create a bool vector 
  vector<bool> printThese(NUMBER_MATCH_COLUMNS, false);
  printThese[SCAN_COL] = true;
  printThese[CHARGE_COL] = true;
  printThese[SPECTRUM_PRECURSOR_MZ_COL] = true;
  printThese[SPECTRUM_NEUTRAL_MASS_COL] = true;
  printThese[PEPTIDE_MASS_COL] = true;
  printThese[DELTA_CN_COL] = true;
  printThese[XCORR_SCORE_COL] = true;
  printThese[XCORR_RANK_COL] = true;
  printThese[MATCHES_SPECTRUM_COL] = true;
  printThese[SEQUENCE_COL] = true;
  printThese[CLEAVAGE_TYPE_COL] = true;
  printThese[PROTEIN_ID_COL] = true;
  printThese[FLANKING_AA_COL] = true;
  //  printThese[UNSHUFFLED_SEQUENCE_COL] = true;

  /* percolator
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new Percolator();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false, printThese);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	percolator score	percolator rank	percolator q-value	percolator PEP	matches/spectrum	sequence	cleavage type	protein id	flanking aa");
*/
  // q-ranker
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new QRanker();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false, printThese);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	q-ranker score	q-ranker q-value	q-ranker PEP	matches/spectrum	sequence	cleavage type	protein id	flanking aa");

  // q-values, decoys
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new ComputeQValues();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false, printThese);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	decoy q-value (xcorr)	decoy PEP (xcorr)	matches/spectrum	sequence	cleavage type	protein id	flanking aa");

  // q-values, weibull
  printThese[PVALUE_COL] = true;
  defaultWriterPtr = new MatchFileWriter(filename);
  defaultApplicationPtr = new ComputeQValues();
  defaultWriterPtr->addColumnNames(defaultApplicationPtr, false, printThese);
  defaultWriterPtr->writeHeader();
  delete defaultWriterPtr;
  delete defaultApplicationPtr;
  defaultWriterPtr = NULL;
  defaultApplicationPtr = NULL;
  testFileLine(filename, "scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	xcorr score	xcorr rank	p-value	Weibull est. q-value	Weibull est. PEP	matches/spectrum	sequence	cleavage type	protein id	flanking aa");
} 


