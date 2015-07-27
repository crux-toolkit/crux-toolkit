/*************************************************************************//**
 * \file MatchFileReader.cpp
 * \brief Object for parsing the tab-delimited files
 ****************************************************************************/

#include "MatchColumns.h"
#include "MatchFileReader.h"
#include "DelimitedFile.h"

#include "model/MatchCollection.h"

using namespace std;

/**
 * \returns a blank MatchFileReader object 
 */
MatchFileReader::MatchFileReader() : DelimitedFileReader(), PSMReader() {
}

/**
 * \returns a MatchFileReader object and loads the tab-delimited
 * data specified by file_name.
 */  
MatchFileReader::MatchFileReader(const char* file_name) : DelimitedFileReader(file_name, true) {
  parseHeader();
}

/** 
 * \returns a MatchFileReader object and loads the tab-delimited
 * data specified by file_name.
 */
MatchFileReader::MatchFileReader(const string& file_name)
  : DelimitedFileReader(file_name, true), PSMReader(file_name) {
  parseHeader();
}

MatchFileReader::MatchFileReader(const string& file_name, Database* database, Database* decoy_database)
  : DelimitedFileReader(file_name, true), PSMReader(file_name, database, decoy_database) {
  parseHeader();
}

MatchFileReader::MatchFileReader(istream* iptr) : DelimitedFileReader(iptr, true, '\t') {
  parseHeader();
}

/**
 * Destructor
 */
MatchFileReader::~MatchFileReader() {
}

/**
 * Open a new file from an existing MatchFileReader.
 */
void MatchFileReader::loadData(
  const char* file_name, ///< new file to open
  bool hasHeader
){
  DelimitedFileReader::loadData(file_name, hasHeader);
  if( hasHeader ){
    parseHeader();
  }
}

/**
 * Open a new file from an existing MatchFileReader.
 */
void MatchFileReader::loadData(
  const string& file_name, ///< new file to open
  bool hasHeader
){
  DelimitedFileReader::loadData(file_name, hasHeader);
  if( hasHeader ){
    parseHeader();
  }
}

/**
 * parses the header and builds the internal hash table
 */
void MatchFileReader::parseHeader() {
  for (int idx = 0; idx < NUMBER_MATCH_COLUMNS; idx++) {
    match_indices_[idx] = findColumn(get_column_header(idx));
  }
}


/**
 * \returns the FLOAT_T value of a cell, checks for infinity
 */
FLOAT_T MatchFileReader::getFloat(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getFloat", get_column_header(col_type));
    return -1;
  }
  return DelimitedFileReader::getFloat(idx);
}

/**
 * \returns the double value of a cell, checks for infinity
 */
double MatchFileReader::getDouble(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  carp(CARP_DETAILED_DEBUG, "reading double from column %s", get_column_header(col_type));
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getDouble", get_column_header(col_type));
    return -1;
  }
  return DelimitedFileReader::getDouble(idx);
}


/**
 * \returns the integer value of a cell, checks for infinity.
 */
int MatchFileReader::getInteger(
  MATCH_COLUMNS_T col_type ///< the column name
) {
  carp(CARP_DETAILED_DEBUG, "Reading integer from column %s", get_column_header(col_type));
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getInteger", get_column_header(col_type));
    return -1;
  }
  return DelimitedFileReader::getInteger(idx);
}

/**
 * \returns the string value of a cell
 */
string MatchFileReader::getString(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  carp(CARP_DETAILED_DEBUG, "Getting string from column %s",get_column_header(col_type));
  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_DEBUG, "column \"%s\" not found for getString", get_column_header(col_type));
    return "";
  }
  return DelimitedFileReader::getString(idx);
}

bool MatchFileReader::empty(
  MATCH_COLUMNS_T col_type ///<the column type
) {
  int idx = match_indices_[col_type];
  if (idx == -1) {
    return true;
  }
  return DelimitedFileReader::getString(idx).empty();
}

/**
 * Fills in the given vector with a bool value indicating if each
 * MATCH_COLUMN_T type is present in the file being read.
 * \returns Argument vector has NUM_MATCH_COLUMN_T values if a
 * valid file is open and header has been parsed, else vector is empty.
 */
void MatchFileReader::getMatchColumnsPresent(
  std::vector<bool>& col_is_present)
{
  col_is_present.clear();

  // has a header been parsed? 
  if( column_names_.empty() ){
    return;
  }
  col_is_present.assign(NUMBER_MATCH_COLUMNS, false);

  for(int col_idx = 0; col_idx < NUMBER_MATCH_COLUMNS; col_idx++){
    col_is_present[col_idx] = (match_indices_[col_idx] > -1);
  }
}

MatchCollection* MatchFileReader::parse() {

  MatchFileReader delimited_result_file(file_path_);
  MatchCollection* match_collection = new MatchCollection();
  match_collection->preparePostProcess();

  match_collection->extendTabDelimited(database_, delimited_result_file, decoy_database_);

  return match_collection;


}

MatchCollection* MatchFileReader::parse(
  const string& file_path,
  Database* database,
  Database* decoy_database
  ) {
  MatchFileReader delimited_result_file(file_path);
  MatchCollection* match_collection = new MatchCollection();
  match_collection->preparePostProcess();
  match_collection->extendTabDelimited(database, delimited_result_file, decoy_database);
  return match_collection;
}

