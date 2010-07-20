/*************************************************************************//**
 * \file MatchFileReader.cpp
 * \brief Object for parsing the tab-delimited files
 ****************************************************************************/

#include "MatchColumns.h"
#include "MatchFileReader.h"
#include "DelimitedFile.h"

using namespace std;

/**
 * \returns a blank MatchFileReader object 
 */
MatchFileReader::MatchFileReader() : DelimitedFileReader() {
  ;
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
MatchFileReader::MatchFileReader(const string& file_name) : DelimitedFileReader(file_name, true) {
  parseHeader();
}

/**
 * Destructor
 */
MatchFileReader::~MatchFileReader() {
}

/**
 * parses the header and builds the internal hash table
 */
void MatchFileReader::parseHeader() {
  
  for (int idx=0;idx<NUMBER_MATCH_COLUMNS;idx++) {
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

    carp(CARP_FATAL,
	 "column \"%s\" not found",
	 get_column_header(col_type));
    return -1;
  } else {
    return DelimitedFileReader::getFloat(idx);
  }
}

/**
 * \returns the integer value of a cell, checks for infinity.
 */
int MatchFileReader::getInteger(
  MATCH_COLUMNS_T col_type ///< the column name
) {

  int idx = match_indices_[col_type];
  if (idx == -1) {

    carp(CARP_FATAL,
	 "column \"%s\" not found",
	 get_column_header(col_type));
    return -1;
  } else {

    return DelimitedFileReader::getInteger(idx);
  }
}


//this is to avoid the return warning
string BLANK_STRING="";

/**
 * \returns the string value of a cell
 */
std::string& MatchFileReader::getString(
  MATCH_COLUMNS_T col_type ///<the column type
) {

  int idx = match_indices_[col_type];
  if (idx == -1) {

    carp(CARP_FATAL,
	 "column \"%s\" not found",
	 get_column_header(col_type));
    return BLANK_STRING;
  } else {

    return DelimitedFileReader::getString(idx);
  }
}

bool MatchFileReader::empty(
  MATCH_COLUMNS_T col_type ///<the column type
) {

  int idx = match_indices_[col_type];
  if (idx == -1) {
    carp(CARP_FATAL,
	 "column \"%s\" not found",
	 get_column_header(col_type));
    return true;
  } else {

    return DelimitedFileReader::getString(idx).empty();
  }
}

/**
 * Gets a vector of strings from cell where the
 * string in the cell has delimiters that are
 * different than the column delimiter. The
 * default delimiter is a comma.
 * Uses the current_row_ as the row index.
 * Clears the integer vector before 
 * populating it.
 */
void MatchFileReader::getStringVectorFromCell(
      MATCH_COLUMNS_T col_type, ///< the column name
      std::vector<std::string>& string_vector, ///<the vector of integers
      char delimiter ///<the delimiter to use
) {

  string& string_ans = getString(col_type);

  //get the list of strings separated by delimiter
  string_vector.clear();
  DelimitedFile::tokenize(string_ans, string_vector, delimiter);
}

