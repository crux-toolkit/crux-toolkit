/*************************************************************************//**
 * \file DelimitedFileReader.cpp
 * \brief Object for parsing the tab-delimited files
 ****************************************************************************/

#include "DelimitedFileReader.h"

#include <fstream>

#include <iostream>
#include <string>

#include "carp.h"
#include "DelimitedFile.h"

using namespace std;

/**
 * \returns a DelimitedFileReader object
 */  
DelimitedFileReader::DelimitedFileReader() {
  num_rows_valid_ = false;
  istream_ptr_ = NULL;
  delimiter_ = '\t';
  owns_stream_ = false;
}

/**
 * \returns a DelimitedFileReader object and loads the tab-delimited
 * data specified by file_name.
 */ 
DelimitedFileReader::DelimitedFileReader(
  const char *file_name, ///< the path of the file to read
  bool has_header, ///< indicates whether the header exists (default true).
  char delimiter ///< the delimiter to use (default tab).
  ) {

  istream_ptr_ = NULL;
  num_rows_valid_ = false;
  delimiter_ = delimiter;
  loadData(file_name, has_header);
}

/** 
 * \returns a DelimitedFileReader object and loads the tab-delimited
 * data specified by file_name.
 */
DelimitedFileReader::DelimitedFileReader(
  const std::string& file_name, ///< the path of the file  to read
  bool has_header, ///< indicates whether the header exists (default true).
  char delimiter ///< the delimiter to use (default tab)
  ) {

  istream_ptr_ = NULL;
  delimiter_ = delimiter;
  loadData(file_name, has_header);
}

/**
 * \returns a DelimitedFileReader object and loads the tab-delimted
 * data specified by istream.
 */
DelimitedFileReader::DelimitedFileReader(
  std::istream* istream_ptr, ///< the stream to be read
  bool has_header, ///<indicates whether header exists
  char delimiter ///< the delimiter to use (default tab)
  ) {

  
  istream_ptr_ = istream_ptr;
  istream_begin_ = istream_ptr->tellg();

  delimiter_ = delimiter;
  has_header_ = has_header;
  owns_stream_ = false;

  loadData();

}


/**
 * Destructor
 */
DelimitedFileReader::~DelimitedFileReader() {

  if (istream_ptr_ != NULL && owns_stream_) {
    delete istream_ptr_;
  }
}

/**
 * \returns the number of rows, assuming a square matrix
 */
unsigned int DelimitedFileReader::numRows() {

  if (!num_rows_valid_) {

    num_rows_ = 0;


    streampos last_pos = istream_ptr_->tellg();

    istream_ptr_->clear();
    istream_ptr_->seekg(istream_begin_, ios::beg);
    
    string temp_str;
    bool has_next = getline(*istream_ptr_,temp_str) != NULL;

    while (has_next) {
      num_rows_++;
      has_next = getline(*istream_ptr_, temp_str) != NULL;
    }
    
    if (has_header_) {
      num_rows_--;
    }
    num_rows_valid_ = true;
    istream_ptr_->clear();
    istream_ptr_->seekg(last_pos);

  }
  return num_rows_;
}
 
/**
 *\returns the number of columns
 */
unsigned int DelimitedFileReader::numCols() {

  return column_names_.size();
}

/*
 *\returns a printable string of the columns available in this file
 */
string DelimitedFileReader::getAvailableColumnsString() {

  ostringstream oss;
  oss << "Available columns:"<<endl;
  for (unsigned int col_idx=0;col_idx<numCols();col_idx++) {
    oss << col_idx << "  " << getColumnName(col_idx) << endl;
  }

  string ans = oss.str();
  return ans;
}

/*
 *\returns the column name header string
 */
string DelimitedFileReader::getHeaderString() {

  if (numCols() == 0) {
    return string("");
  }
  
  ostringstream oss;
  oss << getColumnName(0);
  for (unsigned int col_idx=1;col_idx<numCols();col_idx++) {
    oss << delimiter_ << getColumnName(col_idx);
  }
  
  string ans = oss.str();
  return ans;
}


void DelimitedFileReader::loadData() {

  if (!istream_ptr_->good()) {
    carp(CARP_ERROR, "Stream is not good!");
    carp(CARP_ERROR, "Filename:%s", file_name_.c_str());
    carp(CARP_ERROR, "EOF:%i", istream_ptr_ -> eof());
    carp(CARP_ERROR, "Fail:%i", istream_ptr_ -> fail());
    carp(CARP_ERROR, "Bad:%i", istream_ptr_ -> bad());
    carp(CARP_FATAL, "Exiting....");
  }
  current_row_ = 0;
  num_rows_valid_ = false;
  has_current_ = false;
  column_names_.clear();
  column_mismatch_warned_ = false;
  istream_begin_ = istream_ptr_->tellg(); 

  has_next_ = getline(*istream_ptr_, next_data_string_) != NULL;
  
  if (has_header_) {
    if (has_next_) {
      tokenize(next_data_string_, column_names_, delimiter_);
      has_next_ = getline(*istream_ptr_, next_data_string_) != NULL;
    }
    else {
      carp(CARP_WARNING,"No data/headers found!");
      return;
    }
  }

  if (has_next_) {
    next();
  } 

}

/**
 * clears the current data and column names,
 * parses the header if it exists,
 * reads the file one line at a time while
 * populating the data matrix with the 
 * strings separated by tabs.
 */
void DelimitedFileReader::loadData(
  const char *file_name, ///< the file path
  bool has_header ///< header indicator
  ) {

  file_name_ = string(file_name);
  has_header_ = has_header;

  //special case, if filename is '-', then use standard input.
  if (file_name_ == "-") {
    istream_ptr_ = &cin;
    owns_stream_ = false;
  } else {
    istream_ptr_ = new ifstream(file_name, ios::in);
    owns_stream_ = true;
  }
  loadData();

}

/**
 * loads a tab delimited file
 */
void DelimitedFileReader::loadData(
  const string& file, ///< the file path
  bool has_header ///< header indicator
  ) {

  loadData(file.c_str(), has_header);
}

/**
 * finds the index of a column
 *\returns the column index, -1 if not found.
 */ 
int DelimitedFileReader::findColumn(
  const string& column_name ///< the column name
  ) {

  for (unsigned int col_idx=0;col_idx<column_names_.size();col_idx++) {
    if (column_names_[col_idx] == column_name) {
      return col_idx;
    }
  }
  return -1;
}

/**
 * finds the index of a column
 *\returns the column index, -1 if not found.
 */ 
int DelimitedFileReader::findColumn(
  const char* column_name ///< the column name
) {
  string sname = string(column_name);
  return findColumn(sname);
}

/**
 *\returns the name of the column
 */
const string& DelimitedFileReader::getColumnName(
  unsigned int col_idx ///< the column index
  ) {
  return column_names_.at(col_idx);
}

/**
 *\returns the column_names
 */
const vector<string>& DelimitedFileReader::getColumnNames() {

  return column_names_;
}

/**
 *\returns the current row index
 */
int DelimitedFileReader::getCurrentRowIndex() const {
  return current_row_;
}


/**
 * \returns the current row string
 */
const string& DelimitedFileReader::getString() {
  if (!has_current_) {
    carp(CARP_FATAL, "End of file!");
  }
  
  return current_data_string_;
}

/**
 *\returns the string value of the cell
 */
const string& DelimitedFileReader::getString(
  unsigned int col_idx ///< the column index
  ) {

  if (col_idx >= data_.size()) {
    carp(CARP_FATAL, 
      "col idx:%i is out of bounds! (0,%i,%i)",
      col_idx, 
      (column_names_.size()-1),
      (data_.size()-1));
  }

  return data_.at(col_idx);
}

/** 
 * \returns the string value of the cell.
 */
const string& DelimitedFileReader::getString(
  const char* column_name ///<the column name
  ) {
  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s\n" 
                     "Available Columns:%s\n",
                      column_name, getAvailableColumnsString().c_str());

  }
  return getString(col_idx);
}

/**
 * \returns the value of the cell
 * using the current row
 */ 
template<typename TValue>
TValue DelimitedFileReader::getValue(
  unsigned int col_idx ///< the column index 
  ) {

  const string& string_ans = getString(col_idx);
  TValue type_ans;
  from_string<TValue>(type_ans, string_ans);
  return type_ans;
}

/**
 * \returns the FLOAT_T value of a cell, checks for infinity
 */
FLOAT_T DelimitedFileReader::getFloat(
  unsigned int col_idx ///< the column index
  ) {
  
  const string& string_ans = getString(col_idx);
  if (string_ans == "Inf") {

    return numeric_limits<FLOAT_T>::infinity();
  } else if (string_ans == "-Inf") {

    return -numeric_limits<FLOAT_T>::infinity();
  }
  else {
    return getValue<FLOAT_T>(col_idx);
  }
}

/** 
 * \returns the FLOAT_T value of a cell, checks for infinity.
 */
FLOAT_T DelimitedFileReader::getFloat(
    const char* column_name ///<the column name
) {

  carp(CARP_DETAILED_DEBUG, "getFloat for %s", column_name);
  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s\n" 
                     "Available Columns:%s\n",
                     column_name, getAvailableColumnsString().c_str());
  }
  return getFloat(col_idx);
}

/**
 * \returns the double value of a cell, checks for infinity. 
 */
double DelimitedFileReader::getDouble(
  unsigned int col_idx ///< the column index 
  ) {

  const string& string_ans = getString(col_idx);
  if (string_ans == "") {

    return 0.0;
  } else if (string_ans == "Inf") {

    return numeric_limits<double>::infinity();
  } else if (string_ans == "-Inf") {

    return -numeric_limits<double>::infinity();
  }
  else {

    return getValue<double>(col_idx);
  }
}

/** 
 * \returns the double value of a cell, checks for infinity.
 */
double DelimitedFileReader::getDouble(
  const char* column_name ///<the column name
  ) {

  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s\n" 
                     "Available Columns:%s\n",
                     column_name, getAvailableColumnsString().c_str());
  }
  return getDouble(col_idx);
}

/**
 * \returns the integer value of a cell. 
 */
int DelimitedFileReader::getInteger(
  unsigned int col_idx ///< the column index 
  ) {
  //TODO : check the string for a valid integer.
  return getValue<int>(col_idx);
}

/**
 * \returns the integer value of a cell, checks for infintiy.
 */
int DelimitedFileReader::getInteger(
  const char* column_name ///< the column name
) {

  int col_idx = findColumn(column_name);
  if (col_idx == -1) {
    carp(CARP_FATAL, "Cannot find column %s", column_name);
  }

  return getInteger(col_idx);
}

/**
 * gets an vector of strings from cell where the
 * string in the cell has delimiters that are
 * different than the column delimiter. The
 * default delimiter is a comma
 * uses the current_row_ as the row index.
 * clears the integer vector before 
 * populating it.
 */
void DelimitedFileReader::getStringVectorFromCell(
  const char* column_name, ///< the column name
  std::vector<std::string>& string_vector, ///<the vector of integers
  char delimiter ///<the delimiter to use
  ) {

  const string& string_ans = getString(column_name);

  //get the list of strings separated by delimiter
  string_vector.clear();
  tokenize(string_ans, string_vector, delimiter);
}

/**
 * gets an vector of integers from cell where the
 * string in the cell are integers which are separated
 * by a delimiter which is differnt than the column
 * delimiter.  The default delimiter is a comma
 * uses the current_row_ as the row index.
 * clears the integer vector before 
 * populating it.
 */
void DelimitedFileReader::getIntegerVectorFromCell(
    const char* column_name, ///< the column name
    vector<int>& int_vector, ///<the vector of integers
    char delimiter ///<the delimiter to use
  ) {
  
  //get the list of strings separated by delimiter
  vector<string> string_vector_ans;

  getStringVectorFromCell(column_name, string_vector_ans, delimiter);

  //convert each string into an integer.
  int_vector.clear();

  for (vector<string>::iterator string_iter = string_vector_ans.begin();
    string_iter != string_vector_ans.end();
    ++string_iter) {

    int int_ans;
    from_string<int>(int_ans, *string_iter);
    int_vector.push_back(int_ans);
  }
}

/**
 * gets an vector of doubles from cell where the
 * string in the cell are integers which are separated
 * by a delimiter which is differnt than the column
 * delimiter.  The default delimiter is a comma
 * uses the current_row_ as the row index.
 * clears the double vector before 
 * populating it.
 */
void DelimitedFileReader::getDoubleVectorFromCell(
  const char* column_name, ///< the column name
  std::vector<double>& double_vector, ///<the vector of integers
  char delimiter ///<the delimiter to use
) {
  
  //get the list of strings separated by delimiter
  vector<string> string_vector_ans;

  getStringVectorFromCell(column_name, string_vector_ans, delimiter);

  //convert each string into an integer.
  double_vector.clear();

  for (vector<string>::iterator string_iter = string_vector_ans.begin();
    string_iter != string_vector_ans.end();
    ++string_iter) {

    double double_ans;
    from_string<double>(double_ans, *string_iter);
    double_vector.push_back(double_ans);
  }
}

/*Iterator functions.*/
/**
 * resets the file pointer to the beginning of the file.
 */
void DelimitedFileReader::reset() {
  istream_ptr_->clear();
  istream_ptr_->seekg(istream_begin_, ios::beg);
  loadData();
  
}

/**
 * parses the next line in the file. 
 */
void DelimitedFileReader::next() {

  if (has_next_) {
    current_row_++;
    current_data_string_ = next_data_string_;
    //parse next_data_string_ into data_
    tokenize(current_data_string_, data_, delimiter_);
    //make sure data has the right number of columns for the header.
    if (data_.size() < column_names_.size()) {
      if (!column_mismatch_warned_) {
        carp(CARP_WARNING,
          "Column count %d for line %d is less than header %d",
          data_.size(), current_row_, column_names_.size());
        carp(CARP_WARNING,
          "%s", current_data_string_.c_str());
        carp(CARP_WARNING,
          "Suppressing warnings, other mismatches may exist!");
        column_mismatch_warned_ = true;
      }
      while (data_.size() < column_names_.size()) {
        data_.push_back(string(""));
      }
    }

    //read next line
    has_next_ = getline(*istream_ptr_, next_data_string_) != NULL;
    has_current_ = true;
  } else {
    has_current_ = false;
  }
}

/**
 * \returns whether there are more rows to 
 * iterate through
 */
bool DelimitedFileReader::hasNext() {
  return has_next_ || has_current_;
}
