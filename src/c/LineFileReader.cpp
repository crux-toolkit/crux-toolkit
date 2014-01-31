/*************************************************************************//**
 * \file LineFileReader.cpp
 * \brief Object for parsing the tab-delimited files
 ****************************************************************************/

#include "LineFileReader.h"

#include <fstream>

#include "carp.h"

using namespace std;

/**
 * \returns a LineFileReader object
 */  
LineFileReader::LineFileReader() {

  file_ptr_ = NULL;
}

/**
 * \returns a LineFileReader object and loads the tab-delimited
 * data specified by file_name.
 */  
LineFileReader::LineFileReader(
  const char *file_name ///< the path of the file to read 
  ){

  file_ptr_ = NULL;
  loadData(file_name);
}

/** 
 * \returns a LineFileReader object and loads the tab-delimited
 * data specified by file_name.
 */
LineFileReader::LineFileReader(
    const string& file_name ///< the path of the file  to read
  ){

  file_ptr_ = NULL;
  loadData(file_name);
}

/**
 * Destructor
 */
LineFileReader::~LineFileReader() {

  if (file_ptr_ != NULL) {

    file_ptr_ -> close();
    delete file_ptr_;
  }
}

/**
 * clears the current data and column names,
 * parses the header if it exists,
 * reads the file one line at a time while
 * populating the data matrix with the 
 * strings separated by tabs.
 */
void LineFileReader::loadData(
  const char *file_name ///< the file path
  ) {

  file_name_ = string(file_name);

  file_ptr_ = new fstream(file_name, ios::in);
  current_row_ = 0;


  if (!file_ptr_ -> is_open()) {
    carp(CARP_ERROR, "Opening %s or reading failed", file_name);
    return;
  } else {
    has_next_ = getline(*file_ptr_, next_data_string_) != NULL;
    carp(CARP_DEBUG, "first line:%s",next_data_string_.c_str());
    if (!has_next_) {
      carp(CARP_WARNING,"No data found!");
    } 
  }
}

/**
 * loads a tab delimited file
 */
void LineFileReader::loadData(
  const string& file ///< the file path
  ) {

  loadData(file.c_str());
}

/**
 * \returns the current row string
 */
const string& LineFileReader::current() {
  if (!has_current_) {
    //carp(CARP_FATAL, "End of file!");
  }
  
  return current_data_string_;
}

/*Iterator functions.*/

/**
 * resets the file pointer to the beginning of the file.
 */
void LineFileReader::reset() {

  if (file_ptr_ != NULL) {
    file_ptr_ -> close();
    delete file_ptr_;
  }
  loadData(file_name_);
  
}

/**
 * parses the next line in the file. 
 */
const string& LineFileReader::next() {
  //Do we already have a next line ready?
  if (has_next_) {
    current_row_++;
    current_data_string_ = next_data_string_;
    //read next line
    has_next_ = getline(*file_ptr_, next_data_string_) != NULL;
    has_current_ = true;
  } else {
    has_current_ = false;
  }

  return current_data_string_;
}


int LineFileReader::getCurrentRow() {

  return current_row_;

}

/**
 * \returns whether there are more rows to 
 * iterate through
 */
bool LineFileReader::hasNext() {
  return has_next_;
}


