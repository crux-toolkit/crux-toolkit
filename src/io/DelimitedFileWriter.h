/**
 * \file DelimitedFileWriter.h
 * $Revision: 1.0 $
 * DATE: October 19, 2010
 * AUTHOR: Barbara Frewen
 * \brief Object for writing tab-delimited files.
 * This class writes to files one line at a time, creating a table of
 * data.  The values for the current row are set one column at a time
 * and then the row is written to file, separating columns of data
 * with the delimiter.  The default delimiter is tab, but can be set
 * to any character.  A header row may be defined and printed to the
 * file at any row in the flile.  Once a header has been written,
 * every row after that will have the same number of fields.
 */

#ifndef DELIMITED_FILE_WRITER_H
#define DELIMITED_FILE_WRITER_H

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ios>
#include "parameter.h"

class DelimitedFileWriter {

 protected:
  std::ofstream* file_ptr_; ///< the file to write to
  char delimiter_; ///< separate columns with this character
  std::vector<std::string> column_names_; ///< one entry per column
  std::vector<std::string> current_row_; ///< values for next row to write

 public:
  /**
   * \returns An empty DelimitedFileWriter object.
   */
  DelimitedFileWriter();

  /**
   * \returns A DelimitedFileWriter object with the given file to
   * write to.
   */
  DelimitedFileWriter(const char* filename); // full path of file

  /**
   * Destructor
   */
  virtual ~DelimitedFileWriter();

  /**
   * Writes any existing data, closes any open file and opens the
   * given file.
   */
  virtual void openFile(const char* filename);

  /**
   * Sets the delimeter to separate columns.
   */
  void setDelimiter(char delimiter);

  /**
   * Sets the name of the column at the given index, beginning with
   * zero.  Replaces any existing name.
   */
  void setColumnName(const std::string& name, ///< new name to set
                     unsigned int col_idx);   ///< column to name

  /**
   * Sets the names of all columns, clearing any existing names.
   */
  void setColumnNames(const std::vector<std::string>& names);

  /**
   * Writes any given column headers to file.
   */
  void writeHeader();

  /**
   * Writes the data of the current line to file, clears current data.
   */
  void writeRow();

  // template functions only compile correctly when in the header file

  /**
   * Sets the value for the current row of the given column.
   */
  template<typename ValueType>
  void setColumnCurrentRow
    (unsigned int col_idx, ///< set value for this column
     const ValueType& value, ///< the value to set
     unsigned int precision, ///<written at this precision
     bool fixed_float=true){ ///<use fixed float notation?

    // make sure the current_row_ is long enough
    size_t max_size = std::max((size_t)col_idx + 1, column_names_.size());
    
    while( current_row_.size() < max_size ){
      current_row_.push_back("");
    }
    
    current_row_[col_idx] = to_string(value, precision, fixed_float);
  }

  /**
   * Sets the value for the current row of the given column.
   */
  template<typename ValueType>
  void setColumnCurrentRow
    (unsigned int col_idx, ///< set value for this column
     const ValueType& value){ ///< the value to set
    // make sure the current_row_ is long enough
    size_t max_size = std::max((size_t)col_idx + 1, column_names_.size());
    
    while( current_row_.size() < max_size ){
      current_row_.push_back("");
    }
    
    current_row_[col_idx] = to_string(value);
  }

  /**
   * Turn the given value into a string.  For floating point numbers,
   * use the --precision option value.
   */
  template<typename TValue>
  static std::string to_string(TValue value) {
    return to_string(value, get_int_parameter("precision"));
  }

  /**
   * Turn the given value into a string.  For floating point numbers,
   * use the given precision.
   */
  template<typename TValue>
  static std::string to_string
    (TValue& value,
     int precision,
     bool fixed_float = true) {

    std::ostringstream oss;
    oss << std::setprecision(precision);
    if (fixed_float) {
      oss << std::fixed;
    } else {
      oss.unsetf(std::ios_base::floatfield);
    }
    oss << value;
    std::string out_string = oss.str();
    return (out_string != "-0") ? out_string : "0";
  }

};





















#endif //DELIMITED_FILE_WRITER_H
