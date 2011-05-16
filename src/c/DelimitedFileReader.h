/**
 * \file DelimitedFileReader.h
 * $Revision: 1.00 $ 
 * DATE: Jan 7, 2010
 * AUTHOR: Sean McIlwain
 * \brief Object for reading tab-delimited files.
 * This class generates a table of values. The default delimiter is tab.
 * This class is capable of reading string, integers, and floating point
 * Types from each cell of the table.  This class also provides function
 * for reading a list of integers or string from a cell using a delimiter
 * that is different from the column delimiter (default is comma ',').
 * This class reads the data in line by line
 ****************************************************************************/
#ifndef DELIMITEDFILEREADER_H
#define DELIMITEDFILEREADER_H

#include <limits>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "parameter.h"


class DelimitedFileReader {

 protected:
  /** The data vector is indexed by column, then by row
   *  Ex. data_[0][4] means 0th column, 4th row
   *  data_[0] is a string vector of all data in a column.
   */

  std::string next_data_string_; ///<the next data string.
  std::string current_data_string_; ///<the current data string.
  std::vector<std::string> data_; ///<the current vectorized data.
  std::vector<std::string> column_names_; ///<the column names.

  char delimiter_; ///<the delimiter to use.

  unsigned int current_row_; ///<current row count
  bool has_next_; ///<indicator of whether there is another row to parse
  bool has_current_; ///<indicator of whether the current row is parsed
  bool has_header_; ///<indicator of whether there is a header in the file

  bool owns_stream_; ///<indicator of whether the object owns the stream

  std::istream* istream_ptr_; ///<pointer to the stream itself

  std::streampos istream_begin_; ///<position pointer for the beginning of the stream

  std::string file_name_; ///<file name that the stream is open on.

  bool num_rows_valid_; ///<indicator whether the number of rows is valid
  unsigned int num_rows_; ///<number of rows in the file.

  bool column_mismatch_warned_; ///<indicator of whether the column mismatch warning has been issued

  /**
   * clears the current data and column names,
   * parses the header if it exists,
   * reads the file one line at a time while
   * populating the data matrix with the 
   * strings separated by delimiter.
   */
  void loadData();

  virtual void loadData(
    const char *file_name, ///< the file path
    bool has_header=true ///< header indicator
  );

  virtual void loadData(
    const std::string& file_name, ///< the file path
    bool has_header=true ///< header indicator
  );

 public:
  /**
   * \returns a blank DelimitedFileReader object 
   */
  DelimitedFileReader();
  
  /**
   * \returns a DelimitedFileReader object and loads the tab-delimited
   * data specified by file_name.
   */  
  DelimitedFileReader(
    const char *file_name, ///< the path of the file to read
    bool has_header=true, ///< indicates whether the header exists (default true).
    char delimiter='\t' ///< the delimiter to use (default tab).
  );

  /** 
   * \returns a DelimitedFileReader object and loads the tab-delimited
   * data specified by file_name.
   */
  DelimitedFileReader(
    const std::string& file_name, ///< the path of the file  to read
    bool has_header=true, ///< indicates whether the header exists (default true).
    char delimiter='\t' ///< the delimiter to use (default tab)
  );

  DelimitedFileReader(
    std::istream* istream_ptr, ///< the stream to be read
    bool has_header=true, ///<indicates whether header exists
    char delimiter='\t' ///< the delimiter to use (default tab)
  );

  /**
   * Destructor
   */
  virtual ~DelimitedFileReader();

  /**
   *\returns the number of rows, assuming a square matrix
   */
  unsigned int numRows();

  /**
   *\returns the number of columns
   */
  unsigned int numCols();

  /*
   *\returns a printable string of the columns available in this file
   */
  std::string getAvailableColumnsString();

  /*
   *\returns the column name header string
   */
  std::string getHeaderString();


  /**
   * finds the index of a column
   *\returns the column index, -1 if not found.
   */ 
  int findColumn(
    const std::string& column_name ///< the column name
  );

 /**
   * finds the index of a column
   *\returns the column index, -1 if not found.
   */ 
  int findColumn(
    const char* column_name ///< the column name
  );

  /**
   *\returns the name of the column
   */
  const std::string& getColumnName(
    unsigned int col_idx ///< the column index
  );
  
  /**
   *\returns the column_names
   */
  const std::vector<std::string>& getColumnNames();

  /**
   *\returns the current row index
   */
  int getCurrentRowIndex() const;

  /**
   * \returns the current row string
   */
  const std::string& getString();

  /**
   * \returns the string value of the cell
   * using the current row
   */
  const std::string& getString(
    const char* column_name ///<the column name
  );

  /**
   * \returns the string value of the cell
   * using the current row
   */
  const std::string& getString(
    unsigned int col_idx ///< the column index
  );

  /**
   * \returns the value of the cell
   * using the current row
   */ 
  template<typename TValue>
  TValue getValue(
    unsigned int col_idx ///< the column index
  );

  /**
   * \returns the double value of the cell, checks for infinity
   * uses the current row
   */
  FLOAT_T getFloat(
    const char* column_name ///<the column name
  );

  FLOAT_T getFloat(
    unsigned int col_idx ///<the col index
  );

  /**
   * \returns the double value from cell, checks for infinity
   * uses the current_row_ as the row index
   */
  double getDouble(
    const char* column_name ///<the column name
  );

  double getDouble(
    unsigned int col_idx ///<the col index
  );

  /**
   * get an integer type from cell, checks for infinity.
   * uses the current_row_ as the row index.
   */
  int getInteger(
    const char* column_name ///< the column name
  );

  int getInteger(
    unsigned int col_idx ///<the col index
  );

  /**
   * gets an vector of strings from cell where the
   * string in the cell has delimiters that are
   * different than the column delimiter. The
   * default delimiter is a comma
   * uses the current_row_ as the row index.
   * clears the integer vector before 
   * populating it.
   */
  void getStringVectorFromCell(
    const char* column_name, ///< the column name
    std::vector<std::string>& string_vector, ///<the vector of integers
    char delimiter=',' ///<the delimiter to use
  );

  /**
   * gets an vector of integers from cell where the
   * string in the cell are integers which are separated
   * by a delimiter which is differnt than the column
   * delimiter.  The default delimiter is a comma
   * uses the current_row_ as the row index.
   * clears the integer vector before 
   * populating it.
   */
  void getIntegerVectorFromCell(
    const char* column_name, ///< the column name
    std::vector<int>& int_vector, ///<the vector of integers
    char delimiter=',' ///<the delimiter to use
  );

  /**
   * gets an vector of doubles from cell where the
   * string in the cell are integers which are separated
   * by a delimiter which is differnt than the column
   * delimiter.  The default delimiter is a comma
   * uses the current_row_ as the row index.
   * clears the double vector before 
   * populating it.
   */
  void getDoubleVectorFromCell(
    const char* column_name, ///< the column name
    std::vector<double>& double_vector, ///<the vector of integers
    char delimiter=',' ///<the delimiter to use
  );

  /*Iterator functions.*/
  /**
   * resets the file pointer to the beginning of the file.
   */
  void reset();

  /**
   * parses the next line in the file. 
   */
  void next();

  /**
   * \returns whether there are more rows to 
   * iterate through
   */
  bool hasNext();

  /**
   * converts a datatype to a string
   */
  template<typename TValue>
  static std::string to_string(
    TValue& value ///< the data to convert
  ) {

    std::ostringstream oss;
    oss << std::setprecision(get_int_parameter("precision"));
    oss << value;
    std::string out_string = oss.str();
    return out_string;
  }

};

#endif //DELIMITEDFILEREADER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
