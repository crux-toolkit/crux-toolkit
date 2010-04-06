/**
 * \file DelimitedFile.h
 * $Revision: 1.00 $ 
 * DATE: Jan 7, 2010
 * AUTHOR: Sean McIlwain
 * \brief Object for reading tab-delimited files.
 * This class generates a table of values. The default delimiter is tab.
 * This class is capable of reading string, integers, and floating point
 * Types from each cell of the table.  This class also provides function
 * for reading a list of integers or string from a cell using a delimiter
 * that is different from the column delimiter (default is comma ',').
 ****************************************************************************/
#ifndef DELIMITEDFILE_H
#define DELIMITEDFILE_H

#include <limits>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "parameter.h"


class DelimitedFile {

 protected:
  /** The data vector is indexed by column, then by row
   *  Ex. data_[0][4] means 0th column, 4th row
   *  data_[0] is a string vector of all data in a column.
   */
  std::vector<std::vector<std::string> > data_;
  std::vector<std::string> column_names_;
  unsigned int current_row_; //used for iterating through the table.

  template <typename T>
  void reorderRows(
    std::multimap<T, unsigned int>& sort_indices, 
    BOOLEAN_T ascending);  

 public:
  /**
   * \returns a blank DelimitedFile object 
   */
  DelimitedFile();
  
  /**
   * \returns a DelimitedFile object and loads the tab-delimited
   * data specified by file_name.
   */  
  DelimitedFile(
    const char *file_name, ///< the path of the file to read 
    bool hasHeader=true ///< indicate whether header exists
  );

  /** 
   * \returns a DelimitedFile object and loads the tab-delimited
   * data specified by file_name.
   */
  DelimitedFile(
    const std::string& file_name, ///< the path of the file  to read
    bool hasHeader=true ///< indicates whether header exists
  );

  /**
   * clears the table
   */
  void clear();

  /**
   * Destructor
   */
  virtual ~DelimitedFile();
  
  /**
   *\returns the number of rows, assuming a square matrix
   */
  unsigned int numRows();

  /**
   *\returns the number of columns
   */
  unsigned int numCols();

  /**
   *\returns the number of rows for a column
   */
  unsigned int numRows(
    unsigned int col_idx ///<the column index
  );

  /**
   * clears the current data and column names,
   * parses the header if it exists,
   * reads the file one line at a time while
   * populating the data matrix with the 
   * strings separated by tabs.
   */
  void loadData(
    const char *file_name, ///< the file path
    bool hasHeader=true ///< header indicator
  );

  /**
   * loads a tab delimited file
   */
  void loadData(
    const std::string& file_name, ///< the file path
    bool hasHeader=true ///< header indicator
  );

  /**
   * saves a tab delimited file
   */ 
  void saveData(
    const char* file_name ///< the file path
  );

  /**
   * saves a tab delimited file
   */ 
  void saveData(
    const std::string& file_name ///< the file path
  );

  /**
   * adds a column to the delimited file
   *\returns the column index.
   */
  unsigned int addColumn(
    std::string& column_name ///< the column name
  );

  /**
   * adds a column to the delimited file
   *\returns the column index.
   */
  unsigned int addColumn(
    const char* column_name ///< the column name
  );

  /**
   * adds a vector of columns to the delimited file
   */
  void addColumns(
    std::vector<std::string>& column_names);

  /**
   * adds a column to the delimited file
   */
  unsigned int addColumn();

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
    * adds a row to the delimited file
    *\returns the new row index
    */
  unsigned int addRow();
 
  /**
   *\returns the string vector corresponding to the column
   */
  std::vector<std::string>& getColumn(
    std::string column ///< the column name
  );

  /**
   *\returns the string vector corresponding to the column
   */
  std::vector<std::string>& getColumn(
    unsigned int col_idx ///< the column index
  );

  /**
   *\returns the name of the column
   */
  std::string& getColumnName(
    unsigned int col_idx ///< the column index
  );
  
  /**
   *\returns the column_names
   */
  std::vector<std::string>& getColumnNames();

  /**
   *\returns the string value of the cell
   */
  std::string& getString(
    unsigned int col_idx, ///< the column index
    unsigned int row_idx  ///< the row index
  );

  /** 
   * gets a string value of the cell.
   */
  std::string& getString(
    const char* column_name, ///<the column name
    unsigned int row_idx ///< the row index
  );

  /**
   * gets a string value of the cell
   * uses the current_row_ as the row index
   */
  std::string& getString(
    const char* column_name ///<the column name
  );


  /**
   * sets the string value of the cell
   */
  void setString(
    unsigned int col_idx, ///< the column index
    unsigned int row_idx, ///< the row index
    std::string& value ///< the new value
  );

  /**
   * sets the string value of the cell
   */
  void setString(
    unsigned int col_idx, ///< the column index
    unsigned int row_idx, ///< the row index
    char* value ///< the new value
  );

  /**
   *\returns the data of the cell in the specified type
   */
  template<typename TValue>
  TValue getValue(
    unsigned int col_idx, ///< the column index 
    unsigned int row_idx  ///< the row index
  );

  /**
   * sets the cell value with the specified type
   */
  template<typename TValue>
  void setValue(
    unsigned int col_idx, ///< the column index 
    unsigned int row_idx, ///< the row index
    TValue value ///< the new value
  ) {
      int precision = get_int_parameter("precision");
      std::ostringstream ss;
      ss << std::setprecision(precision) << value;
      std::string svalue = ss.str();
      setString(col_idx, row_idx, svalue);
  }  

  template<typename TValue>
  void setValue(
    const std::string& column_name, ///< the column index
    unsigned int row_idx, ///< the row index
    TValue value ///< the new value
  ) {

    int col_idx = findColumn(column_name);
    if (col_idx == -1) {
      carp(CARP_FATAL,"Column not found %s",column_name.c_str());
    }
    setValue(col_idx, row_idx, value);
  }

  /**
   * gets a double type from cell, checks for infinity. 
   */
  FLOAT_T getFloat(
    unsigned int col_idx, ///< the column index
    unsigned int row_idx ///< the row index
  );

  /** 
   * gets a double type from cell, checks for infinity.
   */
  FLOAT_T getFloat(
    const char* column_name, ///<the column name
    unsigned int row_idx ///< the row index
  );

  /**
   * gets a double value from cell, checks for infinity
   * uses the current_row_ as the row index
   */
  FLOAT_T getFloat(
    const char* column_name ///<the column name
  );

  /**
   * gets a double type from cell, checks for infinity. 
   */
  double getDouble(
    unsigned int col_idx, ///< the column index 
    unsigned int row_idx ///< the row index
  );

  /** 
   * gets a double type from cell, checks for infinity.
   */
  double getDouble(
    const char* column_name, ///<the column name
    unsigned int row_idx ///<the row index
  );

  /**
   * gets a double value from cell, checks for infinity
   * uses the current_row_ as the row index
   */
  double getDouble(
    const char* column_name ///<the column name
  );


 /**
   * gets an integer type from cell, checks for infinity. 
   */
  int getInteger(
    unsigned int col_idx, ///< the column index 
    unsigned int row_idx ///< the row index
  );

  /**
   * get an integer type from cell, checks for infintiy.
   */
  int getInteger(
    const char* column_name, ///< the column name
    unsigned int row_idx ///<the row index
  );

  /**
   * get an integer type from cell, checks for infinity.
   * uses the current_row_ as the row index.
   */
  int getInteger(
    const char* column_name ///< the column name
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
   * sorts the table by a column. Assumes the data type is
   * Float. By default sorts in ascending order.
   */
  void sortByFloatColumn(
    const std::string& column_name, ///< the column name
    BOOLEAN_T ascending = TRUE);
  
  /**
   * sorts the table by a column. Assumes the data type is 
   * integer.
   */
  void sortByIntegerColumn(
    const std::string& column_name, ///< the column name
    BOOLEAN_T ascending = TRUE);

  void sortByIntegerColumn(
    unsigned int col_idx,
    BOOLEAN_T ascending = TRUE);


  
  /**
   * sorts the table by a column. Assumes the data type is 
   * string.
   */
  void sortByStringColumn(
    const std::string& column_name,
    BOOLEAN_T ascending = TRUE);


  void copyToRow(
    DelimitedFile& dest,
    int src_row_idx, 
    int dest_row_idx
  ); 


  /*Iterator functions.*/
  /**
   * resets the current_row_ index to 0.
   */
  void reset();

  /**
   * increments the current_row_, 
   */
  void next();


  /**
   * \returns whether there are more rows to 
   * iterate through
   */
  BOOLEAN_T hasNext();

  /**
   * tokenize a string by delimiter
   */
  static void tokenize(
    const std::string& str,
    std::vector<std::string>& tokens,
    char delimiter = '\t'
  );

  template<typename T>
  static std::string splice(
    const std::vector<T>& elements,
    char delimiter = '\t') {

      if (elements.size() == 0) return "";

      int precision = get_int_parameter("precision");
      std::ostringstream ss;
      ss << std::setprecision(precision);
      
      ss << elements[0];
      for (int idx=1;idx < elements.size();idx++) {
        ss << delimiter << elements[idx];
      }
      std::string out_string = ss.str();
      return out_string;
  }

  /**
   * convert string to data type
   */
  template<typename TValue>  
  static bool from_string(
    TValue& value,
    const std::string& s
    ) {

    std::istringstream iss(s);
    return !(iss >> std::dec >> value).fail();
  }   
  
  /**
   * Allows object to be printed to a stream
   */
  friend std::ostream &operator<< (std::ostream& os, DelimitedFile& delimited_file); 

};






#endif //DELIMITEDFILE_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
