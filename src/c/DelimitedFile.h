/**
 * \file DelimitedFile.h
 * $Revision: 1.00 $ 
 * DATE: Jan 7, 2010
 * AUTHOR: Sean McIlwain
 * \brief Object for reading tab-delimited files.
 ****************************************************************************/
#ifndef DELIMITEDFILE_H
#define DELIMITEDFILE_H

#include <ios>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include "parameter.h"
}


class DelimitedFile {

 protected:
  /** The data vector is indexed by column, then by row
   *  Ex. data_[0][4] means 0th column, 4th row
   *  data_[0] is a string vector of all data in a column.
   */
  std::vector<std::vector<std::string> > data_;
  std::vector<std::string> column_names_;

  /**
   * tokenize a string by delimiter
   */
  void tokenize(
    const std::string& str,
    std::vector<std::string>& tokens,
    char delimiter = '\t'
  );

  /**
   * convert string to data type
   */
  template<typename T>  
  bool from_string(
    T& value,
    const std::string& s
  );

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
   * adds a column to the delimited file
   */
  unsigned int addColumn();

  /**
   * finds the index of a column
   *\returns the column index, -1 if not found.
   */ 
  int findColumn(
    std::string& column_name ///< the column name
  );

 /**
   * finds the index of a column
   *\returns the column index, -1 if not found.
   */ 
  int findColumn(
    const char* column_name ///< the column name
  );

  
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
   *\returns the string value of the cell
   */
  std::string& getString(
    unsigned int col_idx, ///< the column index
    unsigned int row_idx  ///< the row index
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

  /**
   * gets a double type from cell, checks for infinity. 
   */
  double getDouble(
    unsigned int col_idx, ///< the column index 
    unsigned int row_idx ///< the row index
  );

 /**
   * gets an integer type from cell, checks for infinity. 
   */
  int getInteger(
    unsigned int col_idx, ///< the column index 
    unsigned int row_idx ///< the row index
  );

};






#endif //DELIMITEDFILE_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
