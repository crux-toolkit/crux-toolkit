/**
 * \file LineFileReader.h
 * $Revision: 1.00 $ 
 * DATE: Jan 7, 2010
 * AUTHOR: Sean McIlwain
 * \brief Object for reading text files line by line.
 ****************************************************************************/
#ifndef LINEFILEREADER_H
#define LINEFILEREADER_H

#include <limits>
#include <ios>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

//#include "parameter.h"


class LineFileReader {

 protected:
  /** The data vector is indexed by column, then by row
   *  Ex. data_[0][4] means 0th column, 4th row
   *  data_[0] is a string vector of all data in a column.
   */

  std::string next_data_string_; //the next data string.
  std::string current_data_string_; //the current data string.
  unsigned int current_row_; //current row count
  bool has_next_;
  bool has_current_;
  std::fstream* file_ptr_;

  std::string file_name_;

  bool num_rows_valid_;
  unsigned int num_rows_;


 public:
  /**
   * \returns a blank LineFileReader object 
   */
  LineFileReader();
  
  /**
   * \returns a LineFileReader object and loads the tab-delimited
   * data specified by file_name.
   */  
  LineFileReader(
    const char *file_name ///< the path of the file to read 
  );

  /** 
   * \returns a LineFileReader object and loads the tab-delimited
   * data specified by file_name.
   */
  LineFileReader(
    const std::string& file_name ///< the path of the file  to read
  );

  /**
   * Destructor
   */
  virtual ~LineFileReader();

  virtual void loadData(
    const char *file_name ///< the file path
  );

  virtual void loadData(
    const std::string& file_name ///< the file path
  );


  const std::string& current();
  const std::string& next();

  int getCurrentRow();

  bool hasNext();

  void reset();

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


};

#endif //LINEFILEREADER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
