/**
 * \file SortColumn.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Given a delimited file and a column-name, sort the 
 * file.
 *****************************************************************************/
#ifndef SORTCOLUMN_H
#define SORTCOLUMN_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>
#include <vector>

class SortColumn: public CruxApplication {

 protected:
  //parameters
  std::string delimited_filename_; ///<delimited filename to sort
  std::string column_name_string_; ///<column name to sort by
  COLTYPE_T column_type_;     ///<column's type (int, real, string).
  bool ascending_;            ///<ascending/descending order?
  char delimiter_;            ///<file's delimiter
  bool header_;               ///<print out the header?
  unsigned int col_sort_idx_; ///<column index to sort by

  //private methods.
  /**
   * sorts the delimited file based upon the column type parameter
   */
  void sortDelimited(
    DelimitedFile* delimited
  );

  /**
   * merges a list of sorted delimited files and prints 
   * out the resulting sorted file.
   */
  void mergeDelimitedFiles(
    std::vector<std::string>& temp_filenames
  );

  /**
   * /returns the result of comparing two file(s) current row and column: 
   * 1 : file1 > file2
   * -1 : file1 < file2
   * 0 : file1 = file2
   */
  int compare(
    DelimitedFileReader* file1,
    DelimitedFileReader* file2
  );


 public:

  /**
   * \returns a blank SortColumn object
   */
  SortColumn();

  /**
   * Destructor
   */
  ~SortColumn();

  /**
   * main method for SortColumn
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for SortColumn
   */
  virtual std::string getName();

  /**
   * \returns the description for SortColumn
   */
  virtual std::string getDescription();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
