/**
 * \file MatchFileWriter.h
 * DATE: October 26, 2010
 * AUTHOR: Barbara Frewen
 * \brief Object for writing tab-deliminted text files of PSMs (matches).
 * This class is the extension of DelimitedFileWriter where the
 * columns are those in MatchColumns.  Which columns are written to
 * file depend on the COMMAND_TYPE_T.  For some commands, they also
 * depend on the columns found in the input file as given by a
 * MatchFileReader. 
 */

#ifndef MATCH_FILE_WRITER_H
#define MATCH_FILE_WRITER_H

#include "DelimitedFileWriter.h"
#include "MatchColumns.h"
#include "objects.h"
#include <iostream>

class MatchFileWriter : public DelimitedFileWriter {
 protected:
  bool match_to_print_[NUMBER_MATCH_COLUMNS];///< set before writing header
  int match_indices_[NUMBER_MATCH_COLUMNS];///<idx of each MATCH_COLUMN in file
  int match_precision_[NUMBER_MATCH_COLUMNS];///< precision for each column
  bool match_fixed_float_[NUMBER_MATCH_COLUMNS]; ///< do we use fixed format for the float field?
  unsigned int num_columns_; ///< number of columns being printed

  void setPrecision();

 public:

  /**
   * \returns A blank MatchFileWriter object.
   */
  MatchFileWriter();

  /**
   * \returns A MatchFileWriter object and opens a file for writing.
   */
  MatchFileWriter(const char* filename);

  /**
   * Destructor
   */
  ~MatchFileWriter();

  // private?
  /**
   * Defines the columns to print based on the vector of flags
   * indiciating if the MATCH_COLUMN_T should be printed.
   */
  void addColumnNames(const std::vector<bool>& col_is_printed);

  // private?
  /**
   * Adds one column to print.
   */
  void addColumnName(MATCH_COLUMNS_T col_type);

  /**
   * Adds which columns to print based on the COMMAND_TYPE_T. Only for
   * search-for-matches and sequest-search.
   */
  void addColumnNames(CruxApplication* application, bool has_decoys);

  /**
   * Adds which columns to print based on the COMMAND_TYPE_T and a list
   * of columns to print. For all post-search commands.
   */
  void addColumnNames
    (CruxApplication* application, 
     bool has_decoys,
     const std::vector<bool>& cols_to_print);

  /**
   * Write header to file based on 
   */
  virtual void writeHeader();

  /**
   * Set the value in the current row for the given MATCH_COLUMN_T.
   */
  template<typename ValueType>
    void setColumnCurrentRow
    (MATCH_COLUMNS_T col_type,
     const ValueType& value){

    int file_column = match_indices_[col_type];
    // ignore if this column isn't being printed
    if( file_column == -1 ){
      return;
    }
    current_row_.at(file_column) = 
      DelimitedFileWriter::to_string(value, match_precision_[col_type], match_fixed_float_[col_type]);
  }

};

#endif // MATCH_FILE_WRITER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

