/**
 * \file MatchFileReader.h
 * $Revision: 1.00 $ 
 * DATE: June 10, 2010
 * AUTHOR: Sean McIlwain
 * \brief Object for reading tab-delimited files.
 * This class generates a table of values. The default delimiter is tab.
 * This class is capable of reading string, integers, and floating point
 * Types from each cell of the table.  This class also provides function
 * for reading a list of integers or string from a cell using a delimiter
 * that is different from the column delimiter (default is comma ',').
 * MatchFileReader hashes the column_id for quick access.
 **************************************************************************/

#ifndef MATCHFILEREADER_H
#define MATCHFILEREADER_H

#include "DelimitedFileReader.h"
#include "MatchColumns.h"


class MatchFileReader: public DelimitedFileReader {
  protected:
    int match_indices_[NUMBER_MATCH_COLUMNS];

    void parseHeader();

  public:

   /**
    * \returns a blank MatchFileReader object 
    */
    MatchFileReader();

   /**
    * \returns a MatchFileReader object and loads the tab-delimited
    * data specified by file_name.
    */  
    MatchFileReader(
      const char* file_name
    );

    /** 
     * \returns a MatchFileReader object and loads the tab-delimited
     * data specified by file_name.
     */
    MatchFileReader(
      const std::string& file_name
    );

    /**
     * \returns a MatchFileReader object and load the tab-delimited
     * data specified by an input stream
     */
    MatchFileReader(
      std::istream* iptr
    );

    /**
     * Destructor
     */
    virtual ~MatchFileReader();

    /**
     * Open a new file from an existing MatchFileReader.
     */
    virtual void loadData(
      const char* file_name, ///< new file to open
      bool hasHeader = true  ///< parse first line as header
    );

    /**
     * Open a new file from an existing MatchFileReader.
     */
    virtual void loadData(
      const std::string& file_name, ///< new file to open
      bool hasHeader = true         ///< parse first line has header
    );

    /**
     * \returns the FLOAT_T value of a cell, checks for infinity
     */
    FLOAT_T getFloat(
      MATCH_COLUMNS_T col_type ///<the column type
    );
   
    /**
     * \returns the double value of a cell, checks for infinity
     */
    double getDouble(
      MATCH_COLUMNS_T col_type ///< the column type
    );

    /**
     * \returns the integer value of a cell, checks for infinity.
     */
    int getInteger(
      MATCH_COLUMNS_T col_type ///< the column type
    );

    /**
     * gets a string value of the cell
     */
    const std::string& getString(
      MATCH_COLUMNS_T col_type ///<the column type
    );

    /**
     * returns whether the column is empty or not.
     */
    bool empty(
      MATCH_COLUMNS_T col_type ///<the column type
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
      MATCH_COLUMNS_T col_type, ///< the column name
      std::vector<std::string>& string_vector, ///<the vector of integers
      char delimiter=',' ///<the delimiter to use
    );

    /**
     * Fills in the given vector with a bool value indicating if each
     * MATCH_COLUMN_T type is present in the file being read.
     * \returns Argument vector has NUM_MATCH_COLUMN_T values if a
     * valid file is open and header has been parsed, else vector is empty.
     */
    void getMatchColumnsPresent (std::vector<bool>& col_is_present);


    static MatchCollection* parse(
      const char* file_path,
      Database* database,
      Database* decoy_database
    );

};

#endif //MATCHFILEREADER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
