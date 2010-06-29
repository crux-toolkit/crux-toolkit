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

// N.B. Compare to the corresponding list in match_collection.cpp.
enum _match_columns {
  SCAN_COL,
  CHARGE_COL,
  SPECTRUM_PRECURSOR_MZ_COL,
  SPECTRUM_NEUTRAL_MASS_COL,
  PEPTIDE_MASS_COL,
  DELTA_CN_COL,
  SP_SCORE_COL,
  SP_RANK_COL,
  XCORR_SCORE_COL,
  XCORR_RANK_COL,
  PVALUE_COL,
  WEIBULL_QVALUE_COL,
  DECOY_XCORR_QVALUE_COL,
  PERCOLATOR_SCORE_COL,
  PERCOLATOR_RANK_COL,
  PERCOLATOR_QVALUE_COL,
  QRANKER_SCORE_COL,
  QRANKER_QVALUE_COL,
  BY_IONS_MATCHED_COL,
  BY_IONS_TOTAL_COL,
  MATCHES_SPECTRUM_COL,
  SEQUENCE_COL,
  CLEAVAGE_TYPE_COL,
  PROTEIN_ID_COL,
  FLANKING_AA_COL,
  UNSHUFFLED_SEQUENCE_COL,
  ETA_COL,
  BETA_COL,
  SHIFT_COL,
  CORR_COL,
  NUMBER_MATCH_COLUMNS
};

typedef enum _match_columns MATCH_COLUMNS_T;

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
     * Destructor
     */
    virtual ~MatchFileReader();

    /**
     * \returns the FLOAT_T value of a cell, checks for infinity
     */
    FLOAT_T getFloat(
      MATCH_COLUMNS_T col_type ///<the column type
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
    std::string& getString(
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


};

#endif //MATCHFILEREADER_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
