/**
 * \file PMCDelimitedFileWriter.h
 * $Revision$
 * \brief Writes out delimited files from ProteinMatchCollection objects.
 */

#ifndef PROTEINMATCHCOLLECTION_DELIMITED_FILE_WRITER_H_
#define PROTEINMATCHCOLLECTION_DELIMITED_FILE_WRITER_H_

#include <vector>

#include "model/match_objects.h"
#include "DelimitedFile.h"
#include "MatchFileWriter.h"
#include "model/Peptide.h"
#include "model/PeptideMatch.h"
#include "model/PostProcessProtein.h"
#include "model/ProteinMatch.h"
#include "model/ProteinMatchCollection.h"
#include "PSMWriter.h"
#include "model/Spectrum.h"
#include "model/SpectrumMatch.h"

using namespace std;

class PMCDelimitedFileWriter : public MatchFileWriter, public PSMWriter {

 public:
  /**
   * Defines the type of match file being written
   */

  /**
   * Returns an empty PMCDelimitedFileWriter object
   */
  PMCDelimitedFileWriter();

  /**
   * Destructor
   */
  ~PMCDelimitedFileWriter();

  /**
   * Opens a file, writes it, and closes the file
   */
  void writeFile(
    CruxApplication* application, ///< application writing the file
    string filename, ///< name of the file to open
    MATCH_FILE_TYPE type, ///< type of file to be written
    ProteinMatchCollection* collection ///< collection to be written
  );

  /**
   * Closes any open file, then opens a file for the specified type of writing
   */
  void openFile(
    CruxApplication* application, ///< application writing the file
    string filename, ///< name of the file to open
    MATCH_FILE_TYPE type ///< type of file to be written
  );

  /**
   * Closes any open file, if any
   */
  void closeFile();

  /**
   * Writes all three file types using the file stem as a base filename
   */
  void writeAll(
    CruxApplication* application, ///< application writing the files
    ProteinMatchCollection* collection, ///< collection to be written
    string stem ///< filestem to be prepended to the filenames
  );

  void write(
    MatchCollection* collection,
    string database
  );

  /**
   * Writes the data in a ProteinMatchCollection to the currently open file
   */
  void write(
    ProteinMatchCollection* collection ///< collection to be written
  );

  /**
    * Sets whether we should write our output in tab delimited or HTML format. Default = false (tab delimited)
    */
  void setWriteHTML(bool write_html);

 private:

  bool write_html_; // this determines whether we want to write in html format. Default value = false;
  CruxApplication* application_; // pointer to the application using this program

  // function pointer to the appropriate writing function for the current file type
  void (PMCDelimitedFileWriter::*write_function_)(ProteinMatchCollection*);

  /**
   * Set up the columns for a protein match file
   */
  void setUpProteinsColumns(
    CruxApplication* application ///< application writing the file
  );

  /**
   * Write protein match file
   */
  void writeProteins(
    ProteinMatchCollection* collection ///< collection to be written
  );

  /**
   * Set up the columns for a peptide match file
   */
  void setUpPeptidesColumns(
    CruxApplication* application ///< application writing the file
  );

  /**
   * Write peptide match file
   */
  void writePeptides(
    ProteinMatchCollection* collection ///< collection to be written
  );

  /**
   * Set up the columns for a spectrum match file
   */
  void setUpPSMsColumns(
    CruxApplication* application ///< application writing the file
  );

  /**
   * Write spectrum match file
   */
  void writePSMs(
    ProteinMatchCollection* collection ///< collection to be written
  );

  /**
   * Gets the cleavage type as a string
   */
  string getCleavageType();

  /**
   * Adds a column to the header if match has score
   */
  void addScoreColumnIfExists(
    AbstractMatch* match, ///< match to get score from
    SCORER_TYPE_T scoreType, ///< score type to get
    MATCH_COLUMNS_T column ///< column to add the score to
  );

  /**
   * Adds a column to the header if match has rank
   */
  void addRankColumnIfExists(
    AbstractMatch* match, ///< match to get rank from
    SCORER_TYPE_T scoreType, ///< rank type to get
    MATCH_COLUMNS_T column ///< column to add the rank to
  );

  /**
   * Adds the value of a score to the specified column, if it exists in the match
   */
  void addScoreIfExists(
    AbstractMatch* match, ///< match to get score from
    SCORER_TYPE_T scoreType, ///< score type to get
    MATCH_COLUMNS_T column ///< column to add the score to
  );

  /**
   * Adds the value of a rank to the specified column, if it exists in the match
   */
  void addRankIfExists(
    AbstractMatch* match, ///< match to get rank from
    SCORER_TYPE_T scoreType, ///< rank type to get
    MATCH_COLUMNS_T column ///< column to add the rank to
  );

  /**
   * Adds the value of a char* to the specified column, then frees the pointer
   */
  void setAndFree(
    MATCH_COLUMNS_T column, ///< column to add the char* to
    char* value ///< char* to set the column value to
  );

  /**
   * Writes the header in HTML format, this is used for HTMLWriter instead of
   * DelimitedFileWriter's writeHeader()
   */
  void writeHTMLHeader();

  /**
   * Writes a row in HTML format, this is used for HTMLWriter instead of
   * DelimitedFileWriter's writeRow() 
   */
  void writeHTMLRow();

};

#endif // PROTEINMATCHCOLLECTION_DELIMITED_FILE_WRITER_H_

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
