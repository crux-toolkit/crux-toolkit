/**
 * \file PMCDelimitedFileWriter.h
 * $Revision$
 * \brief Writes out delimited files from ProteinMatchCollection objects.
 */

#ifndef PROTEINMATCHCOLLECTION_DELIMITED_FILE_WRITER_H_
#define PROTEINMATCHCOLLECTION_DELIMITED_FILE_WRITER_H_

#include <vector>

#include "crux-file-utils.h"
#include "match_objects.h"
#include "DelimitedFile.h"
#include "MatchFileWriter.h"
#include "Peptide.h"
#include "PeptideMatch.h"
#include "PostProcessProtein.h"
#include "ProteinMatch.h"
#include "ProteinMatchCollection.h"
#include "Spectrum.h"
#include "SpectrumMatch.h"

using namespace std;

class PMCDelimitedFileWriter : public MatchFileWriter {

public:
  /**
   * Defines the type of match file being written
   */
  enum MATCH_FILE_TYPE { NONE, PROTEINS, PEPTIDES, PSMS };

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

  /**
   * Writes the data in a ProteinMatchCollection to the currently open file
   */
  void write(
    ProteinMatchCollection* collection ///< collection to be written
  );

private:
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

};

#endif // PROTEINMATCHCOLLECTION_DELIMITED_FILE_WRITER_H_

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
