#ifndef HTMLWRITER_H
#define HTMLWRITER_H

#include <vector>

#include "PMCDelimitedFileWriter.h"

class HTMLWriter : public PMCDelimitedFileWriter {

public:
  /*
   * Re-implemenation of writeRow from DelimitedFileWriter in HTML format
   */
  void writeRow();

  /*
   * Re-implementation of writeHeader from MatchFileWriter/DelimitedFileWriter in HTML fromat
   */
  void writeHeader();

  /**
   * Writes the data in a ProteinMatchCollection to the currently open file
   */
  void write(
    ProteinMatchCollection* collection ///< collection to be written
  );

  void write(
    MatchCollection* collection,
    std::string database
  );

  /**
   * Write spectrum match file
   */
  void writePSMs(
    ProteinMatchCollection* collection ///< collection to be written
  );

};

#endif
