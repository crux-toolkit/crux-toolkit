/*
Abstract class for a peptide-spectrum match writers
*/

#ifndef PSMWRITER_H
#define PSMWRITER_H

#include "model/MatchCollection.h"
#include "parameter.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

class PSMWriter {

 public:
  // Constructor
  PSMWriter();

  // Destructor
  virtual ~PSMWriter();

  enum MATCH_FILE_TYPE { NONE, PROTEINS, PEPTIDES, PSMS };

  // Pure virtual write, open file, and close file

  virtual void openFile(CruxApplication* application, std::string filename, MATCH_FILE_TYPE type) = 0;

  virtual void write(MatchCollection* collection, std::string database) = 0;

  virtual void closeFile() = 0;
};

#endif
