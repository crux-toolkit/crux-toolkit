/*
Abstract class for a peptide-spectrum match readers
*/

#ifndef PSMREADER_H
#define PSMREADER_H

#include "MatchCollection.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <gflags/gflags.h>
#include <string>

class PSMReader {

public:

  // Constructors
  PSMReader();
  PSMReader(const std::string& file_path);
  PSMReader(const std::string& file_path, Database* database,
    Database* decoy_database=NULL);

  // Destructor
  virtual ~PSMReader();

  // Pure Virtual Parse Functions

  virtual MatchCollection* parse() = 0;

/*  virtual MatchCollection* parse(const char* file_path,
    Database* database,
    Database* decoy_database) = 0;*/

  // Methods usable by all Readers

  void setDatabase(Database* database);

  void setDecoyDatabase(Database* decoy_database);

protected:
  Database* database_;
  Database* decoy_database_;
  std::string file_path_;
};

#endif
