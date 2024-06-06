#ifndef _CNPXSEARCHDATABASE_H
#define _CNPXSEARCHDATABASE_H

#include "NeoPepXMLStructs.h"
#include <iostream>
#include <string>

class CnpxSearchDatabase {
public:
  CnpxSearchDatabase();

  void write(FILE* f, int tabs=-1);

  std::string database_name;
  npxDateTime database_release_date;
  std::string database_release_identifier;
  std::string local_path;
  std::string orig_database_url;
  int size_in_db_entries;
  int size_of_residues;
  std::string type;
  std::string URL;

private:

};

#endif