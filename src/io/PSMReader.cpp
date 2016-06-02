/*
Abstract class for a peptide-spectrum match readers
*/

#include <cstdio>
#include "PSMReader.h"
#include "parameter.h"

using namespace std;

PSMReader::PSMReader()
  : database_(NULL), decoy_database_(NULL) {
}

PSMReader::PSMReader(const string& file_path, Database* database, Database* decoy_database)
  : file_path_(file_path), database_(database), decoy_database_(decoy_database) {
}

PSMReader::~PSMReader() {
}

void PSMReader::setDatabase(Database* database) {
  database_ = database;
}

void PSMReader::setDecoyDatabase(Database* decoy_database) {
  decoy_database_ = decoy_database;
}

