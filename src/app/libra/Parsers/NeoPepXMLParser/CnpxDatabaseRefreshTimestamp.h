#ifndef _CNPXDATABASEREFRESHTIMESTAMP_H
#define _CNPXDATABASEREFRESHTIMESTAMP_H

#include <string>

class CnpxDatabaseRefreshTimestamp {
public:
  CnpxDatabaseRefreshTimestamp();
  CnpxDatabaseRefreshTimestamp(bool b);

  bool present();
  void write(FILE* f);

  std::string database;
  int min_num_enz_term;

private:
  bool active;

};

#endif