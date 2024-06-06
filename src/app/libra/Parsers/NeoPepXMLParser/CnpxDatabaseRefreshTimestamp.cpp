#include "CnpxDatabaseRefreshTimestamp.h"

using namespace std;

CnpxDatabaseRefreshTimestamp::CnpxDatabaseRefreshTimestamp() {
  min_num_enz_term = 0;
  database.clear();
  active = false;
}

CnpxDatabaseRefreshTimestamp::CnpxDatabaseRefreshTimestamp(bool b) {
  min_num_enz_term = 0;
  database.clear();
  active = b;
}

bool CnpxDatabaseRefreshTimestamp::present() {
  return active;
}

void CnpxDatabaseRefreshTimestamp::write(FILE* f) {
  fprintf(f, "<database_refresh_timestamp database=\"%s\"", database.c_str());
  if(min_num_enz_term>0) fprintf(f, " min_num_enz_term=\"%d\"", min_num_enz_term);
  fprintf(f, "/>\n");
}