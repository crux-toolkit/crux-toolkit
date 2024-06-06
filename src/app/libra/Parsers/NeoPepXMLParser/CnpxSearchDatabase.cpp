#include "CnpxSearchDatabase.h"

using namespace std;

CnpxSearchDatabase::CnpxSearchDatabase(){
  database_name.clear();
  database_release_date.clear();
  database_release_identifier.clear();
  local_path.clear();
  orig_database_url.clear();
  size_in_db_entries=0;
  size_of_residues=0;
  type.clear();
  URL.clear();
}

void CnpxSearchDatabase::write(FILE* f, int tabs){
  string el = "search_database";
  if (local_path.empty()) NPXerrMsg(el, "local_path");
  if (type.empty()) NPXerrMsg(el, "type");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<search_database local_path=\"%s\"", local_path.c_str());
  if (!URL.empty()) fprintf(f, " URL=\"%s\"", URL.c_str());
  if (!database_name.empty()) fprintf(f, " database_name=\"%s\"", database_name.c_str());
  if (!orig_database_url.empty()) fprintf(f, " orig_database_url=\"%s\"", orig_database_url.c_str());
  if (database_release_date.date.day>0) fprintf(f, " database_release_date=\"%s\"", database_release_date.write().c_str());
  if (!database_release_identifier.empty()) fprintf(f, " database_release_identifier=\"%s\"", database_release_identifier.c_str());
  if (size_in_db_entries>0) fprintf(f, " size_in_db_entries=\"%d\"", size_in_db_entries);
  if (size_of_residues>0) fprintf(f, " size_of_residues=\"%d\"", size_of_residues);
  fprintf(f, " type=\"%s\"", type.c_str());
  fprintf(f, "/>\n");

}
