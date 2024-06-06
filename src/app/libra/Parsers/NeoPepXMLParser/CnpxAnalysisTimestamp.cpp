#include "CnpxAnalysisTimestamp.h"

using namespace std;

CnpxAnalysisTimestamp::CnpxAnalysisTimestamp(){
  analysis.clear();
  id=0;
  time.clear();
}

void CnpxAnalysisTimestamp::write(FILE* f){

  fprintf(f, "<analysis_timestamp analysis=\"%s\"", analysis.c_str());
  fprintf(f, " time=\"%s\"", time.write().c_str());
  fprintf(f, " id=\"%d\">\n", id);

  if (database_refresh_timestamp.present()) database_refresh_timestamp.write(f);


  fprintf(f, "</analysis_timestamp>\n");

}
