#include "CnpxInputFile.h"

using namespace std;

void CnpxInputFile::write(FILE* f){
  fprintf(f, "<inputfile name=\"%s\"", name.c_str());
  if(directory.size()>0) fprintf(f, " directory=\"%s\"", directory.c_str());
  fprintf(f, "/>\n");
}