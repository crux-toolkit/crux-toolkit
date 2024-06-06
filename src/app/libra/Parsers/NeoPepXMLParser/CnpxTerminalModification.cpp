#include "CnpxTerminalModification.h"

using namespace std;

CnpxTerminalModification::CnpxTerminalModification(){
  description.clear();
  massdiff = 0;
  mass = 0;
  variable.clear();
  terminus.clear();
  protein_terminus.clear();
  symbol.clear();
}

void CnpxTerminalModification::write(FILE* f){

  fprintf(f, "<terminal_modification terminus=\"%s\"", terminus.c_str());
  fprintf(f, " massdiff=\"%.6lf\"", massdiff);
  fprintf(f, " mass=\"%.6lf\"", mass);
  fprintf(f, " variable=\"%s\"", variable.c_str());
  if (symbol.size()>0) fprintf(f, " symbol=\"%s\"", symbol.c_str());
  if (protein_terminus.size()>0) fprintf(f, " protein_terminus=\"%s\"", protein_terminus.c_str());
  if (description.size()>0) fprintf(f, " description=\"%s\"", description.c_str());
  fprintf(f, "/>\n");

}
