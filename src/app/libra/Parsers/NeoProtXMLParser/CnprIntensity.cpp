#include "CnprIntensity.h"

using namespace std;

CnprIntensity::CnprIntensity(){
  channel = -1;
  mz=0;
  ratio=0;
  error=0;
}

void CnprIntensity::write(FILE* f, int tabs){
  //required
  string el = "intensity";
  if (channel<0) NPRerrMsg(el, "channel");

  NPRprintTabs(f, tabs);
  fprintf(f, "<intensity");
  fprintf(f, " channel=\"%d\"", channel);
  fprintf(f, " mz=\"%.5lf\"", mz);
  fprintf(f, " mz=\"%.2lf\"", ratio);
  fprintf(f, " mz=\"%.2lf\"", error);
  fprintf(f, "/>\n");

}
