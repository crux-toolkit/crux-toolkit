#include "CnprFragmentMasses.h"

using namespace std;

CnprFragmentMasses::CnprFragmentMasses(){
  channel = -1;
  mz=0;
}

void CnprFragmentMasses::write(FILE* f, int tabs){
  //required
  string el = "fragment_masses";
  if (channel<0) NPRerrMsg(el, "channel");
  if (mz==0) NPRerrMsg(el, "mz");

  NPRprintTabs(f, tabs);
  fprintf(f, "<fragment_masses");
  fprintf(f, " channel=\"%d\"", channel);
  fprintf(f, " mz=\"%.5lf\"", mz);
  fprintf(f, "/>\n");

}