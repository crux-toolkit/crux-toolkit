#include "CnprAffectedChannel.h"

using namespace std;

CnprAffectedChannel::CnprAffectedChannel(){
  channel = -1;
  correction = 0;
}

void CnprAffectedChannel::write(FILE* f, int tabs){
  //required
  string el = "affected_channel";
  if (channel<0) NPRerrMsg(el, "channel");

  NPRprintTabs(f, tabs);
  fprintf(f, "<affected_channel");
  fprintf(f, " channel=\"%d\"", channel);
  fprintf(f, " correction=\"%.3lf\"", correction);
  fprintf(f, "/>\n");

}