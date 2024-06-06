#include "CnpxAffectedChannel.h"

using namespace std;

CnpxAffectedChannel::CnpxAffectedChannel(){
  channel = -1;
  correction = 0;
}

void CnpxAffectedChannel::write(FILE* f, int tabs){
  string el = "affected_channel";
  if (channel == -1) NPXerrMsg(el, "channel");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<affected_channel");
  fprintf(f, " channel=\"%d\"", channel);
  fprintf(f, " correction=\"%.3lf\"", correction);
  fprintf(f, "/>\n");

}