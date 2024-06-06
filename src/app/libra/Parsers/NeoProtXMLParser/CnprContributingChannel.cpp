#include "CnprContributingChannel.h"

using namespace std;

CnprContributingChannel::CnprContributingChannel(){
  channel=-1;
}

void CnprContributingChannel::write(FILE* f, int tabs){
  //required
  string el = "contributing_channel";
  if (channel<0) NPRerrMsg(el, "channel");

  NPRprintTabs(f, tabs);
  fprintf(f, "<contributing_channel");
  fprintf(f, " channel=\"%d\"", channel);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<affected_channel.size(); j++) affected_channel[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</contributing_channel>\n");
}