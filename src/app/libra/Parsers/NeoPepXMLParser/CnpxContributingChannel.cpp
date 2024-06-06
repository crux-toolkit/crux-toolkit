#include "CnpxContributingChannel.h"

using namespace std;

CnpxContributingChannel::CnpxContributingChannel(){
  channel = -1;
}

void CnpxContributingChannel::write(FILE* f, int tabs){
  string el = "contributing_channel";
  if (channel == -1) NPXerrMsg(el, "contributing_channel");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<contributing_channel");
  fprintf(f, " channel=\"%d\"", channel);
  fprintf(f, ">\n");

  for (size_t i = 0; i<affected_channel.size(); i++) affected_channel[i].write(f, t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</contributing_channel>\n");

}