#include "CnpxIntensity.h"

using namespace std;

CnpxIntensity::CnpxIntensity(){
  channel = -1;
  target_mass = 0;
  absolute = 0;
  normalized = 0;
  reject = false;
}

void CnpxIntensity::write(FILE* f, int tabs){
  string el = "intensity";
  if (channel == -1) NPXerrMsg(el, "channel");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<intensity");
  fprintf(f, " channel=\"%d\"", channel);
  fprintf(f, " target_mass=\"%.3lf\"", target_mass);
  fprintf(f, " absolute=\"%.3lf\"", absolute);
  fprintf(f, " normalized=\"%.3lf\"", normalized);
  if (reject) fprintf(f, " reject=\"1\"");
  fprintf(f, "/>\n");

}