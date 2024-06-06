#include "CnmzRun.h"

using namespace std;

void CnmzRun::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "run";
  if (defaultInstrumentConfigurationRef.empty()) NMZerrMsg(el, "defaultInstrumentConfigurationRef");
  if (id.empty()) NMZerrMsg(el, "id");

  NMZprintTabs(f, tabs);
  fprintf(f, "<run");
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, " defaultInstrumentConfigurationRef=\"%s\"", defaultInstrumentConfigurationRef.c_str());
  if (!defaultSourceFileRef.empty()) fprintf(f, " defaultSourceFileRef=\"%s\"", defaultSourceFileRef.c_str());
  if (!sampleRef.empty()) fprintf(f, " sampleRef=\"%s\"", sampleRef.c_str());
  if (startTimeStamp.date.day>0) fprintf(f, " startTimeStamp=\"%s\"", startTimeStamp.write().c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  spectrumList.write(f, t, iterative);

  if (iterative) return;

  NMZprintTabs(f, tabs);
  fprintf(f, "</run>\n");

}
