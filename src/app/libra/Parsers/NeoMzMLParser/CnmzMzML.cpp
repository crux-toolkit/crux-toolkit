#include "CnmzMzML.h"

using namespace std;

void CnmzMzML::write(FILE* f, int tabs,bool iterative){
  //required
  string el = "mzML";
  if (version.empty()) NMZerrMsg(el, "version");

  NMZprintTabs(f, tabs);
  fprintf(f, "<mzML");
  fprintf(f, " xmlns=\"%s\"", xmlns.c_str());
  fprintf(f, " xmlns:xsi=\"%s\"", xmlns_xsi.c_str());
  fprintf(f, " xsi:schemaLocation=\"%s\"", xsi_schemaLocation.c_str());
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, " version=\"%s\"", version.c_str());
  fprintf(f, ">\n");

  int t=tabs;
  if(t>-1) t++;
  cvList.write(f,t,iterative);
  fileDescription.write(f,t,iterative);
  for (size_t a = 0; a<referencableParamGroupList.size(); a++) referencableParamGroupList[a].write(f, t, iterative);
  if(!sampleList.empty()) sampleList[0].write(f,t,iterative);
  softwareList.write(f,t,iterative);
  instrumentConfigurationList.write(f,t,iterative);
  dataProcessingList.write(f,t,iterative);
  run.write(f,t,iterative);

  if (iterative) return;

  NMZprintTabs(f, tabs);
  fprintf(f, "</mzML>\n");

}
