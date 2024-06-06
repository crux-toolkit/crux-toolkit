#include "CnmzSpectrum.h"

using namespace std;

void CnmzSpectrum::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "spectrum";
  if (id.empty()) NMZerrMsg(el, "id");

  NMZprintTabs(f, tabs);
  fprintf(f, "<spectrum");
  fprintf(f, " index=\"%d\"", index);
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, " defaultArrayLength=\"%d\"", defaultArrayLength);
  if (!dataProcessingRef.empty()) fprintf(f, " dataProcessingRef=\"%s\"", dataProcessingRef.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referencableParamGroupRef.size(); a++) referencableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);
  for (size_t a = 0; a<scanList.size(); a++) scanList[a].write(f, t, iterative);
  for (size_t a = 0; a<precursorList.size(); a++) precursorList[a].write(f, t, iterative);
  for (size_t a = 0; a<binaryDataArrayList.size(); a++) binaryDataArrayList[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</spectrum>\n");


}

