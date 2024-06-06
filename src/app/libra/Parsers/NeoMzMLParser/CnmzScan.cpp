#include "CnmzScan.h"

using namespace std;

void CnmzScan::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "scan";

  NMZprintTabs(f, tabs);
  fprintf(f, "<scan");
  if (!externalSpectrumID.empty()) fprintf(f, " externalSpectrumID=\"%s\"", externalSpectrumID.c_str());
  if (!instrumentConfigurationRef.empty()) fprintf(f, " instrumentConfigurationRef=\"%s\"", instrumentConfigurationRef.c_str());
  if (!sourceFileRef.empty()) fprintf(f, " sourceFileRef=\"%s\"", sourceFileRef.c_str());
  if (!spectrumRef.empty()) fprintf(f, " spectrumRef=\"%s\"", spectrumRef.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referencableParamGroupRef.size(); a++) referencableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);
  for (size_t a = 0; a<scanWindowList.size(); a++) scanWindowList[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</scan>\n");

}
