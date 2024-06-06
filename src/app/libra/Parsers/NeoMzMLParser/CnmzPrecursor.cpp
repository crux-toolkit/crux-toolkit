#include "CnmzPrecursor.h"

using namespace std;

void CnmzPrecursor::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "precursor";

  NMZprintTabs(f, tabs);
  fprintf(f, "<precursor");
  if (!externalSpectrumID.empty()) fprintf(f, " externalSpectrumID=\"%s\"", externalSpectrumID.c_str());
  if (!sourceFileRef.empty()) fprintf(f, " sourceFileRef=\"%s\"", sourceFileRef.c_str());
  if (!spectrumRef.empty()) fprintf(f, " spectrumRef=\"%s\"", spectrumRef.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<isolationWindow.size(); a++) isolationWindow[a].write(f, t, iterative);
  for (size_t a = 0; a<selectedIonList.size(); a++) selectedIonList[a].write(f, t, iterative);
  activation.write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</precursor>\n");

}
