#include "CnmzFileDescription.h"

using namespace std;

void CnmzFileDescription::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "fileDescription";

  NMZprintTabs(f, tabs);
  fprintf(f, "<fileDescription>\n");

  int t = tabs;
  if (t>-1) t++;
  fileContent.write(f,t,iterative);
  for (size_t a = 0; a<sourceFileList.size(); a++) sourceFileList[a].write(f, t, iterative);
  for (size_t a = 0; a<contact.size(); a++) contact[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</fileDescription>\n");

}

