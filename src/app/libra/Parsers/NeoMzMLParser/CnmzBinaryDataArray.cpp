#include "CnmzBinaryDataArray.h"

using namespace std;

void CnmzBinaryDataArray::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "binaryDataArray";

  NMZprintTabs(f, tabs);
  fprintf(f, "<binaryDataArray");
  fprintf(f, " encodedLength=\"%d\"", encodedLength);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  binary.write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</binaryDataArray>\n");


}

