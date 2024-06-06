#include "CnmzBinaryDataArrayList.h"

using namespace std;

void CnmzBinaryDataArrayList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "binaryDataArrayList";

  NMZprintTabs(f, tabs);
  fprintf(f, "<binaryDataArrayList");
  fprintf(f, " count=\"%d\"", count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<binaryDataArray.size(); a++) binaryDataArray[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</binaryDataArrayList>\n");


}

