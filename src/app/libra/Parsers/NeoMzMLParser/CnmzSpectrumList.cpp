#include "CnmzSpectrumList.h"

using namespace std;

void CnmzSpectrumList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "spectrumList";
  if (defaultDataProcessingRef.empty()) NMZerrMsg(el, "defaultDataProcessingRef");

  NMZprintTabs(f, tabs);
  fprintf(f, "<spectrumList");
  if(iterative) {
    fprintf(f, " count=\"");
    fptr=nmzftell(f);
    fprintf(f,"%.10d\"",0);
  } else fprintf(f, " count=\"%d\"",count);
  fprintf(f, " defaultDataProcessingRef=\"%s\"", defaultDataProcessingRef.c_str());
  fprintf(f, ">\n");

  if(iterative) return;

  int t = tabs;
  if (t>-1) t++;
  //for(size_t a=0;a<spectrum.size();a++) spectrum[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</spectrumList>\n");

}

void CnmzSpectrumList::writeUpdate(FILE* f){
  nmzfseek(f, fptr, 0);
  //cout << "Update: " << count << endl;
  fprintf(f, "%.10d", count);
}
