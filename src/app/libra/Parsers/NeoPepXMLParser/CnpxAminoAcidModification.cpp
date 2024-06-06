#include "CnpxAminoAcidModification.h"

using namespace std;

CnpxAminoAcidModification::CnpxAminoAcidModification(){
  aminoacid.clear();
  binary.clear();
  description.clear();
  massdiff=0;
  mass=0;
  variable.clear();
  peptide_terminus.clear();
  protein_terminus.clear();
  symbol.clear();
}

void CnpxAminoAcidModification::write(FILE* f, int tabs){
  string el = "aminoacid_modification";
  if (aminoacid.empty()) NPXerrMsg(el, "aminoacid");
  if (variable.empty()) NPXerrMsg(el, "variable");
  if (massdiff == 0) NPXerrMsg(el, "massdiff");
  if (mass == 0) NPXerrMsg(el, "mass");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<aminoacid_modification aminoacid=\"%s\"", aminoacid.c_str());
  fprintf(f, " massdiff=\"%.6lf\"", massdiff);
  fprintf(f, " mass=\"%.6lf\"", mass);
  fprintf(f, " variable=\"%s\"", variable.c_str());
  if(!peptide_terminus.empty()) fprintf(f, " peptide_terminus=\"%s\"",peptide_terminus.c_str());
  if (!protein_terminus.empty()) fprintf(f, " protein_terminus=\"%s\"", protein_terminus.c_str());
  if (!symbol.empty()) fprintf(f, " symbol=\"%s\"", symbol.c_str());
  if (!binary.empty()) fprintf(f, " binary=\"%s\"", binary.c_str());
  if (!description.empty()) fprintf(f, " description=\"%s\"", description.c_str());
  fprintf(f, "/>\n");

}
