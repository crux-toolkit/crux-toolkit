#include "CnpxSearchSummary.h"

using namespace std;

CnpxSearchSummary::CnpxSearchSummary() {
  search_id = 0;
}

CnpxAminoAcidModification* CnpxSearchSummary::addAminoAcidModification(std::string aminoAcid, double massDiff, double mass, std::string variable){
  CnpxAminoAcidModification s;
  s.aminoacid=aminoAcid;
  s.massdiff=massDiff;
  s.mass=mass;
  s.variable=variable;
  aminoacid_modification.push_back(s);
  return &aminoacid_modification.back();
}

CnpxEnzymaticSearchConstraint* CnpxSearchSummary::addEnzymaticSearchConstraint(std::string enzyme, int maxInternalCleavages, int minTermini){
  CnpxEnzymaticSearchConstraint s;
  s.enzyme=enzyme;
  s.max_num_internal_cleavages=maxInternalCleavages;
  s.min_number_termini=minTermini;
  enzymatic_search_constraint.push_back(s);
  return &enzymatic_search_constraint.back();
}

CnpxParameter* CnpxSearchSummary::addParameter(std::string name, std::string value){
  CnpxParameter s;
  s.name=name;
  s.value=value;
  parameter.push_back(s);
  return &parameter.back();
}

CnpxSearchDatabase* CnpxSearchSummary::addSearchDatabase(std::string localPath, std::string type){
  CnpxSearchDatabase s;
  s.local_path = localPath;
  s.type=type;
  search_database.push_back(s);
  return &search_database.back();
}

void CnpxSearchSummary::write(FILE* f, int tabs){
  size_t i;

  string el = "search_summary";
  if (base_name.empty()) NPXerrMsg(el, "base_name");
  if (search_engine.empty()) NPXerrMsg(el, "search_engine");
  if (precursor_mass_type.empty()) NPXerrMsg(el, "precursor_mass_type");
  if (fragment_mass_type.empty()) NPXerrMsg(el, "fragment_mass_type");
  if (search_id == 0) NPXerrMsg(el, "search_id");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<search_summary base_name=\"%s\"", base_name.c_str());
  fprintf(f, " search_engine=\"%s\"", search_engine.c_str());
  if(!search_engine_version.empty()) fprintf(f, " search_engine_version=\"%s\"", search_engine_version.c_str());
  fprintf(f, " precursor_mass_type=\"%s\"", precursor_mass_type.c_str());
  fprintf(f, " fragment_mass_type=\"%s\"", fragment_mass_type.c_str());
  fprintf(f, " search_id=\"%d\"", search_id);
  fprintf(f, ">\n");

  for (i = 0; i<search_database.size(); i++) search_database[i].write(f,t);
  for(i=0;i<enzymatic_search_constraint.size();i++) enzymatic_search_constraint[i].write(f);
  for(i=0;i<aminoacid_modification.size();i++) aminoacid_modification[i].write(f,t);
  for (i = 0; i<terminal_modification.size(); i++) terminal_modification[i].write(f);
  for (i = 0; i < parameter.size(); i++) parameter[i].write(f,t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</search_summary>\n");
}