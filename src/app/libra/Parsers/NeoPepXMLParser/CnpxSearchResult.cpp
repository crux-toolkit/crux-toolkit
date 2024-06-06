#include "CnpxSearchResult.h"

using namespace std;

CnpxSearchResult::CnpxSearchResult(){
  search_id=0;
}

CnpxSearchHit* CnpxSearchResult::addSearchHit(int hitRank, std::string peptide, std::string protein, int numTotProteins, double calcNeutPepMass, double massdiff){
  CnpxSearchHit  sh;
  sh.hit_rank=hitRank;
  sh.peptide=peptide;
  sh.protein=protein;
  sh.num_tot_proteins=numTotProteins;
  sh.calc_neutral_pep_mass=calcNeutPepMass;
  sh.massdiff=massdiff;
  search_hit.push_back(sh);
  return &search_hit.back();
}

void CnpxSearchResult::write(FILE* f, int tabs){
  size_t i;

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<search_result");
  if (search_id>0) fprintf(f, " search_id=\"%d\"", search_id);
  fprintf(f, ">\n");

  for (i = 0; i<search_hit.size(); i++) search_hit[i].write(f,t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</search_result>\n");

}
