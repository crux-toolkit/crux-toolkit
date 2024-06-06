#include "CnpxUIPSM.h"

using namespace std;

void CnpxUIPSM::setPSM(CnpxSpectrumQuery& s){
  size_t i,j;
  CnpxSearchHit* sh = &s.search_result[0].search_hit[0];

  assumed_charge = s.assumed_charge;
  precursorNeutralMass = s.precursor_neutral_mass;
  retentionTime = s.retention_time_sec;
  scanNumber = s.start_scan;
  spectrumID = s.spectrum;

  if(sh->xlink_type.compare("loop")==0) xlType=1;
  else if (sh->xlink_type.compare("xl")==0) xlType=2;
  else xlType=0;

  calcNeutralMass = sh->calc_neutral_pep_mass;
  peptide = sh->peptide;

  npxUIProtein p;
  proteins.clear();
  p.protein = sh->protein;
  p.nextAA = sh->peptide_next_aa[0];
  p.prevAA = sh->peptide_prev_aa[0];
  proteins.push_back(p);
  for(i=0;i<sh->alternative_protein.size();i++){
    p.protein = sh->alternative_protein[i].protein;
    p.nextAA = sh->alternative_protein[i].peptide_next_aa[0];
    p.prevAA = sh->alternative_protein[i].peptide_prev_aa[0];
    proteins.push_back(p);
  }

  modified=false;
  mod.modifiedPeptide.clear();
  mod.mods.clear();
  if (sh->modification_info.size()>0){
    modified=true;
    mod.modifiedPeptide=sh->modification_info[0].modified_peptide;
    npxUIMod m;
    for(i=0;i<sh->modification_info[0].mod_aminoacid_mass.size();i++){
      m.mass = sh->modification_info[0].mod_aminoacid_mass[i].mass;
      m.pos = sh->modification_info[0].mod_aminoacid_mass[i].position;
      m.diffMass=0;
      if (sh->modification_info[0].mod_aminoacid_mass[i].variable>0) m.diffMass = sh->modification_info[0].mod_aminoacid_mass[i].variable;
      if (sh->modification_info[0].mod_aminoacid_mass[i].staticMass>0) m.diffMass = sh->modification_info[0].mod_aminoacid_mass[i].staticMass;
      mod.mods.push_back(m);
    }
  }

  npxUIScore sc;
  scores.clear();
  for(i=0;i<sh->search_score.size();i++){
    sc.name=sh->search_score[i].name;
    sc.value=sh->search_score[i].value;
    scores.push_back(sc);
  }

  hasPeptideProphet=false;
  hasIProphet=false;
  peptideProphet.probability=0;
  peptideProphet.parameters.clear();
  iProphet.probability=0;
  iProphet.parameters.clear();
  for(i=0;i<sh->analysis_result.size();i++){
    if(sh->analysis_result[i].peptide_prophet_result.present()){
      hasPeptideProphet=true;
      peptideProphet.probability=sh->analysis_result[i].peptide_prophet_result.probability;
      for (j = 0; j<sh->analysis_result[i].peptide_prophet_result.search_score_summary.parameter.size();j++){
        sc.name = sh->analysis_result[i].peptide_prophet_result.search_score_summary.parameter[j].name;
        sc.value = sh->analysis_result[i].peptide_prophet_result.search_score_summary.parameter[j].value;
      }
    }
    if (sh->analysis_result[i].interprophet_result.present()){
      hasIProphet = true;
      iProphet.probability = sh->analysis_result[i].interprophet_result.probability;
      for (j = 0; j<sh->analysis_result[i].interprophet_result.search_score_summary.parameter.size(); j++){
        sc.name = sh->analysis_result[i].interprophet_result.search_score_summary.parameter[j].name;
        sc.value = sh->analysis_result[i].interprophet_result.search_score_summary.parameter[j].value;
      }
    }
  }  

  sh = NULL;
}
