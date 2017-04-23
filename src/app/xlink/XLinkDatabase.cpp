#include "XLinkDatabase.h"
#include "util/modifications.h"
#include "model/ModifiedPeptidesIterator.h"
#include "util/GlobalParams.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

using namespace std;

Database* XLinkDatabase::protein_database_;
XLinkBondMap XLinkDatabase::bondmap_;

vector<vector<Crux::Peptide*> > XLinkDatabase::target_peptides_; ///< all peptides generated with no additional missed cleavage;

vector<vector<Crux::Peptide*> > XLinkDatabase::decoy_peptides_;

std::vector<LinearPeptide> XLinkDatabase::target_linear_peptides_;
std::vector<LinearPeptide> XLinkDatabase::decoy_linear_peptides_;

std::vector<MonoLinkPeptide> XLinkDatabase::target_monolink_peptides_;
std::vector<SelfLoopPeptide> XLinkDatabase::target_selfloop_peptides_;
std::vector<SelfLoopPeptide> XLinkDatabase::decoy_selfloop_peptides_;

std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides_;
std::vector<XLinkablePeptide> XLinkDatabase::decoy_xlinkable_peptides_;

std::vector<XLinkablePeptide> XLinkDatabase::decoy_xlinkable_peptides2_;
std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides2_;

std::vector<XLinkablePeptide> XLinkDatabase::target_xlinkable_peptides_flatten_;
std::vector<XLinkablePeptide> XLinkDatabase::decoy_xlinkable_peptides_flatten_;

bool XLinkDatabase::addPeptideToDatabase(Crux::Peptide* peptide) {
  
  bool added = false;
  vector<int> link_sites;
  set<int> skip;
  
  bool dead = GlobalParams::getXLinkIncludeDeadends();
  bool linear = GlobalParams::getXLinkIncludeLinears();
    
  if (dead || linear) {
    set<int> skip;
    if (peptide->getMissedCleavageSites(skip) <= GlobalParams::getMissedCleavages()) {
      if (peptide->hasMonoLink()) {  //this is implicit
        MonoLinkPeptide mpeptide(peptide);
        mpeptide.getMass(GlobalParams::getIsotopicMass());
        target_monolink_peptides_.push_back(mpeptide);
        added = true;
      } else if (linear) {
        LinearPeptide lpeptide(peptide);
        lpeptide.getMass(GlobalParams::getIsotopicMass());
        target_linear_peptides_.push_back(lpeptide);
        added = true;
      }
    }
  }
  
  if (GlobalParams::getXLinkIncludeInterIntra() ||
      GlobalParams::getXLinkIncludeInter() ||
      GlobalParams::getXLinkIncludeIntra() ||
      GlobalParams::getXLinkIncludeDeadends()) {
    if (peptide->countModifiedAAs() <= GlobalParams::getMaxXLinkMods()) {
      XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites); 
      if (!link_sites.empty()) {
	      XLinkablePeptide xlp(peptide, link_sites);
	      xlp.getMass(GlobalParams::getIsotopicMass());
        target_xlinkable_peptides_.push_back(xlp);
        added = true;
      }
    }
  }
  if (GlobalParams::getXLinkIncludeSelfloops()) {
    XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites, 1);
    if (link_sites.size() > 1) {
      XLinkablePeptide xlp = XLinkablePeptide(peptide, link_sites);
      for (size_t link1_idx = 0; link1_idx<xlp.numLinkSites()-1;link1_idx++) {
	      for (size_t link2_idx = link1_idx+1; link2_idx < xlp.numLinkSites();link2_idx++) {
	        if (bondmap_.canLink(xlp, link1_idx, link2_idx)) {
            //create the candidate.
	          SelfLoopPeptide self_loop(xlp, xlp.getLinkSite(link1_idx), xlp.getLinkSite(link2_idx));
	          if (self_loop.getNumMissedCleavages() <= GlobalParams::getMissedCleavages()) {
	            self_loop.getMass(GlobalParams::getIsotopicMass());
              target_selfloop_peptides_.push_back(self_loop);
              added=true;
	          }
	        }
	      }
      }
    }
  }
  
  
  
  return(added);
}


void XLinkDatabase::initialize() {
  carp(CARP_INFO, "Initializing database");
  //Step one, load the database
  string input_file = Params::GetString("protein fasta file");
  string link_string = Params::GetString("link sites");
  bondmap_ = XLinkBondMap(link_string);
  protein_database_ = NULL;
  int num_protein = prepare_protein_input(input_file, &protein_database_);
    
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  //Step two create all peptides.
  int additional_cleavages = 0;

  bool generate_xlinkable = Params::GetBool("xlink-include-inter-intra") ||
           Params::GetBool("xlink-include-inter") ||
           Params::GetBool("xlink-include-intra");

  if (generate_xlinkable) {
    additional_cleavages = 1;
  }
  
  if (Params::GetBool("xlink-include-selfloops")) {
    additional_cleavages = 2;
  }
  
  //determine the maximum number of missed cleavages can be caused by modifications
  //This includes mono-links and normal mods.
  int blocked_cleavages = 0;
  
  vector<AA_MOD_T*> current_aa_mods;
  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    peptide_mod_get_aa_mods(peptide_mod, current_aa_mods);
    int current_blocked_cleavages = 0;
    for (vector<AA_MOD_T*>::iterator iter = current_aa_mods.begin();
         iter != current_aa_mods.end();
         ++iter) {
      if ((*iter)->getPreventsCleavage()) {
        current_blocked_cleavages++;
      }
    }
    blocked_cleavages = max(blocked_cleavages, current_blocked_cleavages);
  }
  
  int total_missed_cleavages = GlobalParams::getMissedCleavages()+additional_cleavages+blocked_cleavages;

  carp(CARP_INFO, "Missed cleavages from user:%d", GlobalParams::getMissedCleavages());
  carp(CARP_INFO, "Blocked cleavages from xlink/selfloops:%d", additional_cleavages);
  carp(CARP_INFO, "Blocked cleavages from mods/mono-links:%d", blocked_cleavages);
  carp(CARP_INFO, "Total missed cleavages to consider:%d", total_missed_cleavages);
  
  for (size_t idx=0;idx<=total_missed_cleavages;idx++) {
    target_peptides_.push_back(vector<Crux::Peptide*>());
    decoy_peptides_.push_back(vector<Crux::Peptide*>());
  }
  //Generate all possible target peptides
  
  size_t peptide_count = 0;
  size_t used_count = 0;
  
  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    double delta_mass = peptide_mod_get_mass_change(peptide_mod);
    carp(CARP_INFO, "Modification %d has delta mass of %g Da.", mod_idx + 1, delta_mass);

    ModifiedPeptidesIterator* peptide_iterator =
      new ModifiedPeptidesIterator(
        GlobalParams::getMinMass(), 
        GlobalParams::getMaxMass(), 
        peptide_mod, 
        false, 
        protein_database_,
        additional_cleavages+blocked_cleavages);

    //add the targets
    while (peptide_iterator->hasNext()) {
      peptide_count++;
      Crux::Peptide* peptide = peptide_iterator->next();
      int missed_cleavages = peptide->getMissedCleavageSites();
      
      bool added = addPeptideToDatabase(peptide);
      if (added) {
        used_count++;
        target_peptides_[missed_cleavages].push_back(peptide);
      } else {
        delete peptide;
      }
    }
    delete peptide_iterator;
  }
  
  carp(CARP_INFO,"Considered linear %d peptides, of which %d were included in the database.", peptide_count, used_count);
  
  if (Params::GetBool("xlink-include-linears")) {
    carp(CARP_INFO, "  The database contains %d linear peptides.", target_linear_peptides_.size());
    sort(target_linear_peptides_.begin(), target_linear_peptides_.end(), compareLinearPeptideMass);
  }
  
  if (GlobalParams::getXLinkIncludeDeadends()) {
    carp(CARP_INFO, "  The database contains %d mono-link peptides.", target_monolink_peptides_.size());
    sort(target_monolink_peptides_.begin(), target_monolink_peptides_.end(), compareMonoLinkPeptideMass);
  }
  
  if (Params::GetBool("xlink-include-selfloops")) {
    carp(CARP_INFO, "  The database contains %d selfloop peptides.", target_selfloop_peptides_.size());
    sort(target_selfloop_peptides_.begin(), target_selfloop_peptides_.end(), compareSelfLoopPeptideMass);
  }
  
  if (generate_xlinkable) {
    carp(CARP_INFO, "  The database contains %d cross-linkable peptides.", target_xlinkable_peptides_.size());
    sort(target_xlinkable_peptides_.begin(), target_xlinkable_peptides_.end(), compareXLinkablePeptideMass);
    flattenLinkablePeptides(target_xlinkable_peptides_, target_xlinkable_peptides_flatten_);
  }

  carp(CARP_INFO, "Done initializing database");
}

void XLinkDatabase::finalize() {
  carp(CARP_INFO, "Start to finalize database");
  target_linear_peptides_.clear();
  decoy_linear_peptides_.clear();
  target_selfloop_peptides_.clear();
  decoy_selfloop_peptides_.clear();
  target_xlinkable_peptides_.clear();
  decoy_xlinkable_peptides_.clear();
  target_xlinkable_peptides_flatten_.clear();
  for (size_t idx1=0;idx1<target_peptides_.size();idx1++) {
    for (size_t idx2=0;idx2<target_peptides_[idx1].size();idx2++) {
      delete target_peptides_[idx1][idx2];
    }
  }

  for (size_t idx1=0;idx1<decoy_peptides_.size();idx1++) {
    for (size_t idx2=0;idx2<decoy_peptides_[idx1].size();idx2++) {
      delete decoy_peptides_[idx1][idx2];
    }
  }

  Database::freeDatabase(protein_database_);
  carp(CARP_INFO, "Done finalizing database");
}

void XLinkDatabase::findSelfLoops(
  vector<XLinkablePeptide>& linkable_peptides, 
  vector<SelfLoopPeptide>& ans) {

  //Loop through all peptides, trying to find peptides with at least two link sites.
  for (vector<XLinkablePeptide>::iterator iter =
	 linkable_peptides.begin();
       iter != linkable_peptides.end();
       ++iter) {

    if (iter->numLinkSites() > 1) {
      for (size_t link1_idx = 0; link1_idx<iter->numLinkSites()-1;link1_idx++) {
	for (size_t link2_idx = link1_idx+1; link2_idx < iter->numLinkSites();link2_idx++) {
	  if (bondmap_.canLink(*iter, link1_idx, link2_idx)) {
            //create the candidate.
	    SelfLoopPeptide self_loop(*iter, iter->getLinkSite(link1_idx), iter->getLinkSite(link2_idx));
	    if (self_loop.getNumMissedCleavages() <= GlobalParams::getMissedCleavages()) {
	      self_loop.getMass(GlobalParams::getIsotopicMass());
              ans.push_back(self_loop);
	    }
	  }
	}
      }
    }
  }
}

void XLinkDatabase::generateAllSelfLoops(bool decoy) {

  int mc = GlobalParams::getMissedCleavages()+2;
  if (decoy) {
    findSelfLoops(decoy_xlinkable_peptides_, decoy_selfloop_peptides_);
    int mc = GlobalParams::getMissedCleavages()+2;
    generateAllLinkablePeptides(decoy_peptides_[mc], decoy_xlinkable_peptides2_);
    findSelfLoops(decoy_xlinkable_peptides2_, decoy_selfloop_peptides_);
    sort(decoy_selfloop_peptides_.begin(), 
      decoy_selfloop_peptides_.end(), 
      compareSelfLoopPeptideMass);
  } else {
    findSelfLoops(target_xlinkable_peptides_, target_selfloop_peptides_);
    generateAllLinkablePeptides(target_peptides_[mc], target_xlinkable_peptides2_);
    findSelfLoops(target_xlinkable_peptides2_, target_selfloop_peptides_);
    sort(
      target_selfloop_peptides_.begin(), 
      target_selfloop_peptides_.end(), 
      compareSelfLoopPeptideMass);
  }
}

void XLinkDatabase::findLinearPeptides(vector<Crux::Peptide*> &peptides, vector<LinearPeptide>& linear_peptides) {

  for (vector<Crux::Peptide*>::iterator iter = peptides.begin();
       iter != peptides.end();
       ++iter) {
    LinearPeptide lpeptide(*iter);
    lpeptide.getMass(GlobalParams::getIsotopicMass());
    linear_peptides.push_back(lpeptide);
    
  }
}

void XLinkDatabase::generateAllLinears(bool decoy) {
  //NOTE - Right now, we need to call this with decoy=false first.
  
  for (size_t cleavage_idx = 0; 
         cleavage_idx <= GlobalParams::getMissedCleavages(); 
         cleavage_idx++) {
      findLinearPeptides(target_peptides_[cleavage_idx], target_linear_peptides_);
    }
}


void XLinkDatabase::flattenLinkablePeptides(vector<XLinkablePeptide>& xpeptides,
					    vector<XLinkablePeptide>& flattened) {

  for (size_t idx=0;idx < xpeptides.size();idx++) {
    XLinkablePeptide& current = xpeptides[idx];
    for (size_t link1_idx=0;link1_idx < current.numLinkSites(); link1_idx++) {
      XLinkablePeptide onelink(current);
      onelink.getMass(MONO);
      onelink.clearSites();
      onelink.addLinkSite(current.getLinkSite(link1_idx));
      flattened.push_back(onelink);
    }
  }
}

void XLinkDatabase::filterLinkablePeptides(
  vector<XLinkablePeptide>& xpeptides,
  vector<XLinkablePeptide>& filtered_xpeptides
  ) {
  bool filter1 = !Params::GetBool("xlink-include-inter-intra");
  bool filter2 = !Params::GetBool("xlink-include-inter") && !Params::GetBool("xlink-include-intra");
  
  for (size_t idx = 0 ;idx < xpeptides.size();idx++) {
  //quick check to see if peptide can come from multiple protein sources.
  //if true, then this should return XLINK_INTER_INTRA_CANDIDATE
    XLINKMATCH_TYPE_T ctype = 
      XLink::getCrossLinkCandidateType(xpeptides[idx].getPeptide(), 
				xpeptides[idx].getPeptide());
    
    if ((filter1 && ctype == XLINK_INTER_INTRA_CANDIDATE) ||
      (filter2 && ctype != XLINK_INTER_INTRA_CANDIDATE)) {
    } else {
      filtered_xpeptides.push_back(xpeptides[idx]);
    }
  }

  carp(CARP_INFO, "kept %d out of %d", filtered_xpeptides.size(), xpeptides.size());

}

void XLinkDatabase::generateAllLinkablePeptides(
  vector<Crux::Peptide*>& peptides, 
  vector<XLinkablePeptide>& xpeptides) {

  carp(CARP_DEBUG, "XLinkDatabase::generateAllLinkablePeptides: start()");
  carp(CARP_DEBUG, "Number of peptides:%d", peptides.size());
  //Loop through peptides
  vector<int> link_sites;
  for (vector<Crux::Peptide*>::iterator iter = peptides.begin(); 
    iter != peptides.end(); 
       ++iter) {
    Crux::Peptide* peptide = *iter;

    if (peptide->countModifiedAAs() <= GlobalParams::getMaxXLinkMods()) {
      XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites); 
      if (!link_sites.empty()) {
	XLinkablePeptide xlp(peptide, link_sites);
	xlp.getMass(MONO);
        xpeptides.push_back(xlp);
      }
    }
  }
  carp(CARP_DEBUG, "XLinkDatabase::generateAllLinkablePeptides: done()");
}


void XLinkDatabase::print() {
   
    string output_directory = Params::GetString("output-dir");
    if (Params::GetBool("xlink-include-linears")) {
      ostringstream oss;
      oss << output_directory << "/" << "xlink_peptides.linear.txt";
      string temp = oss.str();
      ofstream peptides_file(temp.c_str());

      peptides_file << setprecision(8);

      peptides_file << "mass\tsequence\tprotein id\tmissed cleavages\tunshuffled sequence"<<endl;

      for (int idx=0;idx < target_linear_peptides_.size();idx++) {
     
        peptides_file << target_linear_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << target_linear_peptides_[idx].getSequenceString() << "\t";
        peptides_file << target_linear_peptides_[idx].getProteinIdString() << "\t";
        peptides_file << target_linear_peptides_[idx].getNumMissedCleavages() << "\t";
        peptides_file << "" << endl;
      }
      
      for (int idx=0;idx < decoy_linear_peptides_.size();idx++) {
        peptides_file << decoy_linear_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << decoy_linear_peptides_[idx].getSequenceString() << "\t";
        peptides_file << decoy_linear_peptides_[idx].getProteinIdString() << "\t";
        peptides_file << decoy_linear_peptides_[idx].getNumMissedCleavages() << "\t";
        peptides_file << decoy_linear_peptides_[idx].getPeptide(0)->getUnshuffledSequence() << endl;
      }
      
      peptides_file.flush();
    }
    if (target_monolink_peptides_.size()) {
      ostringstream oss;
      oss << output_directory << "/" << "xlink_peptides.monolink.txt";
      string temp = oss.str();
      ofstream peptides_file(temp.c_str());

      peptides_file << setprecision(8);

      peptides_file << "mass\tsequence\tprotein id\tmissed cleavages\tunshuffled sequence"<<endl;

      for (int idx=0;idx < target_monolink_peptides_.size();idx++) {
     
        peptides_file << target_monolink_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << target_monolink_peptides_[idx].getSequenceString() << "\t";
        peptides_file << target_monolink_peptides_[idx].getProteinIdString() << "\t";
        peptides_file << target_monolink_peptides_[idx].getNumMissedCleavages() << "\t";
        peptides_file << "" << endl;
      }
      peptides_file.flush();
      
    }
    if (Params::GetBool("xlink-include-selfloops")) {
      ostringstream oss;
      oss << output_directory << "/" << "xlink_peptides.selfloops.txt";
      string temp = oss.str();
      ofstream peptides_file(temp.c_str());

      peptides_file << setprecision(8);

      peptides_file << "mass\tsequence\tprotein id\t missed cleavages"<<endl;

      for (int idx=0;idx < target_selfloop_peptides_.size();idx++) {
     
        peptides_file << target_selfloop_peptides_[idx].getMass(MONO) << "\t";
        peptides_file << target_selfloop_peptides_[idx].getSequenceString() << "\t";
        peptides_file << target_selfloop_peptides_[idx].getNumMissedCleavages() << "\t";
        peptides_file << target_selfloop_peptides_[idx].getProteinIdString() << endl;
      }
      peptides_file.flush();
    }      
    
    ostringstream oss;
    oss << output_directory << "/" << "xlink_peptides.linkable.txt";
     
    string temp = oss.str();
    ofstream peptides_file(temp.c_str());

    peptides_file << setprecision(8);

    peptides_file << "mass\tsequence\tlinks"<<endl;

    for (int idx=0;idx < target_xlinkable_peptides_.size();idx++) {
     
      peptides_file << target_xlinkable_peptides_[idx].getMass(MONO) << "\t";
      peptides_file << target_xlinkable_peptides_[idx].getModifiedSequenceString() << "\t";
      
      peptides_file << target_xlinkable_peptides_[idx].getLinkSite(0)+1;
      for (size_t idx2 = 1;idx2 < target_xlinkable_peptides_[idx].numLinkSites();idx2++) {
        peptides_file << "," << target_xlinkable_peptides_[idx].getLinkSite(idx2)+1;
      }
      peptides_file << endl;
    }
    peptides_file.flush();
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearBegin(
  bool decoy
  ) {

  if (decoy) {
    return (decoy_linear_peptides_.begin());
  } else {
    return (target_linear_peptides_.begin());
  }
}
vector<LinearPeptide>::iterator XLinkDatabase::getLinearBegin(
  bool decoy,
  FLOAT_T min_mass
  ) {

  if (decoy) {
    return (lower_bound(decoy_linear_peptides_.begin(),
		        decoy_linear_peptides_.end(),
                        min_mass,
                        compareLinearPeptideMassToFLOAT));
  } else {
    return (lower_bound(target_linear_peptides_.begin(), 
                        target_linear_peptides_.end(), 
                        min_mass, 
                        compareLinearPeptideMassToFLOAT)); 
  }
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearEnd(
  bool decoy
) {

  if (decoy) {
    return (decoy_linear_peptides_.end());
  } else {
    return (target_linear_peptides_.end());
  }
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearEnd(
  bool decoy,
  FLOAT_T max_mass
  ) {
  return(std::upper_bound(getLinearBegin(decoy),
                          getLinearEnd(decoy),
                          max_mass,
                          compareFLOATToLinearPeptideMass));  

 
}

vector<LinearPeptide>::iterator XLinkDatabase::getLinearEnd(
  bool decoy,
  std::vector<LinearPeptide>::iterator& siter,
  FLOAT_T max_mass
) {
   return(std::upper_bound(siter,
                            getLinearEnd(decoy),
                            max_mass,
                            compareFLOATToLinearPeptideMass));  
  
}


vector<MonoLinkPeptide>::iterator XLinkDatabase::getMonoLinkBegin(
  bool decoy,
  FLOAT_T min_mass
  ) {
  return (lower_bound(target_monolink_peptides_.begin(), 
                      target_monolink_peptides_.end(), 
                      min_mass, 
                      compareMonoLinkPeptideMassToFLOAT)); 
}

vector<MonoLinkPeptide>::iterator XLinkDatabase::getMonoLinkBegin(
  bool decoy
) {
  return (target_monolink_peptides_.begin());
}

vector<MonoLinkPeptide>::iterator XLinkDatabase::getMonoLinkEnd(
  bool decoy
) {

  return (target_monolink_peptides_.end());
}

vector<MonoLinkPeptide>::iterator XLinkDatabase::getMonoLinkEnd(
  bool decoy,
  FLOAT_T max_mass
  ) {
  return(std::upper_bound(getMonoLinkBegin(decoy),
                          getMonoLinkEnd(decoy),
                          max_mass,
                          compareFLOATToMonoLinkPeptideMass));
}

vector<MonoLinkPeptide>::iterator XLinkDatabase::getMonoLinkEnd(
  bool decoy,
  std::vector<MonoLinkPeptide>::iterator& siter,
  FLOAT_T max_mass
) {
   return(std::upper_bound(siter,
                            getMonoLinkEnd(decoy),
                            max_mass,
                            compareFLOATToMonoLinkPeptideMass));  
}



vector<SelfLoopPeptide>::iterator XLinkDatabase::getSelfLoopBegin(
  bool decoy,
  FLOAT_T min_mass
  ) {

  if (decoy) {
    return (lower_bound(decoy_selfloop_peptides_.begin(), 
                        decoy_selfloop_peptides_.end(), 
                        min_mass, 
                        compareSelfLoopPeptideMassToFLOAT));
  } else {

    return (lower_bound(target_selfloop_peptides_.begin(), 
                        target_selfloop_peptides_.end(), 
                        min_mass, 
                        compareSelfLoopPeptideMassToFLOAT));
  }
}

vector<SelfLoopPeptide>::iterator XLinkDatabase::getSelfLoopEnd(
  bool decoy
) {

  if (decoy) {
    return(decoy_selfloop_peptides_.end());
  } else {
    return(target_selfloop_peptides_.end());
  }
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableBegin() {
  return(target_xlinkable_peptides_.begin());

}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableBegin(FLOAT_T min_mass) {
  return (lower_bound(target_xlinkable_peptides_.begin(), target_xlinkable_peptides_.end(),
		      min_mass, compareXLinkablePeptideMassToFLOAT));
}
vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableEnd() {
    return(target_xlinkable_peptides_.end());
}

vector<XLinkablePeptide>& XLinkDatabase::getXLinkablePeptides(
  bool decoy
) {
  if (decoy) {
    return (decoy_xlinkable_peptides_);
  } else {
    return(target_xlinkable_peptides_);
  }
}



XLinkBondMap& XLinkDatabase::getXLinkBondMap() {
  return bondmap_;
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenBegin() {
  return(target_xlinkable_peptides_flatten_.begin());
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenBegin(
  bool decoy,
  FLOAT_T min_mass
  ) {
  if (decoy) {
    return(lower_bound(decoy_xlinkable_peptides_flatten_.begin(), 
                       decoy_xlinkable_peptides_flatten_.end(),
	  	     min_mass, compareXLinkablePeptideMassToFLOAT));
  } else {
    return(lower_bound(target_xlinkable_peptides_flatten_.begin(), 
                       target_xlinkable_peptides_flatten_.end(),
	  	     min_mass, compareXLinkablePeptideMassToFLOAT));
  }
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenEnd() {
  return(target_xlinkable_peptides_flatten_.end());
}

vector<XLinkablePeptide>::iterator XLinkDatabase::getXLinkableFlattenEnd(
									 bool decoy,
  FLOAT_T max_mass
  ) {

  if (decoy) {
    return(std::upper_bound(decoy_xlinkable_peptides_flatten_.begin(),
	  	     decoy_xlinkable_peptides_flatten_.end(),
		     max_mass,
		     compareXLinkablePeptideMassToFLOAT2));
  } else {

    return(std::upper_bound(target_xlinkable_peptides_flatten_.begin(),
	  	     target_xlinkable_peptides_flatten_.end(),
		     max_mass,
		     compareXLinkablePeptideMassToFLOAT2));
  }
}


int XLinkDatabase::getNLinkable() {
  return(target_xlinkable_peptides_.size());
}
