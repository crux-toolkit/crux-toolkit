/**
 * \file XLinkPeptide.cpp 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September 2014
 * \brief Object for Defining a crosslinked peptide in an xlink search
 *****************************************************************************/
#include "XLinkPeptide.h"
#include "ModifiedPeptidesIterator.h"
#include "IonSeries.h"
#include "Ion.h"

#include "XLinkablePeptideIterator.h"


#include <iostream>
#include <sstream>

using namespace std;

FLOAT_T XLinkPeptide::linker_mass_ = 0;
set<Crux::Peptide*> XLinkPeptide::allocated_peptides_;
FLOAT_T XLinkPeptide::pmin_ = 0;
bool XLinkPeptide::pmin_set_ = false;

XLinkPeptide::XLinkPeptide() : XLinkMatch() {
  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
  is_decoy_ = false;
  
}

/**
 * Constructor using two linkable peptide objects and locations
 */
XLinkPeptide::XLinkPeptide(
  XLinkablePeptide& peptideA, ///< 1st peptide
  XLinkablePeptide& peptideB, ///< 2nd peptide
  int posA, ///<link pos1
  int posB ///<link pos2
  ) : XLinkMatch() {
  
  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
  
  linked_peptides_.push_back(peptideA);
  linked_peptides_.push_back(peptideB);
  link_pos_idx_.push_back(posA);
  link_pos_idx_.push_back(posB);

  doSort();

}

/**
 * Constructor for using two string peptides and locations
 */
XLinkPeptide::XLinkPeptide(
  char* peptideA, ///< sequence of peptide1
  char* peptideB, ///< sequence of peptide2
  int posA, ///< position of crosslink in peptide1
  int posB ///< position of crosslink in peptide2
  ) : XLinkMatch() {

  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
  XLinkablePeptide A(peptideA);
  linked_peptides_.push_back(A);
  XLinkablePeptide B(peptideB);
  linked_peptides_.push_back(B);
  A.addLinkSite(posA);
  link_pos_idx_.push_back(0);
  B.addLinkSite(posB);
  link_pos_idx_.push_back(0);

  doSort();
}


/**
 * makes sure that sequence1 is smaller in alphanumeric value than
 * sequence 2
 */
void XLinkPeptide::doSort() {

  string seq1 = linked_peptides_[0].getModifiedSequenceString();
  
  string seq2 = linked_peptides_[1].getModifiedSequenceString();


  if (seq1 > seq2) {

    //swap peptides
    swap(linked_peptides_[0], linked_peptides_[1]);
    //swap links
    swap(link_pos_idx_[0], link_pos_idx_[1]);
  }

  seq1 = linked_peptides_[0].getModifiedSequenceString();
  
  seq2 = linked_peptides_[1].getModifiedSequenceString();

  assert(seq1 <= seq2);



}

/**
 * sets the static linker mass variable
 */
void XLinkPeptide::setLinkerMass(
  FLOAT_T linker_mass ///< linker mass
  ) {
  linker_mass_=linker_mass;
}

/**
 * \returns the linker mass
 */
FLOAT_T XLinkPeptide::getLinkerMass() {
  return linker_mass_;
}

/**
 * \returns the link position within each peptide
 */
int XLinkPeptide::getLinkPos(
  int peptide_idx ///< 0 - first peptide, 1 - second peptide
  ) {
  
  return linked_peptides_[peptide_idx].getLinkSite(link_pos_idx_[peptide_idx]);
}

/**
 * \returns whether the cross-link is from peptides from two different
 * proteins
 */
bool XLinkPeptide::isInter() {

  Crux::Peptide* peptide_a = linked_peptides_[0].getPeptide();
  Crux::Peptide* peptide_b = linked_peptides_[1].getPeptide();

  vector<Crux::Protein*> proteins_a;
  vector<Crux::Protein*> proteins_b;


  for (PeptideSrcIterator src_iterator = peptide_a->getPeptideSrcBegin();
    src_iterator != peptide_a->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_a.push_back(src->getParentProtein());
  }

  for (PeptideSrcIterator src_iterator = peptide_b->getPeptideSrcBegin();
    src_iterator != peptide_b->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_b.push_back(src->getParentProtein());
  }
  
  sort(proteins_a.begin(), proteins_a.end());
  sort(proteins_b.begin(), proteins_b.end());


  for (unsigned int idx_a=0;idx_a<proteins_a.size();idx_a++) {

    for (unsigned int idx_b=0;idx_b<proteins_b.size();idx_b++) {
      
      string id_a = proteins_a[idx_a] -> getIdPointer();
      string id_b = proteins_b[idx_b] -> getIdPointer();

      if (id_a.find(id_b) == string::npos && id_b.find(id_a) == string::npos) {
        //Found an instance where the proteins are not equal, therefore it is
        //inter. but still could be intra....
        return true;
      }
    }
  }

  //all proteins equal, so it is definately intra
  return false;

}

/**
 * \returns whether the cross-link is from peptides within the same protein
 */
bool XLinkPeptide::isIntra() {

  Crux::Peptide* peptide_a = linked_peptides_[0].getPeptide();
  Crux::Peptide* peptide_b = linked_peptides_[1].getPeptide();

  vector<Crux::Protein*> proteins_a;
  vector<Crux::Protein*> proteins_b;
  

  for (PeptideSrcIterator src_iterator = peptide_a->getPeptideSrcBegin();
    src_iterator != peptide_a->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_a.push_back(src->getParentProtein());
  }

  for (PeptideSrcIterator src_iterator = peptide_b->getPeptideSrcBegin();
    src_iterator != peptide_b->getPeptideSrcEnd();
    ++src_iterator) {
    PeptideSrc* src = *src_iterator;
    proteins_b.push_back(src->getParentProtein());
  }

  sort(proteins_a.begin(), proteins_a.end());
  sort(proteins_b.begin(), proteins_b.end());


  for (unsigned int idx_a=0;idx_a<proteins_a.size();idx_a++) {

    for (unsigned int idx_b=0;idx_b<proteins_b.size();idx_b++) {
      string id_a = proteins_a[idx_a] -> getIdPointer();
      string id_b = proteins_b[idx_b] -> getIdPointer();
      if (id_a.find(id_b) != string::npos || id_b.find(id_a) != string::npos) {
        //Found an instance where the protein are equal, therefore it is
        //intra. but still also be inter as well....
        return true;
      }
    }
  }

  //all proteins inequal, so it is definately inter
  return false;
}

/**
 * Gets all peptides that are linkable, i.e. have link sites
 */
void XLinkPeptide::addLinkablePeptides(
  double min_mass, ///< min mass of peptides
  double max_mass, ///< max mass of peptides
  Index* index, ///< protein index
  Database* database, ///< protein database
  PEPTIDE_MOD_T* peptide_mod, ///< modifications
  bool is_decoy,  ///< are the peptides decoys
  XLinkBondMap& bondmap, ///< valid crosslink map
  vector<XLinkablePeptide>& linkable_peptides ///< list of linkable peptides -out
  ) {

  ModifiedPeptidesIterator* peptide_iterator =
    new ModifiedPeptidesIterator(
      min_mass, 
      max_mass,
      peptide_mod, 
      is_decoy,
      index, 
      database);

  int max_mod_xlink = get_int_parameter("max-xlink-mods");

  while (peptide_iterator->hasNext()) {
    Crux::Peptide* peptide = peptide_iterator->next();
    vector<int> link_sites;

    if (peptide->countModifiedAAs() > max_mod_xlink) {
      delete peptide;
    } else {
      XLinkablePeptide::findLinkSites(peptide, bondmap, link_sites);

      if (link_sites.size() > 0) {
        XLinkablePeptide xlinkable_peptide(peptide, link_sites);
        xlinkable_peptide.setDecoy(is_decoy);
        linkable_peptides.push_back(xlinkable_peptide);
        XLink::addAllocatedPeptide(peptide);
      } else {
        delete peptide;
      }
    }
  }
  
  delete peptide_iterator;
}

/***
 * adds crosslink candidates by iterating through all possible masses
 */
void XLinkPeptide::addCandidates(
  FLOAT_T min_mass, ///< min mass of crosslink
  FLOAT_T max_mass, ///< max mass of crosslinks
  XLinkBondMap& bondmap, ///< valid crosslink map
  Index* index,  ///< protein index
  Database* database, ///< protein database
  PEPTIDE_MOD_T** peptide_mods, ///< modifications for the peptides
  int num_peptide_mods, ///< number of possible modifications
  XLinkMatchCollection& candidates ///< candidates in/out
  ) {

  if (!pmin_set_) {
    FLOAT_T min_length_mass = get_mass_amino_acid('G', MONO) * 
      (FLOAT_T)get_int_parameter("min-length") + 
      MASS_H2O_MONO;
    pmin_ = max((FLOAT_T)get_double_parameter("min-mass"), min_length_mass);
    pmin_set_ = true;
  }
  FLOAT_T peptide1_min_mass = pmin_;
  FLOAT_T peptide1_max_mass = max_mass-pmin_-linker_mass_;

  for (int mod_idx1 = 0; mod_idx1 < num_peptide_mods; mod_idx1++) {
    PEPTIDE_MOD_T* peptide_mod1 = peptide_mods[mod_idx1];
    XLinkablePeptideIterator iter1_target(peptide1_min_mass, peptide1_max_mass, index, database, peptide_mod1, false, bondmap);
    for (int mod_idx2 = 0; mod_idx2 < num_peptide_mods; mod_idx2++) {
      PEPTIDE_MOD_T* peptide_mod2 = peptide_mods[mod_idx2];
      addCandidates(min_mass, max_mass, bondmap, index, database, peptide_mod2, false, iter1_target, candidates);
    }
  } 

}

/**
 * adds crosslink candidates to the XLinkMatchCollection using
 * the passed in iterator for the 1st peptide
 */
void XLinkPeptide::addCandidates(
  FLOAT_T min_mass, ///< min mass of crosslinks
  FLOAT_T max_mass, ///< max mass of crosslinks
  XLinkBondMap& bondmap, ///< valid crosslink map
  Index* index, ///< protein index
  Database* database, ///< protein database
  PEPTIDE_MOD_T* peptide_mod2, ///<modification of the modified peptides
  bool decoy2, ///< are these going to be decoys?
  XLinkablePeptideIterator& iter1, ///< 1st peptide iterator
  XLinkMatchCollection& candidates ///< candidates -in/out
  ) {

  int max_mod_xlink = get_int_parameter("max-xlink-mods");

  set<XLinkablePeptide> visited;

  while (iter1.hasNext()) {
    XLinkablePeptide pep1 = iter1.next();
    //cerr <<"pep1:"<<pep1.getSequence()<<endl;
    char* seq1 = pep1.getPeptide()->getUnshuffledSequence();
    string seq1_string(seq1);
    std::free(seq1);

    FLOAT_T peptide2_min_mass = min_mass - pep1.getMass() - linker_mass_;
    FLOAT_T peptide2_max_mass = max_mass - pep1.getMass() - linker_mass_;

    assert (peptide2_min_mass <= peptide2_max_mass);
  

    if (peptide2_max_mass >= pmin_) {
      XLinkablePeptideIterator iter2(peptide2_min_mass, peptide2_max_mass, index, database, peptide_mod2, decoy2, bondmap);
      while (iter2.hasNext()) {
        XLinkablePeptide pep2 = iter2.next();
        if (visited.find(pep2) == visited.end()) {
          FLOAT_T mass = pep1.getMass() + pep2.getMass() + linker_mass_;
    
          if ((mass >= min_mass) && (mass <= max_mass)) {
    
            int mods = pep1.getPeptide()->countModifiedAAs() + pep2.getPeptide()->countModifiedAAs();
            if (mods <= max_mod_xlink) {
    
              //for every linkable site, generate the candidate if it is legal.
              for (unsigned int link1_idx=0;link1_idx < pep1.numLinkSites(); link1_idx++) {
                for (unsigned int link2_idx=0;link2_idx < pep2.numLinkSites();link2_idx++) {
                  if (bondmap.canLink(pep1, pep2, link1_idx, link2_idx)) {
                  //create the candidate
                    XLinkMatch* newCandidate = 
                      new XLinkPeptide(pep1, pep2, link1_idx, link2_idx);
                    candidates.add(newCandidate);
                  }
                }
              } // for link1_idx 
            } // if (mods <= max_mod_xlink .. 
          } // if ((mass >= min_mass) && (mass <= max_mass))
   
        } // if (visited.find(pep2) == visited.end()
  
      } // while (iter2.hasNext()) 
  
    } // if (peptide2_max_mass >= pmin_)

    visited.insert(pep1);
  } // while (iter1.hasNext()) 
}  



/**
 * \returns the candidate type
 */
XLINKMATCH_TYPE_T XLinkPeptide::getCandidateType() {

  if (isInter()) {
    if (isIntra()) {
      return XLINK_INTER_INTRA_CANDIDATE;
    } else {
      return XLINK_INTER_CANDIDATE;
    }
  } else if (isIntra()) {
    return XLINK_INTRA_CANDIDATE;
  } else {
    carp(CARP_ERROR, "Something wrong happened!");
  }
  return XLINK_INTER_INTRA_CANDIDATE;
}

/**
 * \returns the sequence string
 */
string XLinkPeptide::getSequenceString() {

  doSort();

  string seq1 = linked_peptides_[0].getModifiedSequenceString();
  
  string seq2 = linked_peptides_[1].getModifiedSequenceString();




  //assert(seq1 <= seq2);

  ostringstream oss;
  oss << seq1 << ", " << 
    seq2 << " (" <<
    (getLinkPos(0)+1) << "," <<
    (getLinkPos(1)+1) << ")";

  string svalue = oss.str();

  return svalue;
}

/**
 * \returns the mass of the xlink peptide
 */
FLOAT_T XLinkPeptide::calcMass(MASS_TYPE_T mass_type) {
  return linked_peptides_[0].getMass(mass_type) + 
    linked_peptides_[1].getMass(mass_type) + 
    linker_mass_;
}

/**
 * \returns a shuffled xlink peptide
 */
XLinkMatch* XLinkPeptide::shuffle() {

  XLinkPeptide* decoy = new XLinkPeptide();
  decoy->setZState(getZState());
  decoy->linked_peptides_.push_back(linked_peptides_[0].shuffle());
  decoy->linked_peptides_.push_back(linked_peptides_[1].shuffle());
  decoy->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy->link_pos_idx_.push_back(link_pos_idx_[1]);

  return (XLinkMatch*)decoy;


}

/**
 * fills the ion series with the predicted ions for the cross linked candidate
 */ 
void XLinkPeptide::predictIons(
  IonSeries* ion_series,  ///< IonSeries object to fill
  int charge, ///< charge state of candidate
  bool first ///< is this the first peptide?
  ) {
  carp(CARP_DEBUG, "predictIons:start");
  MASS_TYPE_T fragment_mass_type = get_mass_type_parameter("fragment-mass");

  char* seq = NULL;
   MODIFIED_AA_T* mod_seq = NULL;
  int link_pos;
  FLOAT_T mod_mass;

  if (first) {
    carp(CARP_DEBUG, "predicting first peptide");
    //predict the ion series from the first peptide
    seq = linked_peptides_[0].getSequence();
    mod_seq = linked_peptides_[0].getModifiedSequence(); 
    link_pos = getLinkPos(0);
    mod_mass = linked_peptides_[1].getMass(fragment_mass_type) + linker_mass_;
  } else {
    carp(CARP_DEBUG, "predicting second peptide"); 
    //predict the ion series for the second peptide 
    seq = linked_peptides_[1].getSequence();
    mod_seq = linked_peptides_[1].getModifiedSequence();
    link_pos = getLinkPos(1);
    mod_mass = linked_peptides_[0].getMass(fragment_mass_type) + linker_mass_;
  }
  carp(CARP_DEBUG, "predicting ions");
  //predict the ion series of the peptide
  ion_series->setCharge(charge); 
  ion_series->update(seq, mod_seq); 
  ion_series->predictIons(); 
 
  carp(CARP_DEBUG, "modifying ions");
  //modify the necessary ions and add to the ion_series   
  for (IonIterator ion_iter = ion_series->begin(); 
    ion_iter != ion_series->end(); 
    ++ion_iter) { 
 
    Ion* ion = *ion_iter; 
 
    unsigned int cleavage_idx = ion->getCleavageIdx(); 
    if (ion->isForwardType()) { 
      if (cleavage_idx > (unsigned int)link_pos) {
        FLOAT_T mass = ion->getMassFromMassZ() + mod_mass;
        ion->setMassZFromMass(mass); 
        if (isnan(ion->getMassZ())) { 
          carp(CARP_FATAL, "NAN3"); 
        } 
      } 
    } else { 
      if (cleavage_idx >= (strlen(seq)-(unsigned int)link_pos)) { 
        FLOAT_T mass = ion->getMassFromMassZ() + mod_mass;
        ion->setMassZFromMass(mass); 
        if (isnan(ion->getMassZ())) { 
          carp(CARP_FATAL, "NAN4"); 
        } 
      } 
    } 
  } 

  free(seq); 
  free(mod_seq);
} 

/**
 * predicts the ions for xlinked peptide
 */
void XLinkPeptide::predictIons(
  IonSeries* ion_series, ///< IonSeries to fill
  int charge ///< charge state of the peptide
  ) {

  MASS_TYPE_T fragment_mass_type = get_mass_type_parameter("fragment-mass"); 

  //predict the ion_series of the first peptide.
  char* seq1 = linked_peptides_[0].getSequence();
  MODIFIED_AA_T* mod_seq1 = linked_peptides_[0].getModifiedSequence();

  ion_series->setCharge(charge);
  ion_series->update(seq1, mod_seq1);
  ion_series->predictIons();

  //iterate through all of the ions, if the ion contains a link, then
  //add the mass of peptide2 + linker_mass.

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;

    unsigned int cleavage_idx = ion->getCleavageIdx();

    if (ion->isForwardType()) {
      if (cleavage_idx > (unsigned int)getLinkPos(0)) {
        FLOAT_T mass = ion->getMassFromMassZ();
        mass += linked_peptides_[1].getMass(fragment_mass_type) + linker_mass_;
        ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN1");
        }
      }
    } else {
      if (cleavage_idx >= (strlen(seq1) - (unsigned int)getLinkPos(0))) {
        FLOAT_T mass = ion->getMassFromMassZ();
        mass += linked_peptides_[1].getMass(fragment_mass_type) + linker_mass_;
        ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN2");
        }
      }
    }
  }

  //predict the ion_series of the second peptide.
  IonConstraint* ion_constraint = ion_series->getIonConstraint();

  IonSeries* ion_series2 = 
      new IonSeries(ion_constraint, charge);
  
  char* seq2 = linked_peptides_[1].getSequence();


  MODIFIED_AA_T* mod_seq2 = 
    linked_peptides_[1].getModifiedSequence();
  ion_series2->setCharge(charge);
  ion_series2->update(seq2, mod_seq2);
  ion_series2->predictIons();

  //modify the necessary ions and add to the ion_series  
  for (IonIterator ion_iter = ion_series2->begin();
    ion_iter != ion_series2->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;

    unsigned int cleavage_idx = ion->getCleavageIdx();
    if (ion->isForwardType()) {
      if (cleavage_idx > (unsigned int)getLinkPos(1)) {
       FLOAT_T mass = ion->getMassFromMassZ();
       mass += linked_peptides_[0].getMass(fragment_mass_type) + linker_mass_;
       ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN3");
        }
      }
    } else {
      if (cleavage_idx >= (strlen(seq2)-(unsigned int)getLinkPos(1))) {
       FLOAT_T mass = ion->getMassFromMassZ();
       mass += linked_peptides_[0].getMass(fragment_mass_type) + linker_mass_;
       ion->setMassZFromMass(mass);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN4");
        }
      }
    }
    ion_series->addIon(ion);
  }
  free(seq1);
  free(seq2);
  free(mod_seq1);
  free(mod_seq2);
  
  delete ion_series2;
}

/**
 * \returns the sequence from the ion
 */
string XLinkPeptide::getIonSequence(
  Ion* ion ///< pointer to the ion
  ) {

  int peptide_idx = 0;

  string ion_sequence = ion->getPeptideSequence();

  if (ion_sequence == linked_peptides_[0].getSequence()) {
    peptide_idx = 0;
  } else {
    peptide_idx = 1;
  }

  unsigned int cleavage_idx = ion->getCleavageIdx();

  bool is_linked = false;

  if (ion->isForwardType()) {
    is_linked = (cleavage_idx > (unsigned int)getLinkPos(peptide_idx)); 
  } else {
    is_linked = (cleavage_idx >= (ion_sequence.length() - getLinkPos(peptide_idx)));
  }

  string subseq;
  if (ion->isForwardType()) {
    subseq = ion_sequence.substr(0, cleavage_idx);
  } else {
    subseq = ion_sequence.substr(ion_sequence.length() - cleavage_idx, ion_sequence.length());
  }

  if (!is_linked) {
    return subseq;
  } else {
    string ans;
    if (peptide_idx == 0) {
      char* seq2 = linked_peptides_[1].getSequence();
      ans = subseq + string(",") + string(seq2);
      free(seq2);
    } else {
      char* seq1 = linked_peptides_[0].getSequence();
      ans = string(seq1) + string(",") + subseq;
      free(seq1);
    }
    return ans;
  }
}

/**
 * \returns the Peptide object for the xlinked peptide
 */
Crux::Peptide* XLinkPeptide::getPeptide(
  int peptide_idx ///< 0 or 1
  ) {
  return linked_peptides_[peptide_idx].getPeptide();
}

/**
 * \returns the number of missed cleavages for the cross-linked peptide
 */
int XLinkPeptide::getNumMissedCleavages() {

  char missed_cleavage_link_site = 'K';
  set<int> skip;

  int link1_site = getLinkPos(0);
  int link2_site = getLinkPos(1);
  
  Crux::Peptide* pep1 = linked_peptides_[0].getPeptide();
  Crux::Peptide* pep2 = linked_peptides_[1].getPeptide();
  
  char *seq1 = pep1->getSequencePointer();
  char *seq2 = pep2->getSequencePointer();

  if (seq1[link1_site] == missed_cleavage_link_site) {
    skip.insert(link1_site);
  }

  int missed1 = pep1->getMissedCleavageSites(skip);
  
  skip.clear();

  if (seq2[link2_site] == missed_cleavage_link_site) {
    skip.insert(link2_site);
  }
  
  int missed2 = pep2->getMissedCleavageSites(skip);

  return max(missed1, missed2);

}

/**
 *\returns whether the cross-linked peptide is modified
 */
bool XLinkPeptide::isModified() {

  return linked_peptides_[0].isModified() || linked_peptides_[1].isModified();
}

/**
 *\returns the protein id string for the xlinked peptide
 */
string XLinkPeptide::getProteinIdString() {

  doSort();

  ostringstream oss;

  Crux::Peptide* peptide = this -> getPeptide(0);

  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide : Null first peptide!");
  } else {
    oss << peptide->getProteinIdsLocations();
  }
  oss << ";";

  peptide = this -> getPeptide(1);

  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide : Null second peptide!");
  } else {
    oss << peptide->getProteinIdsLocations();
  }

  return oss.str();
}

/**
 * \returns the protein id strings where the (X) is the position
 * within the protein
 */
string XLinkPeptide::getProteinIdsXLocations(
  int idx ///< peptide index (0 or 1)
  ) {

  Crux::Peptide* peptide = this ->getPeptide(idx);

  ostringstream oss;

  set<string> xlocations;

  int link_pos = getLinkPos(idx);

  for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
    iter != peptide->getPeptideSrcEnd();
    ++iter) {
    PeptideSrc* peptide_src = *iter;
    Crux::Protein* protein = peptide_src->getParentProtein();
    string protein_id = protein->getIdPointer();
    int peptide_loc = peptide_src->getStartIdx();
    ostringstream protein_loc_stream;
    protein_loc_stream << protein_id << "(" << (peptide_loc + link_pos) << ")";
    xlocations.insert(protein_loc_stream.str());
  }

  set<string>::iterator result_iter = xlocations.begin();
  string result_string = *result_iter;

  while (++result_iter != xlocations.end()) {
    result_string += "," + *result_iter;
  }

  return result_string;
}

/**
 * \returns the protein id string where the (X)s are the positions in
 * the proteins
 */
string XLinkPeptide::getProteinIdXString() {
  doSort();

  ostringstream oss;

  oss << getProteinIdsXLocations(0);
  oss << ";";
  oss << getProteinIdsXLocations(1);
  return oss.str();

}

/**
 * \returns the flanking amino acids for both peptides, separated by ;
 */
string XLinkPeptide::getFlankingAAString() {

  doSort();

  ostringstream oss;

  Crux::Peptide* peptide = this -> getPeptide(0);
  
  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide::getFlankingAAString() : Null first peptide!");
  } else {

    char* flanking_aas = peptide->getFlankingAAs();
    oss << flanking_aas;
    std::free(flanking_aas);
  }

  oss << ";";

  peptide = this->getPeptide(1);

  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide::getFlankingAAString() : Null second peptide!");
  } else {

    char* flanking_aas = peptide->getFlankingAAs();
    oss << flanking_aas;
    std::free(flanking_aas);
  }
  
  return oss.str();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
