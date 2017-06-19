/**
 * \file XLinkPeptide.cpp 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September 2014
 * \brief Object for Defining a crosslinked peptide in an xlink search
 *****************************************************************************/
#include "XLinkPeptide.h"
#include "model/ModifiedPeptidesIterator.h"
#include "model/IonSeries.h"
#include "model/Ion.h"
#include "util/Params.h"
#include "util/GlobalParams.h"
#include "XLinkDatabase.h"
#include "XLinkablePeptideIterator.h"
#include "XLinkablePeptideIteratorTopN.h"

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
  target_ = NULL;
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
  getMass(GlobalParams::getIsotopicMass());
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
  A.addLinkSite(posA);
  linked_peptides_.push_back(A);

  XLinkablePeptide B(peptideB);
  B.addLinkSite(posB);
  linked_peptides_.push_back(B);

  link_pos_idx_.push_back(0);
  link_pos_idx_.push_back(0);
  doSort();
  getMass(GlobalParams::getIsotopicMass());
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
  linker_mass_ = linker_mass;
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
  
  return linked_peptides_.at(peptide_idx).getLinkSite(link_pos_idx_.at(peptide_idx));
}

int XLinkPeptide::getLinkIdx(
			     int peptide_idx ///< 0 - first peptide, 1 -second peptide
			     ) {
  return link_pos_idx_[peptide_idx];
}

bool XLinkPeptide::isDecoy() {

  return linked_peptides_.at(0).isDecoy() || linked_peptides_.at(1).isDecoy();

}

string XLinkPeptide::getDecoyType() {
  string returnValue = "target-target";
  
  if (linked_peptides_.at(0).isDecoy() && linked_peptides_.at(1).isDecoy()) {
    returnValue = "decoy-decoy";
  } else if (!linked_peptides_.at(0).isDecoy() && linked_peptides_.at(1).isDecoy()) {
    returnValue = "target-decoy";
  } else if (linked_peptides_.at(0).isDecoy() && !linked_peptides_.at(1).isDecoy()) {
    returnValue = "decoy-target";
  }
  return(returnValue);
}

/**
 * \returns whether the cross-link is from peptides from two different
 * proteins
 */
bool XLinkPeptide::isInter() {
  return XLink::isCrossLinkInter(linked_peptides_.at(0).getPeptide(), linked_peptides_.at(1).getPeptide());
}

/**
 * \returns whether the cross-link is from peptides within the same protein
 */
bool XLinkPeptide::isIntra() {
  return XLink::isCrossLinkInter(linked_peptides_[0].getPeptide(), linked_peptides_[1].getPeptide());
}

/**
 * \returns whether the cross-link is from peptides that are both within the same protein and from different proteins
 */
bool XLinkPeptide::isInterIntra() {
  return(isInter() && isIntra());
}



/***
 * adds crosslink candidates by iterating through all possible masses
 */
int XLinkPeptide::addCandidates(
  Crux::Spectrum* spectrum,
  FLOAT_T precursor_mass,
  int precursor_charge,
  FLOAT_T min_mass, ///< min mass of crosslink
  FLOAT_T max_mass, ///< max mass of crosslinks
  bool decoy,
  XLinkMatchCollection& candidates ///< candidates in/out
  ) {

  carp(CARP_DEBUG, "XLinkPeptide::addCandidates - precursor:%g", precursor_mass);
  carp(CARP_DEBUG, "XLinkPeptide::addCandidates - min:%g", min_mass);
  carp(CARP_DEBUG, "XLinkPeptide::addCandidates - max:%g", max_mass);

  if (!pmin_set_) {
    pmin_ = XLinkDatabase::getXLinkableBegin()->getMass(GlobalParams::getIsotopicMass());
    pmin_set_ = true;
  }
  FLOAT_T peptide1_min_mass = pmin_;
  FLOAT_T peptide1_max_mass = max_mass-pmin_-linker_mass_;

  carp(CARP_DEBUG, "peptide1_min:%g", peptide1_min_mass);
  carp(CARP_DEBUG, "peptide1_max:%g", peptide1_max_mass);

  if (GlobalParams::getXLinkTopN() > 0) {
    vector<XLinkablePeptide> xlinkable_peptides;
    XLinkablePeptideIteratorTopN xlp_iter(spectrum, 
					  precursor_mass, 
					  peptide1_min_mass, 
					  peptide1_max_mass, 
					  precursor_charge, 
					  decoy);
    while(xlp_iter.hasNext()) {
      xlinkable_peptides.push_back(xlp_iter.next());
    }
    sort(xlinkable_peptides.begin(), xlinkable_peptides.end(), compareXLinkablePeptideMass);
    carp(CARP_DEBUG, "get xcorr");
    for (size_t idx =0;idx<xlinkable_peptides.size();idx++) {
      carp(CARP_DEBUG, "%f", xlinkable_peptides[idx].getXCorr());
    }
    return(addCandidates(min_mass, max_mass, xlinkable_peptides, candidates));
  } else {
    return(addCandidates(min_mass, max_mass, XLinkDatabase::getXLinkablePeptides(decoy), candidates));
  }
}

int XLinkPeptide::addXLinkPeptides(
  XLinkablePeptide& pep1, 
  XLinkablePeptide& pep2,
  XLinkMatchCollection& candidates
  ) {

  XLinkBondMap& bondmap = XLinkDatabase::getXLinkBondMap();
  int num_candidates = 0;
  //for every linkable site, generate the candidate if it is legal.
  for (unsigned int link1_idx=0;link1_idx < pep1.numLinkSites(); link1_idx++) {
    for (unsigned int link2_idx=0;link2_idx < pep2.numLinkSites();link2_idx++) {
      if (bondmap.canLink(pep1, pep2, link1_idx, link2_idx)) {
        //create the candidate
        XLinkMatch* newCandidate = 
          new XLinkPeptide(pep1, pep2, link1_idx, link2_idx);
	if (newCandidate->getNumMissedCleavages() <= GlobalParams::getMissedCleavages()) {
          candidates.add(newCandidate);
          num_candidates++;
	} else {
	  delete newCandidate;
	}
      }
    }
  }
  return(num_candidates);
}

/**
 * adds crosslink candidates to the XLinkMatchCollection using
 * the passed in iterator for the 1st peptide
 */
int XLinkPeptide::addCandidates(
  FLOAT_T min_mass, ///< min mass of crosslinks
  FLOAT_T max_mass, ///< max mass of crosslinks
  vector<XLinkablePeptide>& linkable_peptides, 
  XLinkMatchCollection& candidates ///< candidates -in/out
  ) {

  bool include_inter = GlobalParams::getXLinkIncludeInter();
  bool include_intra = GlobalParams::getXLinkIncludeIntra();
  bool include_inter_intra = GlobalParams::getXLinkIncludeInterIntra();

  int max_mod_xlink = GlobalParams::getMaxXLinkMods();
  
  size_t xpeptide_count = linkable_peptides.size();
  if (xpeptide_count <= 0) { return 0;}

  int num_candidates = 0;
  
  bool done = false;

  for (size_t pep_idx1=0;pep_idx1 < xpeptide_count-1;pep_idx1++) {
    XLinkablePeptide& pep1 = linkable_peptides.at(pep_idx1);
    carp(CARP_DEBUG, "pep_idx1:%d %d %f %s",
         pep_idx1,
         xpeptide_count-1,
         linkable_peptides[pep_idx1].getMass(MONO),
         pep1.getModifiedSequenceString().c_str()
         );
    FLOAT_T pep1_mass = pep1.getMass(MONO);
    FLOAT_T pep2_min_mass = min_mass - pep1_mass - linker_mass_;
    FLOAT_T pep2_max_mass = max_mass - pep1_mass - linker_mass_;
    int start_idx2 = pep_idx1+1;
      
    if (pep1_mass + linker_mass_ + linkable_peptides[start_idx2].getMass(MONO) > max_mass) {
      break;
    }
    for (size_t pep_idx2=start_idx2;pep_idx2 < xpeptide_count;pep_idx2++) {
      
      XLinkablePeptide& pep2 = linkable_peptides[pep_idx2];
      carp(CARP_DEBUG, "pep_idx2:%d %d %f %s",
         pep_idx2,
         xpeptide_count,
         linkable_peptides[pep_idx2].getMass(MONO),
         pep2.getModifiedSequenceString().c_str()
         );
      FLOAT_T current_mass = pep2.getMass(MONO);
      if (current_mass > pep2_max_mass) {
	      if (pep_idx2 == start_idx2) {
	        //done = true;
	      }
	      break;
      }
      if (current_mass >= pep2_min_mass) {
        XLINKMATCH_TYPE_T ctype = 
          XLink::getCrossLinkCandidateType(pep1.getPeptide(), pep2.getPeptide());
            
        if ((include_intra && ctype == XLINK_INTRA_CANDIDATE) || 
            (include_inter_intra && ctype == XLINK_INTER_INTRA_CANDIDATE) ||
            (include_inter && ctype == XLINK_INTER_CANDIDATE)) {
		            carp(CARP_DEBUG, "considering %s %s", pep1.getModifiedSequenceString().c_str(), pep2.getModifiedSequenceString().c_str());
	
              int mods = pep1.getPeptide()->countModifiedAAs() + pep2.getPeptide()->countModifiedAAs();
              if (mods <= max_mod_xlink) {
		            carp(CARP_DEBUG, "considering2 %s %s", pep1.getModifiedSequenceString().c_str(), pep2.getModifiedSequenceString().c_str());
                num_candidates += addXLinkPeptides(pep1, pep2, candidates);
              } // if (mods <= max_mod_xlink ..     
          }
	      
        }
      }
      if (done) {
        break;
      }
    
  }
   
  carp(CARP_DEBUG, "Done searching");
//  delete []tested;
  return(num_candidates);
}
  



/**
 * \returns the candidate type
 */
XLINKMATCH_TYPE_T XLinkPeptide::getCandidateType() {

  return(XLink::getCrossLinkCandidateType(linked_peptides_.at(0).getPeptide(), linked_peptides_.at(1).getPeptide()));

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

string XLinkPeptide::getUnshuffledSequence() {
  
  if (is_decoy_) {
    if (target_ == NULL) {
      carp(CARP_FATAL, "null target?!?");
    }
    return(target_ -> getSequenceString());
  } else {
    return (getSequenceString());
  }
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

void XLinkPeptide::shuffle(vector<XLinkMatch*>& decoys) {

  //cerr<<"tt:"<<getSequenceString()<<endl;
  
  XLinkablePeptide d1 = linked_peptides_[0].shuffle();
  XLinkablePeptide d2 = linked_peptides_[1].shuffle();

  XLinkPeptide* decoy_ff = new XLinkPeptide();
  decoy_ff->linked_peptides_.push_back(d1);
  decoy_ff->linked_peptides_.push_back(d2);
  decoy_ff->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy_ff->link_pos_idx_.push_back(link_pos_idx_[1]);
  decoy_ff->is_decoy_ = true;
  decoy_ff->target_ = this;
  decoy_ff->setZState(getZState());
  
  //cerr<<"decoy_ff:"<<decoy_ff->getSequenceString()<<endl;
  
  XLinkPeptide* decoy_tf = new XLinkPeptide();
  decoy_tf->linked_peptides_.push_back(linked_peptides_[0]);
  decoy_tf->linked_peptides_.push_back(d2);
  decoy_tf->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy_tf->link_pos_idx_.push_back(link_pos_idx_[1]);
  decoy_tf->is_decoy_ = true;
  decoy_tf->target_ = this;
  decoy_tf->setZState(getZState());

  //cerr <<"decoy_tf:"<<decoy_tf->getSequenceString()<<endl;
  
  XLinkPeptide* decoy_ft = new XLinkPeptide();
  decoy_ft->linked_peptides_.push_back(d1);
  decoy_ft->linked_peptides_.push_back(linked_peptides_[1]);
  decoy_ft->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy_ft->link_pos_idx_.push_back(link_pos_idx_[1]);
  decoy_ft->is_decoy_ = true;
  decoy_ft->target_ = this;
  decoy_ft->setZState(getZState());

  //cerr <<"decoy_ft:"<<decoy_ft->getSequenceString()<<endl;

  for (size_t idx=0;idx<NUMBER_MASS_TYPES;idx++) {
    if (mass_calculated_[idx]) {
      decoy_ff->mass_calculated_[idx] = true;
      decoy_ff->mass_[idx] = mass_[idx];
      decoy_ft->mass_calculated_[idx] = true;
      decoy_ft->mass_[idx] = mass_[idx];
      decoy_tf->mass_calculated_[idx] = true;
      decoy_tf->mass_[idx] = mass_[idx];
    }
  }
  decoys.push_back(decoy_ff);
  decoys.push_back(decoy_tf);
  decoys.push_back(decoy_ft);
}

/**
 * fills the ion series with the predicted ions for the cross linked candidate
 */ 
void XLinkPeptide::predictIons(
  IonSeries* ion_series,  ///< IonSeries object to fill
  int charge, ///< charge state of candidate
  bool first ///< is this the first peptide?
  ) {
  carp(CARP_DEBUG, "predictIons:start %i %i", charge, first?1:0);

  MASS_TYPE_T fragment_mass_type = GlobalParams::getFragmentMass();

  if (first) {
    linked_peptides_[0].predictIons(
      ion_series, charge, getLinkIdx(0), 
      linker_mass_ + linked_peptides_[0].getMass(fragment_mass_type)); 
  } else {
    linked_peptides_[1].predictIons(
      ion_series, charge, getLinkIdx(1),
      linker_mass_ + linked_peptides_[1].getMass(fragment_mass_type));
  }
} 

/**
 * predicts the ions for xlinked peptide
 */
void XLinkPeptide::predictIons(
  IonSeries* ion_series, ///< IonSeries to fill
  int charge ///< charge state of the peptide
  ) {

  MASS_TYPE_T fragment_mass_type = GlobalParams::getFragmentMass();
  FLOAT_T delta_mass0 =  linked_peptides_[0].getMass(fragment_mass_type) + linker_mass_;
  FLOAT_T delta_mass1 =  linked_peptides_[1].getMass(fragment_mass_type) + linker_mass_;

  //predict the ion_series of the first peptide.
  linked_peptides_[0].predictIons(ion_series, charge, getLinkIdx(0), delta_mass1, true);

  //predict the ion_series of the second peptide.
  linked_peptides_[1].predictIons(ion_series, charge, getLinkIdx(1), delta_mass0, false);
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
      const char* seq2 = linked_peptides_[1].getSequence();
      ans = subseq + string(",") + string(seq2);
    } else {
      const char* seq1 = linked_peptides_[0].getSequence();
      ans = string(seq1) + string(",") + subseq;
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
XLinkablePeptide& XLinkPeptide::getXLinkablePeptide(
  int peptide_idx ///< 0 or 1
  ) {
  return linked_peptides_[peptide_idx];
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

  return missed1 + missed2;

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
  }
  
  vector<string> prot1 = peptide -> getProteinIds();
  string prefix = "";
  if (linked_peptides_[0].isDecoy()) {
    prefix = Params::GetString("decoy-prefix");
  }
  
  oss << prefix << prot1[0];
  for (size_t idx = 1 ; idx < prot1.size() ; idx++) {
    oss << "," << prefix << prot1[idx];
  }
  oss << ";";

  peptide = this -> getPeptide(1);

  if (peptide == NULL) {
    carp(CARP_FATAL, "XLinkPeptide : Null second peptide!");
  }
  vector<string> prot2 = peptide -> getProteinIds();
  prefix = "";
  if (linked_peptides_[1].isDecoy()) {
    prefix = Params::GetString("decoy-prefix");
  }
  oss << prefix << prot2[0];
  for (size_t idx = 1;idx < prot2.size() ; idx++) {
    oss << "," << prefix << prot2[idx];
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
