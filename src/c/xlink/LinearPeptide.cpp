
#include "LinearPeptide.h"
#include "XLink.h"
#include "ModifiedPeptidesIterator.h"
#include "IonSeries.h"

#include <iostream>

using namespace std;

/**
 * Default constructor
 */
LinearPeptide::LinearPeptide() {
  peptide_ = NULL;
  sequence_ = NULL;
}

/**
 * Constructor with sequence
 */
LinearPeptide::LinearPeptide(
  char* sequence ///< sequence string
  ) {
  peptide_ = NULL;
  sequence_ = sequence;
}

/**
 * Constructor from a Crux Peptide
 */
LinearPeptide::LinearPeptide(
  Crux::Peptide* peptide ///< peptide object
  ) {
  peptide_ = peptide;
  sequence_ = NULL;
}

/**
 *Add candidates to the XLinkMatchCollection that are linear
 */
void LinearPeptide::addCandidates(
  FLOAT_T min_mass,  ///< min mass
  FLOAT_T max_mass,  ///< max mass
  Index* index,  ///< protein index
  Database* database, ///< protein database
  PEPTIDE_MOD_T** peptide_mods, ///< modifications peptide can take
  int num_peptide_mods, ///< Number of possible peptide mods
  XLinkMatchCollection& candidates ///< Vector of candidate -inout
  ) {

  int max_missed_cleavages = get_int_parameter("missed-cleavages");

  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    double delta_mass = peptide_mod_get_mass_change(peptide_mod);
    //
    ModifiedPeptidesIterator* peptide_iterator =
      new ModifiedPeptidesIterator(min_mass - delta_mass, max_mass - delta_mass, peptide_mod, 
        false, index, database);

    //add the targets
    while (peptide_iterator->hasNext()) {
      Crux::Peptide* peptide = peptide_iterator->next();
      XLinkMatch* new_candidate = new LinearPeptide(peptide);
      if (new_candidate->getNumMissedCleavages() <= max_missed_cleavages) {
        //cerr <<"Adding Linear Peptide:"<<new_candidate -> getSequenceString()<<" "<<new_candidate->getMass()<<endl;
        candidates.add(new_candidate);
        XLink::addAllocatedPeptide(peptide);
      } else {
        delete new_candidate;
        delete peptide;
      }
    }
    delete peptide_iterator;


    //add the decoys
    peptide_iterator = new
      ModifiedPeptidesIterator(min_mass - delta_mass, max_mass - delta_mass, peptide_mod,
        true, index, database);

    while (peptide_iterator->hasNext()) {
      Crux::Peptide* peptide = peptide_iterator->next();
      XLinkMatch* new_candidate = new LinearPeptide(peptide);
      if (new_candidate->getNumMissedCleavages() <= max_missed_cleavages) {
        candidates.add(new_candidate);
        XLink::addAllocatedPeptide(peptide);
      } else {
        delete new_candidate;
        delete peptide;
      }
    }

    delete peptide_iterator;


  }
}

/**
 * returns the candidate type, either a deadlink or a linear candidate
 */
XLINKMATCH_TYPE_T LinearPeptide::getCandidateType() {
  if (isModified()) {
    return DEADLINK_CANDIDATE;
  } else {
    return LINEAR_CANDIDATE;
  }
}

/**
 * \returns the sequence of the peptide in string format
 */
string LinearPeptide::getSequenceString() {
  ostringstream oss;
  if (peptide_ == NULL) {
    oss << sequence_;

  } else {
    char* seq = peptide_->getModifiedSequenceWithMasses(MOD_MASSES_SEPARATE);
    oss << seq;
    free(seq);
  }

  oss << " ()";
  return oss.str();

}

/**
 * \returns the mass of the peptide
 */
FLOAT_T LinearPeptide::calcMass(
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ) {
  
  if (peptide_ == NULL) {
    return Crux::Peptide::calcSequenceMass(sequence_, mass_type);
  } else {
    return peptide_->calcModifiedMass(mass_type);
  }
}

/**
 *\returns a shuffled version of the peptide
 */
XLinkMatch* LinearPeptide::shuffle() {

  Crux::Peptide* decoy_peptide = new Crux::Peptide(peptide_);

  decoy_peptide->transformToDecoy();

  XLink::addAllocatedPeptide(decoy_peptide);

  LinearPeptide* decoy = new LinearPeptide(decoy_peptide);

  return (XLinkMatch*)decoy;
}

/**
 * predicts the ions for this peptide
 */
void LinearPeptide::predictIons(
  IonSeries* ion_series, ///< ion series to place the ions
  int charge ///< charge state of the peptide
  ) {

  char* seq = NULL;
  MODIFIED_AA_T* mod_seq = NULL;
  if (peptide_ == NULL) {
    seq = my_copy_string(sequence_);
    convert_to_mod_aa_seq(seq, &mod_seq); 
  } else {
    seq = peptide_->getSequence();
    mod_seq = peptide_->getModifiedAASequence();
  }
  ion_series->setCharge(charge);
  ion_series->update(seq, mod_seq);
  ion_series->predictIons();
  free(seq);
  free(mod_seq);

}

/**
 *\returns the ion sequence as a string
 */
string LinearPeptide::getIonSequence(
  Ion* ion ///< ion object
  ) {

  string seq_str = string(sequence_);

  int cleavage_idx = ion->getCleavageIdx();
  if (ion->isForwardType() == B_ION) {
    return seq_str.substr(0,cleavage_idx);
  } else {
    return seq_str.substr(seq_str.length()-cleavage_idx,seq_str.length());
  }
}

/**
 * \returns the peptide for this match
 */
Crux::Peptide* LinearPeptide::getPeptide(int peptide_idx) {
  if (peptide_idx == 0) {
    return peptide_;
  } else {
    return NULL;
  }
}

/**
 *\returns the number of missed cleavages
 */
int LinearPeptide::getNumMissedCleavages() {
  set<int> skip;
  return peptide_->getMissedCleavageSites(skip);
}

/**
 * \returns whether this peptide is modified by a variable mod
 */
bool LinearPeptide::isModified() {

  return peptide_->isModified();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
