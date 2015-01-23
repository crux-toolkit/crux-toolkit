/**
 * \file LinkedIonSeries.cpp 
 * AUTHOR: Sean McIlwain and Paul Draghicescu
 * CREATE DATE: 18 January 2012
 * \brief Object for keeping track of a (non-)crosslinked peptide ions.
 *****************************************************************************/
#include "LinkedIonSeries.h"
#include "xhhc_scorer.h"

#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

#include <iostream>

/**
 * GLOBALS - TODO GET RID OF THESE!
 */
FLOAT_T* mass_matrix = NULL;
int mass_matrix_size = -1;


/**
 * \returns a blank LinkedIonSeries object
 */
LinkedIonSeries::LinkedIonSeries() {
  charge_ = 0;
  fragment_mass_type_ = get_mass_type_parameter("fragment-mass");
}

/**
 * \returns an ion series object assigning the charge.
 */
LinkedIonSeries::LinkedIonSeries(
  int charge ///< The maximum charge of the ion series.
  ) {

  charge_ = charge;
  fragment_mass_type_ = get_mass_type_parameter("fragment-mass");
}

/**
 * Default desctructor
 */
LinkedIonSeries::~LinkedIonSeries() {

}

/**
 * \returns charge for this linked ion series
 */
int LinkedIonSeries::getCharge() {
  return charge_;
}

/**
 * \returns all the ions in this series
 */
vector<LinkedPeptide>& LinkedIonSeries::getIons() {
  return all_ions_;
}

/**
 * \returns the number of ions within this series
 */
int LinkedIonSeries::getSize() {
  return all_ions_.size();
}

/**
 * removes all ions from this series
 */
void LinkedIonSeries::clear() {
  all_ions_.clear();
}

/**
 * prints out tab delimited information about the ion series
 */
void LinkedIonSeries::print() {

  LinkedPeptide::sortByMZ(all_ions_, fragment_mass_type_);

  cout << "m/z\ttype\tcharge\tion" << endl;
  string ion_type;
  for (vector<LinkedPeptide>::iterator ion = all_ions_.begin(); 
    ion != all_ions_.end(); 
    ++ion) {

    if (ion->getIonType() == B_ION) 
      ion_type = "B_ION";
    else 
      ion_type = "Y_ION";
    
    cout << ion->getMZ(fragment_mass_type_) << "\t" << 
            ion_type << "\t" << 
            ion->getCharge() << "\t" << 
            *ion << endl;
  }
}

/** 
 * cleaves linked_peptide at all positions, adding b and y ions
 */
void LinkedIonSeries::addLinkedIons(
  LinkedPeptide& linked_peptide, ///< The linked peptide
  SPLITTYPE_T split_type ///< Which peptide to split (A,B, or BOTH)
  ) {

  linked_peptide.setCharge(charge_);

  fragments_.clear();
  // split the precursor at every cleavage site

  switch(split_type) {
    case SPLITTYPE_A:
      linked_peptide.splitA(fragments_);
      break;
    case SPLITTYPE_B:
      linked_peptide.splitB(fragments_);
      break;
    case SPLITTYPE_BOTH:
      linked_peptide.split(fragments_);
      break;
    case SPLITTYPE_INVALID:
    case NUMBER_SPLITTYPES:
      carp(CARP_FATAL, "SPLITTYPE invalid:%d",(int)split_type);
  }
  
  for (vector<pair<LinkedPeptide, LinkedPeptide> >::iterator ion_pair = 
    fragments_.begin(); 
    ion_pair != fragments_.end(); 
    ++ion_pair) {

    // if b-ion and not a neutral loss
    if (ion_pair->first.getCharge() != 0) {
      ion_pair->first.setIonType(B_ION); 
      ion_pair->first.getMZ(fragment_mass_type_);
      all_ions_.push_back(ion_pair->first);
    }
    // if y-ion and not a neutral loss
    if (ion_pair->second.getCharge() != 0) {
      ion_pair->second.setIonType(Y_ION); 
      ion_pair->second.getMZ(fragment_mass_type_);
      all_ions_.push_back(ion_pair->second);
    }
  }
  
}

/**
 * gets the Total number of by ions in this series.
 */ 
int LinkedIonSeries::getTotalBYIons() {
   
  int ans = 0;
  vector<LinkedPeptide>::iterator ion_iter;

  for (ion_iter = all_ions_.begin(); 
    ion_iter != all_ions_.end(); 
     ++ion_iter) {
	
    if (ion_iter -> getIonType() == B_ION || ion_iter -> getIonType() == Y_ION) {
      ans++;
    }
  }
  return ans;
}

/**
 * adds an ion to the observed list, used in conjunction with 
 * getObservableIons
 */
void LinkedIonSeries::addIonBin(
  map<int, bool>& observed, ///< the observed bin vector -in/out
  int& ions, ///< the number of ions -in/out
  int& ions_bin, ///< the number of binned ions -in/out
  FLOAT_T mz, ///< m/z of the ion
  FLOAT_T bin_width, ///< bin-width
  FLOAT_T min_mz, ///< max m/z for the range
  FLOAT_T max_mz, ///< min m/z for the range
  bool add_flanks ///< add flanks?
  ) {

  if ((mz >= min_mz) && (mz <= max_mz)) {
    ions++;
    int bin_idx = (int)(mz / bin_width + 0.5);
    if (observed.find(bin_idx) == observed.end()) {
      observed[bin_idx] = true;
      ions_bin++;
    }
    if (add_flanks) {
      
      bin_idx = bin_idx - 1;
      FLOAT_T flank_mz = bin_idx * bin_width;
      if (flank_mz >= min_mz) {
        ions++;
        if (observed.find(bin_idx) == observed.end()) {
          observed[bin_idx] = true;
          ions_bin++;
        }
      }
      bin_idx = bin_idx + 2;
      flank_mz = bin_idx * bin_width;
      if (flank_mz <= max_mz) {
        ions++;
        if (observed.find(bin_idx) == observed.end()) {
          observed[bin_idx] = true;
          ions_bin++;
        }
      }
    }
  }
}

/**
 * \returns the total number of ions observable within a binned
 * range
 */
int LinkedIonSeries::getObservableIons(
  FLOAT_T min_mz, ///< the minimum m/z -in
  FLOAT_T max_mz, ///< the maximum m/z -in
  FLOAT_T bin_width, ///< the width of the bins -in
  int& ions_observable, ///<number of ions observable -out
  int& ions_observable_bin ///<number of bins observable -out
  ) {
  
  ions_observable = 0;
  ions_observable_bin = 0;

  map<int, bool> observed;

  for (vector<LinkedPeptide>::iterator ion_iter = all_ions_.begin();
    ion_iter != all_ions_.end();
    ++ion_iter) {

    double mz = ion_iter -> getMZ(MONO);

    if (ion_iter -> getIonType() == B_ION || 
      ion_iter -> getIonType() == Y_ION) {
      addIonBin(observed, ions_observable, ions_observable_bin, mz, bin_width, min_mz, max_mz, true);
      double h2O_mz = ion_iter -> getMZ(MONO) - (MASS_H2O_MONO/ ion_iter -> getCharge());
      addIonBin(observed, ions_observable, ions_observable_bin, h2O_mz, bin_width, min_mz, max_mz, false);
      double nh3_mz = ion_iter -> getMZ(MONO) - (MASS_NH3_MONO/ ion_iter -> getCharge());
      addIonBin(observed, ions_observable, ions_observable_bin, nh3_mz, bin_width, min_mz, max_mz, false);
      double co_mz = ion_iter -> getMZ(MONO) - (MASS_CO_MONO / ion_iter -> getCharge());
      addIonBin(observed, ions_observable, ions_observable_bin, co_mz, bin_width, min_mz, max_mz, false);
    }
  }
  return 0;
}

/**
 * \returns the total number of b-y ions observable within a binned
 * range
 */  
int LinkedIonSeries::getObservableBYIons(
  FLOAT_T min_mz, ///< the minimum m/z -in 
  FLOAT_T max_mz, ///< the maximum m/z -in
  FLOAT_T bin_width, ///< the width of the bins -in
  int &by_observable, ///<number of ions observable -out
  int &by_observable_bin ///<number of bins observable -out
  ) {

  by_observable = 0;
  by_observable_bin = 0;

  map<int, bool> observed;

  vector<LinkedPeptide>::iterator ion_iter;

  for (ion_iter = all_ions_.begin(); 
       ion_iter != all_ions_.end(); 
       ++ion_iter) {
    if (ion_iter -> getMZ(MONO) >= min_mz && ion_iter -> getMZ(MONO) <= max_mz) {
      if (ion_iter ->getIonType() == B_ION || ion_iter -> getIonType() == Y_ION) {
        by_observable++;
        int bin_idx = 
          (int)(ion_iter->getMZ(MONO) / bin_width + 0.5);
        if (observed.find(bin_idx) == observed.end()) {
          observed[bin_idx] = true;
          by_observable_bin++;
        }
      }
    }
  }
  return by_observable_bin;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
