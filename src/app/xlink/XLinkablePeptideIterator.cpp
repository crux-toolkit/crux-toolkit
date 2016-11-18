/**
 * \file XLinkablePeptideIterator.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Iterator for xlinkable peptides
 *****************************************************************************/

#include "XLinkablePeptideIterator.h"
#include "XLink.h"
#include "LinearPeptide.h"
#include "util/GlobalParams.h"
#include <iostream>
#include "XLinkDatabase.h"

using namespace std;

/**
 * constructor that sets up the iterator
 */
XLinkablePeptideIterator::XLinkablePeptideIterator(
    double min_mass, ///< min mass of candidates
    double max_mass, ///< max mass of candidates
    Database* database, ///< protein database
    PEPTIDE_MOD_T** peptide_mods, ///<current peptide mod
    int num_peptide_mods,
    bool is_decoy, ///< generate decoy candidates
    XLinkBondMap& bondmap ///< map of valid links
    ) {

  is_decoy_ = is_decoy;

  iter_ = XLinkDatabase::getXLinkableBegin(min_mass);
  eiter_ = XLinkDatabase::getXLinkableEnd();

  min_mass_ = min_mass;
  max_mass_ = max_mass;
  has_next_ = iter_ != eiter_ && iter_->getMass(GlobalParams::getIsotopicMass()) <= max_mass_;
}

/**
 * Destructor
 */
XLinkablePeptideIterator::~XLinkablePeptideIterator() {
  ;
}

/**
 * queues the next linkable peptide
 */
void XLinkablePeptideIterator::queueNextPeptide() {
  
  iter_++;
  has_next_ = iter_ != eiter_ && iter_->getMass(GlobalParams::getIsotopicMass()) <= max_mass_;
}

/**
 *\returns whether there is another linkable peptide
 */
bool XLinkablePeptideIterator::hasNext() {

  return has_next_;
}

/**
 *\returns the next peptide
 */
XLinkablePeptide& XLinkablePeptideIterator::next() {

  if (!has_next_) {
    carp(CARP_WARNING, "next called on empty iterator!");
  }
  
  XLinkablePeptide& ans = *iter_;
  queueNextPeptide();
  return ans;
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
