/**
 * \file XLinkablePeptideIterator.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Iterator for xlinkable peptides
 *****************************************************************************/

#include "XLinkablePeptideIterator.h"
#include "XLink.h"
#include <iostream>

using namespace std;

/**
 * constructor that sets up the iterator
 */
XLinkablePeptideIterator::XLinkablePeptideIterator(
    double min_mass, ///< min mass of candidates
    double max_mass, ///< max mass of candidates
    Index* index, ///< protein index
    Database* database, ///< protein database
    PEPTIDE_MOD_T* peptide_mod, ///<current peptide mod
    bool is_decoy, ///< generate decoy candidates
    XLinkBondMap& bondmap ///< map of valid links
    ) {

  is_decoy_ = is_decoy;

  bondmap_ = bondmap;

  peptide_iterator_ =     
    new ModifiedPeptidesIterator(
      min_mass, 
      max_mass,
      peptide_mod, 
      is_decoy,
      index, 
      database);
  queueNextPeptide();

  

}

/**
 * Destructor
 */
XLinkablePeptideIterator::~XLinkablePeptideIterator() {
  delete peptide_iterator_;
}

/**
 * queues the next linkable peptide
 */
void XLinkablePeptideIterator::queueNextPeptide() {

  has_next_ = false;
  int max_mod_xlink = get_int_parameter("max-xlink-mods");
  while (peptide_iterator_->hasNext() && !has_next_) {

    Crux::Peptide* peptide = peptide_iterator_->next();
    //cerr << "peptide is:"<<peptide->getSequence()<<endl;  
    if (peptide->countModifiedAAs() <= max_mod_xlink) {
      XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites_);
      if (link_sites_.size() > 0) {
        has_next_ = true;
        current_ = XLinkablePeptide(peptide, link_sites_);
        current_.setDecoy(is_decoy_);
        XLink::addAllocatedPeptide(peptide);
      }
    } 
    if (!has_next_) {
      delete peptide;
    }
  }
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
XLinkablePeptide XLinkablePeptideIterator::next() {

  if (!has_next_) {
    carp(CARP_WARNING, "next called on empty iterator!");
  }

  XLinkablePeptide ans = current_;
  queueNextPeptide();
  return ans;
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
