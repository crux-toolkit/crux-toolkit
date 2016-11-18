/**
 * \file XLinkablePeptideIterator.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Iterator for xlinkable peptides
 *****************************************************************************/

#ifndef XLINKABLEPEPTIDEITERATORTOPN_H
#define XLINKABLEPEPTIDEITERATORTOPN_H

#include "XLinkablePeptideIterator.h"
#include "XLinkablePeptide.h"
#include "XLinkBondMap.h"
#include "XLinkScorer.h"
#include <queue>
#include <vector>


class XLinkablePeptideIteratorTopN: public XLinkablePeptideIterator {

 protected:

  //std::priority_queue<XLinkablePeptide, std::vector<XLinkablePeptide>, CompareXCorr> scored_xlp_;  
  std::vector<XLinkablePeptide*> scored_xlp_; ///< sorted by highest XCorr score.
  int current_count_;
  int top_n_; ///<set by kojak-top-n
  bool has_next_; ///< is there a next candidate
  bool is_decoy_; ///< are we getting decoys

  /**
   * queues the next linkable peptide
   */
  void queueNextPeptide(); 
  
  void scorePeptides(
    XLinkScorer& scorer,
    FLOAT_T precursor_mass,
    std::vector<XLinkablePeptide>::iterator& biter,
    std::vector<XLinkablePeptide>::iterator& eiter
  );

 public:

  /**
   * constructor that sets up the iterator
   */
  XLinkablePeptideIteratorTopN(
    Crux::Spectrum* spectrum,
    FLOAT_T precursor_mass,
    FLOAT_T min_mass, ///< min mass of candidates
    FLOAT_T max_mass, ///< max mass of candidates
    int precursor_charge, ///< Charge of precursor
    bool is_decoy ///< generate decoy candidates
    );

  /**
   * Destructor
   */
  virtual ~XLinkablePeptideIteratorTopN();

  /**
   * \returns whether there is another linkable peptide
   */
  bool hasNext();

  /**
   *\returns the next linkable peptide
   */
  XLinkablePeptide& next();
  
  XLinkablePeptide* nextPtr();

};


#endif

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
