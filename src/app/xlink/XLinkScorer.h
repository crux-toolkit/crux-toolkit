/**
 * \file XLinkScorer.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for scoring xlink products.  Tries to optimize the scoring time
 *****************************************************************************/

#ifndef XLINKSCORER_H_
#define XLINKSCORER_H_
#include "objects.h"
#include "XLinkMatch.h"

class XLinkScorer {
 protected:
  Crux::Spectrum* spectrum_; ///< spectrum object
  Scorer* scorer_xcorr_; ///< xcorr scorer object
  Scorer* scorer_sp_; ///< sp scorer object
  IonConstraint* ion_constraint_xcorr_; ///< xcorr constraint
  IonConstraint* ion_constraint_sp_; ///< sp constraint
  XLinkMatch* candidate_; ///< current candidate
  int charge_; ///< current charege
  IonSeries* ion_series_xcorr_; ///< current ion series xcorr
  IonSeries* ion_series_sp_; ///< current ion series sp
  bool compute_sp_; ///< calculate sp score
 
  /**
   * initializes the object with the spectrum
   * charge and appropriate objects
   */
  void init(
    Crux::Spectrum* spectrum, ///< spectrum to score
    int charge, ///< charge state
    bool compute_sp ///< are we scoring sp?
    );

 public:
  /**
   * default constructor
   */
  XLinkScorer();

  /**
   * Constructor for spectrum and charge
   */
  XLinkScorer(
    Crux::Spectrum* spectrum, ///< spectrum to score
    int charge ///< charge state
    );
  
  /**
   * Constructor that allows compute sp to be set
   */
  XLinkScorer(
    Crux::Spectrum* spectrum, ///spectrum to score
    int charge, ///< charge state
    bool compute_sp ///< scoring sp?
    );

  /**
   * destructor, free up the objects
   */
  virtual ~XLinkScorer();

  /**
   * \returns the xcorr score and sets the xcorr and sp scores to the match
   */
  FLOAT_T scoreCandidate(XLinkMatch* candidate);


};


#endif
/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
