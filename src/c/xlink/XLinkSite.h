/**
 * \file XLinkSite.h
 * $Revision: 1.00 $
 * \brief Object for representing a cross-linkable site on a peptide.
 *********************************************************************/
#ifndef XLINKSITE_H
#define XLINKSITE_H

#include <string>

#include "Peptide.h"

/**
 * \class XLinkSite
 * \brief object respresneting a cross-linkable site on a peptide.
 */
class XLinkSite {

 protected:
  XLINK_SITE_T type_; ///<The type of the cross-link
  char aa_;           ///<Amino acid type (for XLINKSITE_AA types).
 
 public:

  /**
   * Default constructor
   */
  XLinkSite();

  /**
   * Default destructor
   */
  virtual ~XLinkSite();

  /**
   * Constructor that determines the type from the supplied site string 
   */
  XLinkSite(
    std::string& site_string ///<string desscribing the cross-linkable site
    );

  /**
   * \returns whether the peptide contains this site at the supplied sequence index.
   */
  bool hasSite(
    Crux::Peptide* peptide, ///<peptide object pointer 
    int idx             ///<sequence index
    ) const;


  bool hasSite(
    std::string& protein_sequence,
    int idx
  ) const;
    
  /**
   * \returns whether this xlinksite is equal to the passed in xlinksite
   */
  bool operator == (
    XLinkSite& xlink_site_obj ///<XLinkSite to compare to
    );

  /**
   * \returns whether this xlinksite is less than the passed in xlinksite
   */
  bool operator < (
    const XLinkSite& xlink_site_obj ///<XLinkSite to compare to
    ) const;

};

#endif
/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */

