/**
 * \file XLink.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Utility functions for search-for-xlinks
 *****************************************************************************/
#ifndef XLINK_H_
#define XLINK_H_
#include "objects.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>
#include <string>


class XLinkMatch;
class XLinkMatchCollection;


void get_min_max_mass(
  FLOAT_T precursor_mz, 
  int charge, 
  bool use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass);



namespace XLink {

/**
 * \returns protein ids with start locations marked with (X)
 */
std::string get_protein_ids_locations(
  Crux::Peptide* peptide ///< peptide to get locations from
  );

/**
 * add a peptide to the list of allocated peptides
 */
void addAllocatedPeptide(
  Crux::Peptide* peptide ///< peptide to add
  );

/**
 * delete all peptides that are allocated
 */
void deleteAllocatedPeptides();

};
  
#endif

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
