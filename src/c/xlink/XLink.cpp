/**
 * \file XLink.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Utility functions for search-for-xlinks
 *****************************************************************************/
#include "XLink.h"


#include <sstream>

using namespace std;

namespace XLink {

set<Crux::Peptide*> allocated_peptides_; ///< tracker for allocated peptides

/**
 * add a peptide to the list of allocated peptides
 */
void addAllocatedPeptide(
  Crux::Peptide* peptide ///< peptide to add
  ) {

  allocated_peptides_.insert(peptide);
}

/**
 * delete all peptides that are allocated
 */
void deleteAllocatedPeptides() {
  for (set<Crux::Peptide*>::iterator iter =
    allocated_peptides_.begin();
    iter != allocated_peptides_.end();
    ++iter) {
  
  delete *iter;

  }
  allocated_peptides_.clear();
}


void get_protein_ids_locations(
  Crux::Peptide *peptide,  ///< peptide
  set<string>& protein_ids_locations ///< location strings
  ) {

  std::ostringstream protein_field_stream;

  for (PeptideSrcIterator peptide_src_iterator =
    peptide->getPeptideSrcBegin();
    peptide_src_iterator != peptide->getPeptideSrcEnd();
    ++peptide_src_iterator) {
    
    PeptideSrc* peptide_src = *peptide_src_iterator;
    Crux::Protein* protein = peptide_src->getParentProtein();
    char* protein_id = protein->getId();
    int peptide_loc = peptide_src->getStartIdx();
    std::ostringstream protein_loc_stream;
    protein_loc_stream << protein_id << "(" << peptide_loc << ")";
    free(protein_id);
    protein_ids_locations.insert(protein_loc_stream.str());
  }

}

string get_protein_ids_locations(
  Crux::Peptide* peptide ///< peptide to get protein strings from
  ) {
  set<string> protein_ids_locations;

  get_protein_ids_locations(peptide, protein_ids_locations);
  set<string>::iterator result_iter = protein_ids_locations.begin();

  string protein_field_string = *result_iter;

  while(++result_iter != protein_ids_locations.end()) {
    protein_field_string += "," + *result_iter;
  }

  return protein_field_string;

}





};

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */

