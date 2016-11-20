/**
 * \file XLink.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Utility functions for search-for-xlinks
 *****************************************************************************/
#include "XLink.h"
#include "util/Params.h"

#include <sstream>
#include <iostream>
using namespace std;

namespace XLink {

set<Crux::Peptide*> allocated_peptides_; ///< tracker for allocated peptides

bool testInterIntraKeep(
  Crux::Peptide *pep1,
  Crux::Peptide *pep2
  ) {
  return(testInterIntraKeep(pep1, pep2, Params::GetBool("xlink-include-intra"), 
    Params::GetBool("xlink-include-inter"),
    Params::GetBool("xlink-include-inter-intra")));
}

bool testInterIntraKeep(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2,
  bool include_intra,
  bool include_inter,
  bool include_inter_intra
  ) {
  //  carp(CARP_INFO, "testInterIntraKeep:start");
  if (include_intra && 
      include_inter && 
      include_inter_intra) {
    //no need to test
    return true;
  }

  XLINKMATCH_TYPE_T cross_link_type = getCrossLinkCandidateType(pep1, pep2);

  if (cross_link_type == XLINK_INTER_INTRA_CANDIDATE  && include_inter_intra) {
    return true;
  } else if (cross_link_type == XLINK_INTER_CANDIDATE && include_inter) {
    return true;
  } else if (cross_link_type == XLINK_INTRA_CANDIDATE && include_intra) {
    return true;
  }
  return false;
}

XLINKMATCH_TYPE_T getCrossLinkCandidateType(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2
  ) {

  bool is_intra = false;
  bool is_inter = false;
  size_t id1 = -1;
  size_t id2 = -1;

  for (PeptideSrcIterator src_iterator1 = pep1->getPeptideSrcBegin();
    src_iterator1 != pep1->getPeptideSrcEnd();
    ++src_iterator1) {
    PeptideSrc* src1 = *src_iterator1;
    size_t id1t = src1->getParentProtein()->getProteinIdx();
    if (id1 == -1) {
      id1 = id1t;
    } else if (id1 != id1t) {
      return XLINK_INTER_INTRA_CANDIDATE;
    }
    for (PeptideSrcIterator src_iterator2 = pep2->getPeptideSrcBegin();
      src_iterator2 != pep2->getPeptideSrcEnd();
      ++src_iterator2) {
      PeptideSrc* src2 = *src_iterator2;
      size_t id2t = src2->getParentProtein()->getProteinIdx();
      if (id2 == -1) {
        id2 = id2t;
      } else if (id2 != id2t) {
        return XLINK_INTER_INTRA_CANDIDATE;
      }
      if (id1 == id2) {
        is_intra = true;
        if (is_inter) {
          return(XLINK_INTER_INTRA_CANDIDATE);
        }
      } else {
        is_inter = true;
        if (is_intra) {
          return(XLINK_INTER_INTRA_CANDIDATE);
        }
      }
    }
  }

  if (is_intra && !is_inter) {
    return (XLINK_INTRA_CANDIDATE);
  } else if (is_inter & !is_intra) {
    return(XLINK_INTER_CANDIDATE);
  } else {
    carp(CARP_FATAL, "Internal error at getCrossLinkCandidateType");
  }
  return(INVALID_CANDIDATE);
}

/**
 * \returns whether two proposed peptides would contain an intra-protein crosslink
 */
bool isCrossLinkIntra(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2
  ) {

  return(getCrossLinkCandidateType(pep1, pep2) == XLINK_INTRA_CANDIDATE);
}

/**
 * \returns whether two proposed peptides would contain an inter-protein crosslink
 */
bool isCrossLinkInter(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2
  ) {
  
  return(getCrossLinkCandidateType(pep1, pep2) == XLINK_INTER_CANDIDATE);
}

bool isCrossLinkInterIntra(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2
  ) {

  return(getCrossLinkCandidateType(pep1, pep2) == XLINK_INTER_INTRA_CANDIDATE);
}


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
  carp(CARP_DEBUG, "deleting %d peptides", allocated_peptides_.size());
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

  string prefix = "";
  if (peptide->isDecoy()) {
    prefix = Params::GetString("decoy-prefix");
  }
  
  for (PeptideSrcIterator peptide_src_iterator =
    peptide->getPeptideSrcBegin();
    peptide_src_iterator != peptide->getPeptideSrcEnd();
    ++peptide_src_iterator) {
    
    PeptideSrc* peptide_src = *peptide_src_iterator;
    Crux::Protein* protein = peptide_src->getParentProtein();
    string& protein_id = protein->getIdPointer();
    int peptide_loc = peptide_src->getStartIdx();
    std::ostringstream protein_loc_stream;
    protein_loc_stream << prefix << protein_id << "(" << peptide_loc << ")";
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

} // namespace XLink


/*
 * Local Variables:
 * mode: c 
 * c-basic-offset: 2
 * End:
 */

