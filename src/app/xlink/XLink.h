/**
 * \file XLink.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Utility functions for search-for-xlinks
 *****************************************************************************/
#ifndef XLINK_H_
#define XLINK_H_
#include "model/objects.h"
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
 * \returns whether the proposed xlink peptides would survive the inter/intra filters set by the user
 */
bool testInterIntraKeep(Crux::Peptide* pep1, Crux::Peptide* pep2);

bool testInterIntraKeep(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2,
  bool include_intra,
  bool include_inter,
  bool include_inter_intra
  );

/**
 * \ returns the crosslink candidate type between two peptides
 */
XLINKMATCH_TYPE_T getCrossLinkCandidateType(
  Crux::Peptide* pep1, 
  Crux::Peptide* pep2
  );

/**
 * \returns whether two proposed peptides would contain an inter-protein crosslink
 */
bool isCrossLinkInter(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2
  );

/**
 * \returns whether two proposed peptides would contain an intra-protein crosslink
 */
bool isCrossLinkIntra(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2
  );

/**
 * \returns whether two propsed peptides would contain both an intra and inter protein cross link
 */ 
bool isCrossLinkInterIntra(
  Crux::Peptide* pep1,
  Crux::Peptide* pep2
);

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

} // namespace XLink

#endif

/*
 * Local Variables:
 * mode: c 
 * c-basic-offset: 2
 * End:
 */

