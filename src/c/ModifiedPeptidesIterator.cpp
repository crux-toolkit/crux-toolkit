/**
 * \file GenerateModifiedPeptidesIterator.h
 * \brief An object to return candidate peptides from a database or index.
 */
#include "ModifiedPeptidesIterator.h"
#include "SpectrumZState.h"

using namespace std;
using namespace Crux;

/**
 * Constructor that sets all fields with the given values.  Will
 * return peptides in a window around the given m/z or mass with the
 * given modifications applied.
 */
ModifiedPeptidesIterator::ModifiedPeptidesIterator(
  double mz,               ///< Spectrum precrusor mz
  SpectrumZState& zstate,  ///< Target mz of peptides
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  bool is_decoy,  ///< generate decoy peptides
  Index* index,      ///< Index from which to draw peptides OR
  Database* dbase    ///< Database from which to draw peptides
)
{
  peptide_source_ = new GeneratePeptidesIterator(getMinMaxMass(mz, zstate, 
                                                               pmod),
                                                 is_decoy, dbase, index);
  peptide_modification_ = pmod;
  temp_peptide_list_ = new_empty_list();
  max_aas_modified_ = get_int_parameter("max-aas-modified");

  initialize();
  carp(CARP_DETAILED_DEBUG, 
       "After initialize, has next is %d and next peptide is null? %d",
       hasNext(), (next_peptide_ == NULL));
}

/**
 * Constructor for returning all peptides in the index or database
 * that fall within the constraints defined in parameter.cpp.
 */
ModifiedPeptidesIterator::ModifiedPeptidesIterator(
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  Index* index,      ///< Index from which to draw peptides OR
  Database* dbase    ///< Database from which to draw peptides
){
  peptide_source_ = new GeneratePeptidesIterator(getMinMaxMass(), 
                                                 false, // not decoy
                                                 dbase, index);

  peptide_modification_ = pmod;
  temp_peptide_list_ = new_empty_list();
  max_aas_modified_ = get_int_parameter("max-aas-modified");
  initialize();
}

/**
 * Constructor for returnign all peptides in the index or database
 * that fall within the mass range
 */
ModifiedPeptidesIterator::ModifiedPeptidesIterator(
  double min_mass,    ///< min-mass of peptides
  double max_mass,    ///< max-mass of peptides
  PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
  bool is_decoy, ///< generate decoy peptides
  Index* index,     ///< Index from which to draw peptides OR
  Database* dbase   ///< Database from which to draw peptides
) {

  peptide_source_ = new GeneratePeptidesIterator(
    pair<FLOAT_T,FLOAT_T>(min_mass, max_mass),
    is_decoy,
    dbase, index);

  peptide_modification_ = pmod;
  temp_peptide_list_ = new_empty_list();
  max_aas_modified_ = get_int_parameter("max-aas-modified");
  initialize();
}


/**
 * Destructor that frees any peptides not yet returned via next().
 */
ModifiedPeptidesIterator::~ModifiedPeptidesIterator(){
  delete_linked_list(temp_peptide_list_);
  delete next_peptide_;
  delete peptide_source_;
}

/**
 * Use parameter.cpp values to determine the range of masses to
 * consider.
 * \returns A pair with the minimum and maximum mass for generating
 * peptides.  
 */
pair<FLOAT_T,FLOAT_T>ModifiedPeptidesIterator::getMinMaxMass()
{
  return pair<FLOAT_T,FLOAT_T>(get_double_parameter("min-mass"),
                               get_double_parameter("max-mass"));
}

/**
 * Use parameter.cpp for window type and window size to find the mass
 * window around the given m/z or mass.
 * \returns A pair with the minimum and maximum masses to search.
 */
pair<FLOAT_T,FLOAT_T> ModifiedPeptidesIterator::getMinMaxMass(
  double mz, ///< precursor mz for peptide window
  SpectrumZState& zstate, ///< charge/mass pair for peptide window
  PEPTIDE_MOD_T* pmod) ///< peptide mod with the delta mass for peptides
{
  WINDOW_TYPE_T precursor_window_type = 
    get_window_type_parameter("precursor-window-type");
  double window = get_double_parameter("precursor-window");
  double min_mass = 0;
  double max_mass = 0;

  // get the mass difference
  double delta_mass = peptide_mod_get_mass_change(pmod);

  if (precursor_window_type == WINDOW_MASS) {
    double mass = zstate.getNeutralMass() - delta_mass;
    min_mass = mass - window;
    max_mass = mass + window;
  } else if (precursor_window_type == WINDOW_MZ) {
    double min_mz = mz - window;
    double max_mz = mz + window;
    min_mass = (min_mz - MASS_PROTON) * (double)zstate.getCharge() - delta_mass;
    max_mass = (max_mz - MASS_PROTON) * (double)zstate.getCharge() - delta_mass;
  } else if (precursor_window_type == WINDOW_PPM) {
    double mass = zstate.getNeutralMass() - delta_mass;
    min_mass = mass / (1.0 + window * 1e-6);
    max_mass = mass / (1.0 - window * 1e-6);
    carp(CARP_DEBUG,"mass:%f charge:%i min_mass:%f max_mass:%f",
         mass, zstate.getCharge(), min_mass, max_mass);
  } else {
    carp(CARP_FATAL,"Invalid window type");
  }
  pair<FLOAT_T,FLOAT_T> min_max(min_mass, max_mass);
  return min_max;
}

/**
 * \brief Queue next_peptide to be returned.
 *
 * Takes unmodified peptides from a base class and applies the
 * modification to each.
 * 
 * Filters out peptides that can't be modified. Stores multiple
 * versions of same peptide that can have the mod applied in more than
 * one way (ie on different residues) in temp_peptide_list_.  Deletes
 * elements from list as used. Sets next_peptide to NULL when there are
 * no more to be returned.
 * \returns True if there is a peptide to return or false if there are
 * none.
 */
bool ModifiedPeptidesIterator::queueNextPeptide(){

  carp(CARP_DETAILED_DEBUG, "ModifiedPeptidesIterator queueNextPeptide");

  // first, try getting next from the temp list
  if( ! is_empty_linked_list( temp_peptide_list_)){
    carp(CARP_DETAILED_DEBUG,"Queue is getting next peptide from temp list");
    next_peptide_ = (Peptide*)pop_front_linked_list(temp_peptide_list_);

    return true;
  }

  // second, try getting next from iterator
  if( ! peptide_source_->hasNext() ){
    carp(CARP_DETAILED_DEBUG, 
        "ModifiedPeptidesIterator queue has no more peptides in the generator");
    next_peptide_ = NULL;
    return false;   // no more peptides for this iterator
  }

  // else, get the next unmodified peptide
  Peptide* unmod_peptide = peptide_source_->next();
  
  IF_CARP_DETAILED_DEBUG(
    char* debugseq = unmod_peptide->getSequence();
    carp(CARP_DETAILED_DEBUG, "Next peptide in pep_gen is %s", debugseq);
    free(debugseq);
  )

  // apply modifications, discard peptides that can't be modified

  // keep looking until a peptide can be modified or we run out of peptides
  carp(CARP_DETAILED_DEBUG, "Queue is looking for modifyable peptide");
  while( unmod_peptide != NULL &&
         ! is_peptide_modifiable(unmod_peptide, peptide_modification_) ){ 
    delete unmod_peptide;
    unmod_peptide = peptide_source_->next();
  }

  if( unmod_peptide == NULL ){ 
    // none of the remaining peptides were modifiable
    carp(CARP_DETAILED_DEBUG, "Skipped all remaining peptides in generator");
    next_peptide_ = NULL;
    return false;
  }

  IF_CARP_DETAILED_DEBUG(
    char* umodseq = unmod_peptide->getSequence();
    carp(CARP_DETAILED_DEBUG, "Iterator is modifying peptide %s",
         umodseq);
    free(umodseq);
  )
  modify_peptide(unmod_peptide, 
                 peptide_modification_, 
                 temp_peptide_list_,
                 max_aas_modified_ );
  // this put a copy in the list, get rid of the original
  delete unmod_peptide;

  if( is_empty_linked_list(temp_peptide_list_) ){
    carp(CARP_DETAILED_DEBUG, "Modifier didn't return any peptides");
    next_peptide_ = NULL;
    return false;
  }

  // now set next_peptide to the first in the list and move list forward
  next_peptide_ = (Peptide*)pop_front_linked_list(temp_peptide_list_);
  
  IF_CARP_DETAILED_DEBUG(
    char* seq = next_peptide_->getModifiedSequenceWithMasses(MOD_MASS_ONLY);
    carp(CARP_DETAILED_DEBUG, "Queue set next peptide as %s", seq);
    free(seq);
  )

  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

