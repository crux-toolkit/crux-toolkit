/**
 * \file GenerateModifiedPeptidesIterator.h
 * \brief An object to return candidate peptides from a database or index.
 */
#ifndef MODIFIEDPEPTIDESITERATOR_H 
#define MODIFIEDPEPTIDESITERATOR_H 

#include <vector>
#include "GeneratePeptidesIterator.h"
#include "peptide_modifications.h"
#include "SpectrumZState.h"

class ModifiedPeptidesIterator : public PeptideIterator
{
 protected:
  GeneratePeptidesIterator* peptide_source_;///< gets unmod peptides from index
  PEPTIDE_MOD_T* peptide_modification_; ///< the modification to apply
  LINKED_LIST_T* temp_peptide_list_;///< modified peptides ready to return
  int max_aas_modified_;///< modify no more than this many aas per peptides

  std::pair<FLOAT_T,FLOAT_T> getMinMaxMass();
  std::pair<FLOAT_T,FLOAT_T> getMinMaxMass(
    double mz, 
    SpectrumZState& zstate, 
    PEPTIDE_MOD_T* pmod);

  virtual bool queueNextPeptide();

 public:
  /**
   * Constructor for returning peptides in a window around the given
   * m/z.
   */
  ModifiedPeptidesIterator(
    double mz,               ///< Spectrum precrusor mz
    SpectrumZState& zstate,  ///< Target mz of peptides
    PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
    bool is_decoy,  ///< generate decoy peptides
    Index* index,      ///< Index from which to draw peptides OR
    Database* dbase    ///< Database from which to draw peptides
  );

  /**
   * Constructor for returning all peptides in the index or database
   * that fall within the constraints defined in parameter.cpp.
   */
  ModifiedPeptidesIterator(
    PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
    Index* index,      ///< Index from which to draw peptides OR
    Database* dbase    ///< Database from which to draw peptides
   );

  /**
   * Constructor for returnign all peptides in the index or database
   * that fall within the mass range
   */
  ModifiedPeptidesIterator(
    double min_mass,    ///< min-mass of peptides
    double max_mass,    ///< max-mass of peptides
    PEPTIDE_MOD_T* pmod, ///< Peptide mod to apply
    bool is_decoy, ///< generate decoy peptides
    Index* index,     ///< Index from which to draw peptides OR
    Database* dbase   ///< Database from which to draw peptides
  );


  ~ModifiedPeptidesIterator();
};

#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

