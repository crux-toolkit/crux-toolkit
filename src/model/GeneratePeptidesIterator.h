/**
 * \file GeneratePeptidesIterator.h
 * \brief An object to return candidate peptides from a database.
 */
#ifndef GENERATEPEPTIDESITERATOR_H 
#define GENERATEPEPTIDESITERATOR_H 

#include <utility>
#include "PeptideIterator.h"
#include "Database.h"
#include "PeptideConstraint.h"

class GeneratePeptidesIterator : public PeptideIterator
{
 protected:
  PeptideIterator* iterator_;

  virtual bool queueNextPeptide();

 public:
  GeneratePeptidesIterator();

  GeneratePeptidesIterator(
    std::pair<FLOAT_T,FLOAT_T> min_max_mass, ///< max target mass of peptides
    bool is_decoy,  ///< generate target or decoy peptides
    Database* database,  ///< database to provide peptides
    int additional_missed_cleavages = 0
  );

  ~GeneratePeptidesIterator();

};

#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

