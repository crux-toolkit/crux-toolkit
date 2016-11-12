/**
 * \file GeneratePeptidesIterator.h
 * \brief An object to return candidate peptides from a database.
 */
#include <iostream>
#include "GeneratePeptidesIterator.h"
#include "DatabasePeptideIterator.h"

using namespace std;

GeneratePeptidesIterator::GeneratePeptidesIterator()
{}

GeneratePeptidesIterator::GeneratePeptidesIterator(
  pair<FLOAT_T,FLOAT_T> min_max_mass, ///< precursor m/z of spectrum
  bool is_decoy,  ///< generate target or decoy peptides
  Database* database, ///< database to provide peptides
  int additional_missed_cleavages
  )
{
  PeptideConstraint* constraint = PeptideConstraint::newFromParameters();
  constraint->setMinMass(min_max_mass.first);
  constraint->setMaxMass(min_max_mass.second);
  constraint->setNumMisCleavage(constraint->getNumMisCleavage()+additional_missed_cleavages);



  // Check that database exists
  if (database == NULL){
    carp(CARP_FATAL, "Cannot generate peptides when database is NULL.");
  }

	iterator_ = new DatabasePeptideIterator(database,
																					constraint,
																					true, // store all peptides
																					is_decoy);

  initialize();
  PeptideConstraint::free(constraint);
}

GeneratePeptidesIterator::~GeneratePeptidesIterator()
{
  delete iterator_;
}

bool GeneratePeptidesIterator::queueNextPeptide(){
  if( iterator_ == NULL ){
    carp(CARP_FATAL, 
         "GeneratePeptidesIterator is missing its database iterator.");
  }

  if( iterator_->hasNext() ){
    carp(CARP_DETAILED_DEBUG, 
         "GeneratePeptidesIterator has a peptide, returning true");
    next_peptide_ = iterator_->next();
    return true;
  }
  // else

  next_peptide_ = NULL;
  return false;
}
