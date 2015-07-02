/**
 * \file IonFilteredIterator.cpp 
 * \brief An object to iterate over ion objects that meet constraints
 * in the IonSeries.
 *****************************************************************************/
#include "IonFilteredIterator.h"

/**
 * sets up the iterator for next iteration.
 * 
 *\returns true if successfully sets up the ion_filtered_iterator for next iteration
 */
bool IonFilteredIterator::setup()
{
  bool found = false;
  Ion* current_ion = NULL;

  // iterate over ions until discovers the first ion that meets the ion constraint
  while((!found) && (ion_iterator_ != ion_iterator_end_)){
    // get next ion
    current_ion = *ion_iterator_;
    
    // check if the current ion satisfies the ion_constraint for the iterator
    if(constraint_->isSatisfied(current_ion)){
      found = true;
      ion_ = current_ion;
    }
    ++ion_iterator_;
  }
  
  has_next_ = found;

  return true;
}

/**
 * Only copies in the constraint as pointer
 * Instantiates a new ion_filtered_iterator object from ion_series.
 * \returns a IonFilteredIterator object.
 */
IonFilteredIterator::IonFilteredIterator(
  IonSeries* ion_series, ///< ion_series to iterate -in
  IonConstraint* constraint  ///< ion_constraint which returned ions satisfy
  )
{

  
  // set constraint, ion_series
  constraint_ = constraint;
  ion_series_ = ion_series;
  has_next_ = false;

  // set the working array of ions
  if(constraint_->getIonType() == ALL_ION ||
     constraint_->getIonType() == BY_ION ||
     constraint_->getIonType() == BYA_ION){

    ion_iterator_ = ion_series->begin();
    ion_iterator_end_ = ion_series->end();
  }
  else{

    ion_iterator_ = ion_series->getSpecificIons(constraint->getIonType()).begin();
    ion_iterator_end_ = ion_series->getSpecificIons(constraint->getIonType()).end();
  }
  
  // initialize iterator
  setup();

}        

/**
 * The constraint is NOT freed from the iterator.
 * Frees an allocated ion_filtered_iterator object.
 */
IonFilteredIterator::~IonFilteredIterator() 
{

}

/**
 * The basic iterator function has_next.
 */
bool IonFilteredIterator::hasNext()
{
  return has_next_;
}

/**
 * The basic iterator function next.
 */
Ion* IonFilteredIterator::next()
{

  Ion* next_ion = NULL;
  
  // check if a ion is present to return
  if(!has_next_){
    carp(CARP_FATAL, "index out of bounds for ion_filtered_iterator");
  }
  
  next_ion = ion_;
  
  // re-initialize iterator
  setup();

  return next_ion;
}
