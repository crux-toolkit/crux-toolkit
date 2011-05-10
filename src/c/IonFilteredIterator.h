/**
 * \file IonFilteredIterator.h
 ****************************************************************************/
#ifndef IONFILTEREDITERATOR_H
#define IONFILTEREDITERATOR_H

#include "IonSeries.h"
#include "IonConstraint.h"


/**
 * \class ion_filtered_iterator
 * \brief An object to iterate over ion objects that meet constraint in
 * the ion_series 
 */
class IonFilteredIterator{
 protected:
  IonSeries* ion_series_; ///< the ion series that the ion we are iterating
  IonConstraint* constraint_; ///< constraints which the ions obey
  bool has_next_; ///< the boolean which the iterator has a next ion
  int ion_idx_; ///< the current ion that is being returned 
  Ion* ion_; ///< the next ion to return when called upon

  IonIterator ion_iterator_;
  IonIterator ion_iterator_end_;
  /**
   * sets up the iterator for next iteration.
   * 
   *\returns true if successfully sets up the ion_filtered_iterator for next iteration
   */
  bool setup();

 public:

  /**
   * Only copies in the constraint as pointer
   * Instantiates a new ion_filtered_iterator object from ion_series.
   * \returns an IonFilteredIterator object.
   */
  IonFilteredIterator(
    IonSeries* ion_series, ///< ion_series to iterate -in
    IonConstraint* constraint  ///< ion_constraint which returned ions satisfy
    );
    
  /**
   * The constraint is NOT freed from the iterator.
   * Frees an allocated ion_filtered_iterator object.
   */
  virtual ~IonFilteredIterator();

  /**
   * The basic iterator function has_next.
   */
  bool hasNext();
  
  /**
   * The basic iterator function next.
   */
  Ion* next();

};


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
