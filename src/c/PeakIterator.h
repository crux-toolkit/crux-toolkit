/**
 * \file PeakIterator.h 
 * $Revision: $
 * \brief Object for returning peaks in a spectrum.
 *****************************************************************************/
#ifndef PEAK_ITERATOR_H
#define PEAK_ITERATOR_H 

#include "Spectrum.h"

/**
 * \class PeakIterator
 * \brief Object to iterate over the peaks in a spectrum.
 */
class PeakIterator {
 protected:
  Spectrum* spectrum_; ///< The spectrum whose peaks to iterate over. 
  int  peak_index_;    ///< The index of the current peak

 public:
  /**
   * Constructor that initializes the PeakIterator to a particular
   * Spectrum.
   */
  PeakIterator(Spectrum* spectrum);

  /**
   * Default destructor.
   */
  ~PeakIterator();

  /**
   * \returns TRUE if there are additional peaks to iterate over,
   * FALSE if not. 
   */
  bool has_next();

  /**
   * \returns The next peak object in the spectrum, in order of m/z.
   */
  PEAK_T* next();

  /**
   *  Resets the iterator to the first element
   */
  void reset();

};

#endif // PEAK_ITERATOR_H 
