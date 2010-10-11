/**
 * \file PeakIterator.cpp
 * DATE: September 16, 2010
 * AUTHOR: Barbara Frewen
 * \brief Object for returning peaks in a Spectrum.
 */

#include "PeakIterator.h"

using namespace std;

/**
 * Constructor that initializes the PeakIterator to a particular
 * Spectrum.
 */
PeakIterator::PeakIterator(Spectrum* spectrum) : peak_index_(0)
{
  this->spectrum_ = spectrum;
}

/**
 * Default destructor.  Does not free the Spectrum
 */
PeakIterator::~PeakIterator(){}

/**
 * \returns TRUE if there are additional peaks to iterate over,
 * FALSE if not. 
 */
bool PeakIterator::has_next()
{
  return (peak_index_ < spectrum_->getNumPeaks());
}

/**
 * \returns The next peak object in the spectrum, in order of m/z.
 */
PEAK_T* PeakIterator::next()
{
  if( peak_index_ > spectrum_->getNumPeaks() ){
    return NULL;
  }
  PEAK_T* next_peak = spectrum_->peaks_[peak_index_];
  peak_index_++;
  return next_peak;
}

/**
 *  Resets the iterator to the first element
 */
void PeakIterator::reset()
{
  peak_index_ = 0; 
}
