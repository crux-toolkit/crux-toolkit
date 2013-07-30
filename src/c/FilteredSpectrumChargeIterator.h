/**
 * \file FilteredSpectrumChargeIterator.h 
 * \brief Object for representing many spectra.
 *****************************************************************************/
#ifndef FILTEREDSPECTRUMCHARGEITERATOR_H
#define FILTEREDSPECTRUMCHARGEITERATOR_H

//#include "objects.h"
#include <vector>

#include "Spectrum.h"
#include "SpectrumCollection.h"

class FilteredSpectrumChargeIterator {
 protected:
  Crux::SpectrumCollection* spectrum_collection_;///< spectra to iterate over
  bool has_next_;  ///< is there a spec that passes criteria
  int spectrum_index_; ///< The index of the current spectrum
  std::vector<SpectrumZState> zstates_;        ///< Array of possible zstates to search
  int zstate_index_;    ///< The index of the z-state for the current spectrum
  double min_mz_;       ///< return only spec above this mz
  double max_mz_;      ///< return only spec below this mz
  int search_charge_;   ///< which z to search, 0 for all
  int min_peaks_;       ///< minimum number of peaks a spec must have

  /**
   * \brief Sets up an iterator with the next spectrum that complies
   * with the constraints.  Sets has_next to FALSE when there are no
   * more spectra in the collection that pass.  Increments
   * spectrum_index and charge_index.
   */
  void queueNextSpectrum();

 public:
 /**
 * Instantiates a new spectrum_iterator object from
 * spectrum_collection.  This iterator returns unique spectrum-charge
 * pairs (e.g.a spectrum to be searched as +2 and +3 is returned once
 * as +2 and once as +3).  The charge is returned by setting the int
 * pointer in the argument list.  The iterator also filters spectra by
 * mass so that none outside the spectrum-min-mz--spectrum-max-mz
 * range (as defined in parameter.c).  The iterator also filters by
 * minimum number of peaks.
 * \returns a SPECTRUM_ITERATOR_T object.
 */
  FilteredSpectrumChargeIterator(
    Crux::SpectrumCollection* spectrum_collection
  );
  
  /**
   * Frees an filtered_spectrum_charge_iterator object.
   */
  ~FilteredSpectrumChargeIterator();

  /**
   * The basic iterator function has_next.
   */
  bool hasNext();

  /**
   * The basic iterator function next.
   */
  Crux::Spectrum* next(SpectrumZState& zstate);

};




/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
