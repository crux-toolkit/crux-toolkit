/*************************************************************************//**
 * \file FilteredSpectrumChargeIterator.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * \brief code to support working with collection of multiple spectra
 ****************************************************************************/
#include "FilteredSpectrumChargeIterator.h"

#include "parameter.h"

using namespace std;
using namespace Crux;

/**
 * \brief Sets up an iterator with the next spectrum that complies
 * with the constraints.  Sets has_next to FALSE when there are no
 * more spectra in the collection that pass.  Increments
 * spectrum_index and charge_index.
 */
void FilteredSpectrumChargeIterator::queueNextSpectrum() {

  Spectrum* spec = NULL;

  // Are there any more charge states for this spectrum?
  if( zstate_index_ < (int)zstates_.size()-1 ){
    zstate_index_++;
  }
  // Are there any more spectra?
  else if( spectrum_index_ < spectrum_collection_->getNumSpectra()-1){

    spectrum_index_++;
    spec = spectrum_collection_->spectra_[spectrum_index_];
    // first free any existing zstates in the iterator
    zstates_.clear();
    zstates_ = spec->getZStatesToSearch();
    zstate_index_ = 0;
  }else{ // none left
    has_next_ = false;
    return;
  }

  // Does the current pass?
  spec = spectrum_collection_->spectra_[spectrum_index_];
  int this_charge = -1;
  if (zstate_index_ < (int)zstates_.size()) {
    this_charge = zstates_[zstate_index_].getCharge();
  }
  double mz = spec->getPrecursorMz();
  int num_peaks = spec->getNumPeaks();

  // Warn if we skip due to too few peaks.
  if ( num_peaks < min_peaks_ ) {
    carp(CARP_INFO, "Skipping scan %d with %d (< %d) peaks.",
	 spectrum_index_, num_peaks, min_peaks_);
  }

  if( search_charge_ == 0 || search_charge_ == this_charge ){
    if( mz >= min_mz_ && mz <= max_mz_
        && num_peaks >= min_peaks_ ){
      // passes all tests
      has_next_ = true;
      return;
    }
  }
  
  // try the next spectrum
  queueNextSpectrum();

}


 /**
 * Instantiates a new spectrum_iterator object from
 * spectrum_collection.  This iterator returns unique spectrum-charge
 * pairs (e.g.a spectrum to be searched as +2 and +3 is returned once
 * as +2 and once as +3).  The charge is returned by setting the int
 * pointer in the argument list.  The iterator also filters spectra by
 * mass so that none outside the spectrum-min-mass--spectrum-max-mass
 * range (as defined in parameter.c).  The iterator also filters by
 * minimum number of peaks.
 * \returns a SPECTRUM_ITERATOR_T object.
 */
FilteredSpectrumChargeIterator::FilteredSpectrumChargeIterator(
  SpectrumCollection* spectrum_collection
) {

  spectrum_collection_ = spectrum_collection;  
  has_next_ = false;
  spectrum_index_ = -1;
  zstate_index_ = -1;
  min_mz_ = get_double_parameter("spectrum-min-mass");
  max_mz_ = get_double_parameter("spectrum-max-mass");
  min_peaks_ = get_int_parameter("min-peaks");

  const char* charge_str = get_string_parameter_pointer("spectrum-charge");
  if( strcmp( charge_str, "all") == 0){
    search_charge_ = 0;
  }else{
    search_charge_ = atoi(charge_str);
  }

  // queue next spectrum
  carp(CARP_DEBUG,"Queueing next spectrum");
  queueNextSpectrum();
  carp(CARP_DEBUG,"Done Queueing next spectrum");

}
  
FilteredSpectrumChargeIterator::~FilteredSpectrumChargeIterator() {

}

bool FilteredSpectrumChargeIterator::hasNext() {
  return has_next_;
}

Spectrum* FilteredSpectrumChargeIterator::next(SpectrumZState& zstate) {

  carp(CARP_DEBUG,"FilteredSpectrumChargeIterator::next()");
  Spectrum* next_spectrum = NULL;
  if (has_next_) {
    next_spectrum = 
      spectrum_collection_->spectra_[spectrum_index_];

    zstate = zstates_[zstate_index_];
    carp(CARP_DEBUG,"Queueing next spectrum");
    queueNextSpectrum();
    carp(CARP_DEBUG,"Done queueing next spectrum");
  }
  return next_spectrum;

}

  




