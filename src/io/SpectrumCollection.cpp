/**
 * \file SpectrumCollection.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 11 April 2011
 * \brief Abstract class for accessing spectra from a file.
 */
#include "SpectrumCollection.h" 
#include "model/ProteinIndex.h" 
#include "model/Peak.h"
#include "util/utils.h"
#ifndef _MSC_VER
#include "unistd.h"
#endif
#include "parameter.h"
#include <cerrno>
#include <cstring>
#include "io/carp.h"
#include "util/WinCrux.h"
#include <iostream>

using namespace std;
using namespace Crux;

namespace Crux {

/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Resolves any relative paths.  Confirms that file exists.
 */
SpectrumCollection::SpectrumCollection (
  const string& filename ///< The spectrum collection filename. 
  ) 
: filename_(filename), is_parsed_(false), num_charged_spectra_(0) {
#if DARWIN
  char path_buffer[PATH_MAX];
  char* absolute_path_file =  realpath(filename.c_str(), path_buffer);
#else
  char* absolute_path_file =  realpath(filename.c_str(), NULL);
#endif
  if (absolute_path_file == NULL) {
    carp(CARP_FATAL, "Error from spectrum file '%s'. (%s)",
         filename.c_str(), strerror(errno)); 
  }
  
  if(access(absolute_path_file, F_OK)) {
    carp(CARP_FATAL, "File %s could not be opened\n", absolute_path_file);
  }
  filename_ = absolute_path_file;

#ifndef DARWIN
  free(absolute_path_file);
#endif
}

/**
 * Copy constructor.  Creates copies of each spectrum for the new
 * collection.
 */
SpectrumCollection::SpectrumCollection(
  SpectrumCollection& old_collection
  ) : filename_(old_collection.filename_),
      is_parsed_(old_collection.is_parsed_),
      num_charged_spectra_(old_collection.num_charged_spectra_) {
  // copy spectra
  for (SpectrumIterator spectrum_iterator = old_collection.begin();
    spectrum_iterator != old_collection.end();
    ++spectrum_iterator) {

    Spectrum* old_spectrum = *spectrum_iterator;
    Spectrum* new_spectrum = new Spectrum(*old_spectrum);
    this->addSpectrumToEnd(new_spectrum);
  }
} 

/**
 * Default destructor
 */
SpectrumCollection::~SpectrumCollection() {

  for (SpectrumIterator spectrum_iterator = this->begin();
    spectrum_iterator != this->end();
    ++spectrum_iterator) {
    delete *spectrum_iterator;    
  }
  spectra_.clear();
}  

/**
 * \returns the begining of the spectra deque
 */
SpectrumIterator SpectrumCollection::begin() {
  return spectra_.begin();
}

/**
 * \returns the end of the spectra deque
 */
SpectrumIterator SpectrumCollection::end() {
  return spectra_.end();
}

/*
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.
 * \returns The newly allocated Spectrum or NULL if scan number not found.
 */
Spectrum* SpectrumCollection::getSpectrum(
  int first_scan      ///< The first scan of the spectrum to retrieve -in
) {
  Spectrum* spectrum = new Spectrum();
  bool success = getSpectrum(first_scan, spectrum);
  if (success) {
    return spectrum;
  } else {
    delete spectrum;
    return NULL;
  }
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.  Removes any existing information in the
 * given spectrum.
 * \returns True if the spectrum was allocated, false on error.
 */
bool SpectrumCollection::getSpectrum(
  int first_scan,      ///< The first scan of the spectrum to retrieve -in
  Spectrum* spectrum   ///< Put the spectrum info here
) {
  map<int, Spectrum*>::const_iterator i = spectraByScan_.find(first_scan);
  if (i != spectraByScan_.end()) {
    spectrum->copyFrom(i->second);
    return true;
  }
  for (SpectrumIterator spectrum_iterator = this->begin();
       spectrum_iterator != this->end();
       ++spectrum_iterator) {
    if ((*spectrum_iterator)->getFirstScan() == first_scan) {
      spectrum->copyFrom(*spectrum_iterator);
      return true;
    }
  }
  return false;
}


/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum to the end of the spectra array
 * should only be used when the adding in increasing scan num order
 * when adding in random order should use add_spectrum
 * spectrum must be heap allocated
 */
void SpectrumCollection::addSpectrumToEnd(
  Spectrum* spectrum ///< spectrum to add to spectrum_collection -in
  ) {
  // set spectrum
  spectra_.push_back(spectrum);
  num_charged_spectra_ += spectrum->getNumZStates();
}

/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum in correct order into the spectra array
 * spectrum must be heap allocated
 */
void SpectrumCollection::addSpectrum(
  Spectrum* spectrum ///< spectrum to add to spectrum_collection -in
  ) {
    
  unsigned int add_index = 0;

  // find correct location
  // TODO -- replace with binary search if necessary.
  for(; add_index < spectra_.size(); ++add_index) {
    if(spectra_[add_index]->getFirstScan() >
       spectrum->getFirstScan()) {
      break;
    }
  }

  spectra_.insert(spectra_.begin()+add_index, spectrum);

  num_charged_spectra_ += spectrum->getNumZStates();
}


// FIXME maybe a faster way? can't perform binary search since we must know the array index
/**
 * Removes a spectrum from the spectrum_collection.
 */
void SpectrumCollection::removeSpectrum(
  Spectrum* spectrum ///< spectrum to be removed from spectrum_collection -in
) {
  int scan_num = spectrum->getFirstScan();
  unsigned int spectrum_index = 0;
  
  // find where the spectrum is located in the spectrum array
  for(; spectrum_index < spectra_.size(); ++spectrum_index) {
    if(scan_num == spectra_[spectrum_index]->getFirstScan() ) {
      break;
    }
  }
  
  num_charged_spectra_ -= spectrum->getNumZStates();

  delete spectra_[spectrum_index];
  spectra_[spectrum_index] = NULL;
  spectra_.erase(spectra_.begin() + spectrum_index);
} 

/**
 * \returns A pointer to the name of the file from which the spectra
 * were parsed.
 */
const char* SpectrumCollection::getFilename() {
  return filename_.c_str();
}

/**
 * \returns The current number of spectrum in the
 * spectrum_collection.  Zero if the file has not yet been parsed.
 */
int SpectrumCollection::getNumSpectra() {
  return spectra_.size();
}

/**
 * \returns The current number of spectra assuming differnt
 * charge(i.e. one spectrum with two charge states are counted as two
 * spectra) in the spectrum_collection.  Zero if the file has not been
 * parsed.
 */
int SpectrumCollection::getNumChargedSpectra() {
  return num_charged_spectra_;
}


/**
 * \returns True if the spectrum_collection file has been parsed.
 */
bool SpectrumCollection::getIsParsed() {
  return is_parsed_;
}

} // namespace Crux

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
