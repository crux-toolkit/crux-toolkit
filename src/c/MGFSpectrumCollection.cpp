/**
 * \file MGFSpectrumCollection.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 April 2011
 * \brief Class for accessing spectra in .mgf (Mascot general format) files.
 */
#include "MGFSpectrumCollection.h"
#include "parameter.h"

/**
 * Constructor sets filename and initializes member variables.
 */
MGFSpectrumCollection::MGFSpectrumCollection(
 const char* filename ///< The spectrum collection filename. 
) : SpectrumCollection(filename), cur_spec_number_(0)
{}

/**
 * Parses all the spectra from file designated by the filename member
 * variable.
 * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
 */
bool MGFSpectrumCollection::parse(){
  if(is_parsed_){
    return false;
  }

  // get a list of scans to include if requested
  const char* range_string = get_string_parameter("scan-number");
  int first_scan = get_first_in_range_string(range_string);
  int last_scan = get_last_in_range_string(range_string);
  if( first_scan == -1 || last_scan == -1 ){
    carp(CARP_FATAL, "The scan number range '%s' is invalid. "
         "Must be of the form <first>-<last>.", range_string);
  }
  
  // open file
  FILE* file = fopen(filename_.c_str(), "r");
  if( file == NULL ){
    carp(CARP_ERROR, "File %s could not be opened",
         filename_.c_str());
    return false;
  }

  // read the first spectrum
  Spectrum* spectrum = Spectrum::newSpectrumMgf(file, ++cur_spec_number_, 
                                                filename_.c_str() );
  while(spectrum){
    // is this a scan to include? if not skip it
    if( spectrum->getFirstScan() < first_scan ){
      delete spectrum;
      spectrum = Spectrum::newSpectrumMgf(file, ++cur_spec_number_, 
                                          filename_.c_str());
      continue;
    } 
    // are we past the last scan?
    if( spectrum->getFirstScan() > last_scan ){
      delete spectrum;
      break;
    }

    this->addSpectrumToEnd(spectrum); // keep this spectrum

    // get next
    spectrum = Spectrum::newSpectrumMgf(file, ++cur_spec_number_,
                                        filename_.c_str());

  }

  fclose(file);
  is_parsed_ = true;

  return true;
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.
 * \returns The newly allocated Spectrum or NULL if scan number not found.
 */
Spectrum* MGFSpectrumCollection::getSpectrum(
  int first_scan      ///< The first scan of the spectrum to retrieve -in
){
  (void)first_scan;
  carp(CARP_ERROR, "Cannot select spectrum by scan number for MGF files.");
  return NULL;
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.  Removes any existing information in the
 * given spectrum.
 * \returns True if the spectrum was allocated, false on error.
 */
bool MGFSpectrumCollection::getSpectrum(
  int first_scan,      ///< The first scan of the spectrum to retrieve -in
  Spectrum* spectrum   ///< Put the spectrum info here
){
  (void)first_scan;
  (void)spectrum;
  carp(CARP_ERROR, "Cannot select spectrum by scan number for MGF files.");
  return false;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
