/**
 * \file PWIZSpectrumCollection.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 21 April 2011
 * \brief Class to read spectra files using the proteowizard library.
 */
#include "PWIZSpectrumCollection.h" 
#include "crux-utils.h"
#include "parameter.h"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"

/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does not parse file. 
 */
PWIZSpectrumCollection::PWIZSpectrumCollection(
  const char* filename   ///< The spectrum collection filename.
 ) : SpectrumCollection(filename){

}

/**
 * Parses all the spectra from file designated by the filename member
 * variable.
 * \returns True if the spectra are parsed successfully. False if otherwise.
 */
bool PWIZSpectrumCollection::parse() {

  // spectrum_collection has already been parsed
  if(is_parsed_){
    return false;
  }

  pwiz::msdata::MSDataFile* reader = new pwiz::msdata::MSDataFile(filename_);
  if( reader == NULL ){
    carp(CARP_FATAL, "PWIZSpectrumCollection unable to open '%s'.", 
         filename_.c_str());
  }
  carp(CARP_DEBUG, "Using proteowizard to parse spectra.");

  // todo see getFilters() below
  //  pwiz::analysis::SpectrumListFactory::wrap(*reader, spectrum_filters);

  // TODO add first/last scan to base class
  // get a list of scans to include if requested
  const char* range_string = get_string_parameter("scan-number");
  int first_scan = -1;
  int last_scan = -1;

  get_range_from_string(
    range_string,
    first_scan,
    last_scan);

  if( first_scan == -1 || last_scan == -1 ){
    carp(CARP_FATAL, "The scan number range '%s' is invalid. "
         "Must be of the form <first>-<last>.", range_string);
  }

  // get info for translating identifiers into scan numbers
  pwiz::msdata::CVID native_id_format = pwiz::msdata::id::getDefaultNativeIDFormat(*reader);
  native_id_format = native_id_format;

  // look at all spectra in file
  pwiz::msdata::SpectrumListPtr all_spectra = reader->run.spectrumListPtr;
  bool get_peaks = true;
  
  int num_spec = all_spectra->size();
  carp(CARP_INFO, "PWIZ:Number of spectra:%i",num_spec);
  for(int spec_idx = 0; spec_idx < num_spec; spec_idx++){
    carp(CARP_DETAILED_DEBUG, "Parsing spectrum index %d.", spec_idx);
    pwiz::msdata::SpectrumPtr spectrum = all_spectra->spectrum(spec_idx, 
                                                               get_peaks);

    // skip if not ms2
    if( spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>() != 2 ){
      continue;
    }

    // check that scan number is in range
    int scan_number = pwiz::msdata::id::valueAs<int>(spectrum->id, "scan");
    carp(CARP_INFO,"found scan:%i",scan_number);
    if( scan_number < first_scan ){
      continue;
    } else if( scan_number > last_scan ){
      break;
    }

    Spectrum* crux_spectrum = new Spectrum();
    crux_spectrum->parsePwizSpecInfo(spectrum);

    this->addSpectrumToEnd(crux_spectrum);
    //?delete spectrum?;
  }

  delete reader;

  return true;
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.  Removes any existing information in
 * the given spectrum. 
 * \returns True if the spectrum was allocated, false on error.
 */
bool PWIZSpectrumCollection::getSpectrum(
  int first_scan,      ///< The first scan of the spectrum to retrieve -in
  Spectrum* spectrum   ///< Put the spectrum info here
  )
{
  (void)first_scan;
  (void)spectrum;
  carp(CARP_DEBUG, "Using proteowizard to parse spectrum");
  /*
  MSToolkit::MSReader* mst_reader = new MSToolkit::MSReader();
  MSToolkit::Spectrum* mst_spectrum = new MSToolkit::Spectrum();
  bool parsed = false;

  mst_reader->readFile(filename_.c_str(), *mst_spectrum, first_scan);

  if(mst_spectrum->getScanNumber() != 0) {
    spectrum->parseMstoolkitSpectrum(mst_spectrum,
                                     filename_.c_str());
    parsed = true;
  } else {
    carp(CARP_ERROR, "Spectrum %d does not exist in file", first_scan);
    parsed = false;
  }
  delete mst_spectrum;
  delete mst_reader;
  return parsed;
  */
  return false;
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.
 * \returns The spectrum data from file or NULL.
 */
Spectrum* PWIZSpectrumCollection::getSpectrum(
  int first_scan      ///< The first scan of the spectrum to retrieve -in
  )
{
  (void)first_scan;
  carp(CARP_DEBUG, "Using proteowizard to parse spectrum");
  /*
  MSToolkit::MSReader* mst_reader = new MSToolkit::MSReader();
  MSToolkit::Spectrum* mst_spectrum = new MSToolkit::Spectrum();
  Spectrum* return_spec = NULL;  

  mst_reader->readFile(filename_.c_str(), *mst_spectrum, first_scan);

  if(mst_spectrum->getScanNumber() != 0) {
    return_spec = new Spectrum();
    return_spec->parseMstoolkitSpectrum(mst_spectrum,
                                        filename_.c_str()); 
  } else {
    carp(CARP_ERROR, "Spectrum %d does not exist in file", first_scan);
    return_spec = NULL;
  }
  delete mst_spectrum;
  delete mst_reader;
  return return_spec;
  */
  return NULL;
}

/*
void getFilters(){  
  // select only the MS2 level spectra and only those in the given scan range
  ostringstream string_builder;
  vector<string> spectrum_filters;
  spectrum_filters.push_back("msLevel 2");

  string_builder << "scanNumber [";
  string_builder << first_scan_;
  string_builder << ",";
  string_builder << last_scan_;
  string_buildter << "]";
  spectrum_filters.push_back(string_builder.c_str());
  // where filters include scan range, ms2 only
}
*/
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
