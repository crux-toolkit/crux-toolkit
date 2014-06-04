/**
 * \file PWIZSpectrumCollection.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 21 April 2011
 * \brief Class to read spectra files using the proteowizard library.
 */
#include "PWIZSpectrumCollection.h" 
#include "crux-utils.h"
#include "parameter.h"
#include <iostream>
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#ifdef _MSC_VER
#include "pwiz/data/msdata/DefaultReaderList.hpp"
//#include "pwiz/data/vendor_readers/ABI/Reader_ABI.hpp"
//#include "pwiz/data/vendor_readers/ABI/T2D/Reader_ABI_T2D.hpp"
#include "pwiz/data/vendor_readers/Agilent/Reader_Agilent.hpp"
#include "pwiz/data/vendor_readers/Bruker/Reader_Bruker.hpp"
//#include "pwiz/data/vendor_readers/Shimadzu/Reader_Shimadzu.hpp"
#include "pwiz/data/vendor_readers/Thermo/Reader_Thermo.hpp"
#include "pwiz/data/vendor_readers/Waters/Reader_Waters.hpp"
#endif


using namespace std;
/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does not parse file. 
 */
PWIZSpectrumCollection::PWIZSpectrumCollection(
  const char* filename   ///< The spectrum collection filename.
 ) : SpectrumCollection(filename){
#ifdef _MSC_VER
  pwiz::msdata::DefaultReaderList readerList;
  //readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_ABI));
  //readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_ABI_T2D));
  readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Agilent));
  readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Bruker));
  //readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Shimadzu));
  readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Thermo));
  readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Waters));
  reader_ = new pwiz::msdata::MSDataFile(filename_, &readerList);
#else
  reader_ = new pwiz::msdata::MSDataFile(filename_);
#endif
  if( reader_ == NULL ){
    carp(CARP_FATAL, "PWIZSpectrumCollection unable to open '%s'.", 
         filename_.c_str());
  }
}

PWIZSpectrumCollection::~PWIZSpectrumCollection() {
  delete reader_;
}

/**
 * Parses the first/last scan from the title
 * \returns whether this was successful.
 * For MGF files that place their scan numbers in the title string
 */
bool PWIZSpectrumCollection::parseFirstLastScanFromTitle(
  string& scan_title_str,
  int& first_scan,
  int& last_scan
  ) {

  first_scan = -1;
  last_scan = -1;
  vector<string> scan_title_tokens;
  tokenize(scan_title_str, scan_title_tokens, '.');
  bool success = false;
  //make sure we have enough tokens and that the last token is dta.
  if ((scan_title_tokens.size() >= 4) && (scan_title_tokens.back().find("dta") == 0)) {
    carp(CARP_DETAILED_DEBUG, "Attempting to parse title:%s", scan_title_str.c_str());
    size_t n = scan_title_tokens.size();

    int title_charge;
    int title_first_scan;
    int title_last_scan;
    //try to parse the first scan, last scan, and charge from the title, keeping track
    //of whether we were successful.

    success = from_string(title_charge, scan_title_tokens[n-2]);
    success &= from_string(title_last_scan, scan_title_tokens[n-3]);
    success &= from_string(title_first_scan, scan_title_tokens[n-4]);

    if (success) {
      //okay we parsed the three numbers, fill in the results.
      carp(CARP_DETAILED_DEBUG, "Title first scan:%i", title_first_scan);
      carp(CARP_DETAILED_DEBUG, "Title last scan:%i" ,title_last_scan);
      carp(CARP_DETAILED_DEBUG, "Title charge:%i", title_charge);
      first_scan = title_first_scan;
      last_scan = title_last_scan;
    }
  }
  return success;
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

  carp(CARP_DEBUG, "Using proteowizard to parse spectra.");

  // todo see getFilters() below
  //  pwiz::analysis::SpectrumListFactory::wrap(*reader, spectrum_filters);

  // TODO add first/last scan to base class
  // get a list of scans to include if requested
  const char* range_string = get_string_parameter_pointer("scan-number");
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
  pwiz::msdata::CVID native_id_format = 
    pwiz::msdata::id::getDefaultNativeIDFormat(*reader_);

  // look at all spectra in file
  pwiz::msdata::SpectrumListPtr all_spectra = reader_->run.spectrumListPtr;
  bool get_peaks = true;
  
  int num_spec = all_spectra->size();
  carp(CARP_DEBUG, "PWIZ:Number of spectra:%i",num_spec);
  bool assign_new_scans = false;
  int scan_counter = 0;
  for(int spec_idx = 0; spec_idx < num_spec; spec_idx++){
    carp(CARP_DETAILED_DEBUG, "Parsing spectrum index %d.", spec_idx);
    pwiz::msdata::SpectrumPtr spectrum;
    
    try {
      spectrum = all_spectra->spectrum(spec_idx, get_peaks);
    } catch (boost::bad_lexical_cast) {
      carp(CARP_FATAL, "boost::bad_lexical_cast occured while parsing spectrum.\nDoes your spectra contain z-lines?");
    }
    // skip if not ms2
    if( spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>() != 2 ){
      continue;
    }

    // check that scan number is in range
    int scan_number_begin, scan_number_end;
    if (!assign_new_scans) {
      string ms_peak_list_scans = spectrum->cvParam(pwiz::msdata::MS_peak_list_scans).value;
      string ms_spectrum_title = spectrum->cvParam(pwiz::msdata::MS_spectrum_title).value;
      carp(CARP_DEBUG, "ms_peak_list_scans:%s", ms_peak_list_scans.c_str());
      carp(CARP_DEBUG, "ms_spectrum_title:%s", ms_spectrum_title.c_str());
      if (ms_peak_list_scans.empty() || !get_first_last_scan_from_string(ms_peak_list_scans, scan_number_begin, scan_number_end)) {
        if (ms_spectrum_title.empty() || !parseFirstLastScanFromTitle(ms_spectrum_title, scan_number_begin, scan_number_end)) {
          string scan_value = pwiz::msdata::id::translateNativeIDToScanNumber(
          native_id_format, spectrum->id);
          carp(CARP_DEBUG, "scan_value:%s", scan_value.c_str());
          if (scan_value.empty() || !get_range_from_string<int>(
            scan_value.c_str(), scan_number_begin, scan_number_end)) {
              assign_new_scans = true;
              carp(CARP_ERROR, "Pwiz parser could not determine scan numbers "
                         "for this file, assigning new scan numbers.");
          } else {
            carp(CARP_DEBUG, "found scan:%i-%i from native id", scan_number_begin, scan_number_end);
          }
        } else {
          carp(CARP_DEBUG, "found scan:%i-%i from ms_spectrum_title", scan_number_begin, scan_number_end);
        }
      } else {
        carp(CARP_DEBUG, "found scan:%i-%i from ms_peak_list_scans", scan_number_begin, scan_number_end);
      }
    }
    if (assign_new_scans) {
      scan_number_begin = ++scan_counter;
      scan_number_end = scan_number_begin;
    }
    carp(CARP_DEBUG, "found scan:%i %i-%i", scan_number_begin, first_scan, last_scan);
    if( scan_number_end < first_scan ){
      continue;
    } else if( scan_number_begin > last_scan ){
      break;
    }

    Crux::Spectrum* crux_spectrum = new Crux::Spectrum();
    crux_spectrum->parsePwizSpecInfo(spectrum, scan_number_begin, scan_number_end);

    this->addSpectrumToEnd(crux_spectrum);
    //?delete spectrum?;
  }

  is_parsed_ = true;

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
  Crux::Spectrum* spectrum   ///< Put the spectrum info here
  )
{
  parse();
  return SpectrumCollection::getSpectrum(first_scan, spectrum);
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.
 * \returns The spectrum data from file or NULL.
 */
Crux::Spectrum* PWIZSpectrumCollection::getSpectrum(
  int first_scan      ///< The first scan of the spectrum to retrieve -in
  )
{

  parse();
  return SpectrumCollection::getSpectrum(first_scan);
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
