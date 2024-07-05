/**
 * \file PWIZSpectrumCollection.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 21 April 2011
 * \brief Class to read spectra files using the proteowizard library.
 */
#include "PWIZSpectrumCollection.h" 
#include "util/crux-utils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "parameter.h"
#include <iostream>
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#if defined (_MSC_VER) &&  defined(INCLUDE_VENDOR_LIBRARIES)
#include "pwiz/data/msdata/DefaultReaderList.hpp"
//#include "pwiz/data/vendor_readers/ABI/Reader_ABI.hpp"
//#include "pwiz/data/vendor_readers/ABI/T2D/Reader_ABI_T2D.hpp"
#include "pwiz/data/vendor_readers/Agilent/Reader_Agilent.hpp"
#include "pwiz/data/vendor_readers/Bruker/Reader_Bruker.hpp"
#include "pwiz/data/vendor_readers/Shimadzu/Reader_Shimadzu.hpp"
#include "pwiz/data/vendor_readers/Thermo/Reader_Thermo.hpp"
//#include "pwiz/data/vendor_readers/Waters/Reader_Waters.hpp"
#endif


using namespace std;
/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does not parse file. 
 */
PWIZSpectrumCollection::PWIZSpectrumCollection(
  const string& filename   ///< The spectrum collection filename.
) : SpectrumCollection(filename) {
#if defined(_MSC_VER) && defined(INCLUDE_VENDOR_LIBRARIES)
  pwiz::msdata::DefaultReaderList readerList;
  //readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_ABI));
  //readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_ABI_T2D));
  readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Agilent));
  //readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Bruker));
  readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Shimadzu));
  readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Thermo));
  //readerList.push_back(pwiz::msdata::ReaderPtr(new pwiz::msdata::Reader_Waters));
  carp(CARP_DETAILED_INFO, "Support for vendor specific formats enabled.");  
  try {
     reader_ = new pwiz::msdata::MSDataFile(filename_, &readerList);
  } catch (const runtime_error& error) {
    carp(CARP_FATAL, "Unable to parse spectrum file %s. Error: %s.", filename_.c_str(), error.what());
  }
#else
  carp(CARP_DETAILED_INFO, "Support for vendor specific formats not enabled.");  
  try {
    reader_ = new pwiz::msdata::MSDataFile(filename_);
  }
  catch (const runtime_error& error) {
    carp(CARP_FATAL, "Unable to parse spectrum file %s. Error: %s.", filename_.c_str(), error.what());  
  }
#endif
  if( reader_ == NULL ) {
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
 * Assumes format style is either TITLE=ScaffoldIDNumber_853_o100720_3prot_01.10257.10257.3.dta
 * or
 * TITLE=Yp_D27_P0_D3.271.271.1 File:"Yp_D27_P0_D3.raw", NativeID:"controllerType=0 controllerNumber=1 scan=271"
 * Note that input string does not contain "TITLE="
 */
bool PWIZSpectrumCollection::parseFirstLastScanFromTitle(
  string& scan_title_str,
  int& first_scan,
  int& last_scan
  ) {
  int title_charge;
  int title_first_scan;
  int title_last_scan;

  first_scan = -1;
  last_scan = -1;
  vector<string> scan_title_tokens = StringUtils::Split(scan_title_str, '.');
  bool success = false;
  // make sure we have enough tokens and that the last token is dta.
  // FOr TITLE=ScaffoldIDNumber_853_o100720_3prot_01.10257.10257.3.dta format
  if ((scan_title_tokens.size() >= 4) && (scan_title_tokens.back().find("dta") == 0)) {
    carp(CARP_DETAILED_DEBUG, "Attempting to parse title:%s", scan_title_str.c_str());
    size_t n = scan_title_tokens.size();

    // try to parse the first scan, last scan, and charge from the title, keeping track
    // of whether we were successful.
    success = StringUtils::TryFromString(scan_title_tokens[n-2], &title_charge);
    success &= StringUtils::TryFromString(scan_title_tokens[n-3], &title_last_scan);
    success &= StringUtils::TryFromString(scan_title_tokens[n-4], &title_first_scan);

    if (success) {
      // okay we parsed the three numbers, fill in the results.
      carp(CARP_DETAILED_DEBUG, "Title first scan:%i", title_first_scan);
      carp(CARP_DETAILED_DEBUG, "Title last scan:%i" , title_last_scan);
      carp(CARP_DETAILED_DEBUG, "Title charge:%i", title_charge);
      first_scan = title_first_scan;
      last_scan = title_last_scan;
    }
  }

  // For TITLE=Yp_D27_P0_D3.271.271.1 File:"Yp_D27_P0_D3.raw", NativeID:"controllerType=0 controllerNumber=1 scan=271" format
  if (success == false) {
    scan_title_tokens = StringUtils::Split(scan_title_str, " File:");
	string first_token = scan_title_tokens[0];
    vector<string> first_token_tokens = StringUtils::Split(first_token, ".");

    if (first_token_tokens.size() >= 4) {
      carp(CARP_DETAILED_DEBUG, "Attempting to parse title:%s", scan_title_str.c_str());
      size_t n = first_token_tokens.size();

      success = StringUtils::TryFromString(first_token_tokens[n-1], &title_charge);
      success &= StringUtils::TryFromString(first_token_tokens[n-2], &title_last_scan);
      success &= StringUtils::TryFromString(first_token_tokens[n-3], &title_first_scan);

      if (success) {
        carp(CARP_DETAILED_DEBUG, "Title first scan:%i", title_first_scan);
        carp(CARP_DETAILED_DEBUG, "Title last scan:%i" , title_last_scan);
        carp(CARP_DETAILED_DEBUG, "Title charge:%i", title_charge);
        first_scan = title_first_scan;
        last_scan = title_last_scan;
      }
    }
  }
  return success;
}

/**
 * Parses all the spectra from file designated by the filename member
 * variable.
 * \returns True if the spectra are parsed successfully. False if otherwise.
 */
bool PWIZSpectrumCollection::parse(int ms_level, bool dia_mode) {
  // spectrum_collection has already been parsed
  if(is_parsed_) {
    return false;
  }
  carp(CARP_DEBUG, "ms level: %d", ms_level );

  carp(CARP_DEBUG, "Using proteowizard to parse spectra.");

  // todo see getFilters() below
  //  pwiz::analysis::SpectrumListFactory::wrap(*reader, spectrum_filters);

  // TODO add first/last scan to base class
  // get a list of scans to include if requested
  string range_string = Params::GetString("scan-number");
  int first_scan = -1;
  int last_scan = -1;

  // added by Yang
  int curr_ms1_scan = -1;


  get_range_from_string(range_string, first_scan, last_scan);

  if (first_scan == -1 || last_scan == -1) {
    carp(CARP_FATAL, "The scan number range '%s' is invalid. "
         "Must be of the form <first>-<last>.", range_string.c_str());
  }

  // get info for translating identifiers into scan numbers
  pwiz::msdata::CVID native_id_format =
    pwiz::msdata::id::getDefaultNativeIDFormat(*reader_);

  // look at all spectra in file
  pwiz::msdata::SpectrumListPtr all_spectra = reader_->run.spectrumListPtr;
  
  int num_spec = all_spectra->size();
  carp(CARP_DEBUG, "PWIZ:Number of spectra:%i", num_spec);
  bool assign_new_scans = false;
  int scan_counter = 0;
  for (int spec_idx = 0; spec_idx < num_spec; spec_idx++) {
    carp(CARP_DETAILED_DEBUG, "Parsing spectrum index %d.", spec_idx);
    pwiz::msdata::SpectrumPtr spectrum;
    try {
      spectrum = all_spectra->spectrum(spec_idx, true);
    } catch (boost::bad_lexical_cast) {
      carp(CARP_FATAL, "boost::bad_lexical_cast occured while parsing spectrum.\n"
                       "Do your spectra contain z-lines?");
    }

    int curr_ms_level = spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>();
    carp(CARP_DEBUG, "ms_level=%d\t", curr_ms_level );

    // skip if no peaks or not ms2
    // if (spectrum->defaultArrayLength < 1 || spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>() != 2) { continue; }
    // if (spectrum->defaultArrayLength < 1 || spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>() != 2) { continue; }


    // check that scan number is in range
    int scan_number_begin, scan_number_end;
    if (!assign_new_scans) {
      string ms_peak_list_scans = spectrum->cvParam(pwiz::msdata::MS_peak_list_scans).value;
      string ms_spectrum_title = spectrum->cvParam(pwiz::msdata::MS_spectrum_title).value;
      carp(CARP_DETAILED_DEBUG, "ms_peak_list_scans:%s", ms_peak_list_scans.c_str());
      carp(CARP_DETAILED_DEBUG, "ms_spectrum_title:%s", ms_spectrum_title.c_str());
      if (ms_peak_list_scans.empty() || !get_first_last_scan_from_string(ms_peak_list_scans, scan_number_begin, scan_number_end)) {
        if (ms_spectrum_title.empty() || !parseFirstLastScanFromTitle(ms_spectrum_title, scan_number_begin, scan_number_end)) {
          string scan_value = pwiz::msdata::id::translateNativeIDToScanNumber(native_id_format, spectrum->id);
          carp(CARP_DETAILED_DEBUG, "scan_value:%s", scan_value.c_str());
          if (scan_value.empty() || !get_range_from_string<int>(
            scan_value.c_str(), scan_number_begin, scan_number_end)) {
              assign_new_scans = true;
              carp(CARP_ERROR, "Proteowizard parser could not determine scan numbers "
                         "for this file, assigning new scan numbers.");
          } else {
            carp(CARP_DETAILED_DEBUG, "found scan:%i-%i from native id", scan_number_begin, scan_number_end);
          }
        } else {
          carp(CARP_DETAILED_DEBUG, "found scan:%i-%i from ms_spectrum_title", scan_number_begin, scan_number_end);
        }
      } else {
        carp(CARP_DETAILED_DEBUG, "found scan:%i-%i from ms_peak_list_scans", scan_number_begin, scan_number_end);
      }
      if (scan_number_begin == 0) {
        // PWiz assigns scan numbers starting from 0 if they are missing. In this case, we re-assign starting from 1 below.
        carp_once(CARP_INFO, "Parser could not determine scan numbers for this "
                             "file, using ordinal numbers as scan numbers.");
        assign_new_scans = true;
      }
    }
    if (assign_new_scans) {
      carp_once(CARP_WARNING,
           "Proteowizard parser could not determine scan numbers "
           "for this file. Assigning new scan numbers.");
      scan_number_begin = ++scan_counter;
      scan_number_end = scan_number_begin;
    }
    carp(CARP_DETAILED_DEBUG, "found scan:%i %i %i-%i", scan_number_begin, scan_number_end, first_scan, last_scan);
    if( scan_number_end < first_scan ) {
      continue;
    } else if( scan_number_begin > last_scan ) {
      break;
    }

    // added by Yang
    if (spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>() == 1) {
    	curr_ms1_scan = scan_number_begin;
    	if (scan_number_begin != scan_number_end) { carp(CARP_FATAL, "scan_number_begin %d should equal to scan_number_end %d.", scan_number_begin, scan_number_end); }
    }
    // skip if no peaks or ms_level doesn't match
    if (spectrum->defaultArrayLength < 1 || spectrum->cvParam(pwiz::msdata::MS_ms_level).valueAs<int>() != ms_level) { continue; }


    Crux::Spectrum* crux_spectrum = new Crux::Spectrum();
    if (crux_spectrum->parsePwizSpecInfo(spectrum, scan_number_begin, scan_number_end, dia_mode)) {
    	// added by Yang
    	crux_spectrum->setMS1Scan(curr_ms1_scan);
    	carp(CARP_DETAILED_DEBUG, "curr_ms1_scan: %d ", curr_ms1_scan );

    	addSpectrumToEnd(crux_spectrum);
    	spectraByScan_[scan_number_begin] = crux_spectrum;
    } else {
    	delete crux_spectrum;
    }
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
  ) {
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
  ) {
  parse();
  return SpectrumCollection::getSpectrum(first_scan);
}

/*
void getFilters() {  
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
