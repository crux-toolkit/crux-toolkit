#include "MSToolkitSpectrumConverter.h"
#include "MSReader.h"
#include "MSToolkitSpectrumCollection.h"
#include "model/Spectrum.h"
#include "parameter.h"
#include "util/Params.h"
#include "util/crux-utils.h"
#include "MSToolkitSpectrumConverter.h"
#include "SpectrumRecordWriter.h"
using namespace std;
MSToolkitSpectrumConverter::MSToolkitSpectrumConverter(const std::string &filename) : SpectrumConverter(filename) {}
bool MSToolkitSpectrumConverter::parse(vector<pb::Spectrum> &all_spectra, int ms_level,
           bool dia_mode) {
  // get a list of scans to include if requested
  string range_string = Params::GetString("scan-number");
  int first_scan;
  int last_scan;

  bool success = get_range_from_string(range_string, first_scan, last_scan);

  if (!success) {
    carp(CARP_FATAL,
         "The scan number range '%s' is invalid. "
         "Must be of the form <first>-<last>.",
         range_string.c_str());
  }

  carp(CARP_DEBUG, "Using mstoolkit to parse spectra.");

  MSToolkit::MSReader *mst_reader = new MSToolkit::MSReader();
  MSToolkit::Spectrum *mst_spectrum = new MSToolkit::Spectrum();

  // only read ms2 scans
  mst_reader->setFilter(MSToolkit::MS2);
  // read first spectrum
  std::string filename = SpectrumConverter::GetFilename();
  bool ret = mst_reader->readFile(filename.c_str(), *mst_spectrum);
  if (!ret) {
    carp(CARP_FATAL, "MSToolkit: Error reading spectra file: %s",
         filename.c_str());
  }

  while (mst_spectrum->getScanNumber() != 0) {
    // is this a scan to include? if not skip it
    if (mst_spectrum->getScanNumber() < first_scan) {
      mst_reader->readFile(NULL, *mst_spectrum);
      continue;
    }
    // are we past the last scan?
    if (mst_spectrum->getScanNumber() > last_scan) {
      break;
    }
    Crux::Spectrum *parsed_spectrum = new Crux::Spectrum();
    if (parsed_spectrum->parseMstoolkitSpectrum(mst_spectrum,
                                                filename.c_str())) {
      parsed_spectrum->sortPeaks(_PEAK_LOCATION);
      vector<pb::Spectrum> pb_spectra = SpectrumRecordWriter::getPbSpectra(parsed_spectrum);
      for (vector<pb::Spectrum>::const_iterator j = pb_spectra.begin(); j != pb_spectra.end(); ++j) {
        assert(j->has_neutral_mass());
        all_spectra.push_back(*j);
      }
    } else {
      delete parsed_spectrum;
    }

    mst_reader->readFile(NULL, *mst_spectrum);
  }
  delete mst_spectrum;
  delete mst_reader;

  return true;
}