#include "SpectrumRecordSpectrumCollection.h"
#include "app/tide/records.h"
#include "app/tide/spectrum_collection.h"
#include "util/FileUtils.h"

using namespace std;

SpectrumRecordSpectrumCollection::SpectrumRecordSpectrumCollection(
  const string& filename
): SpectrumCollection(filename) {
}

SpectrumRecordSpectrumCollection::~SpectrumRecordSpectrumCollection() {
}

bool SpectrumRecordSpectrumCollection::IsSpectrumRecordFile(const string& filename) {
  if (!FileUtils::Exists(filename) || FileUtils::IsDir(filename)) {
    return false;
  }
  pb::Header header;
  HeadedRecordReader reader(filename, &header);
  return header.file_type() == pb::Header::SPECTRA;
}

bool SpectrumRecordSpectrumCollection::parse() {
  if (is_parsed_) {
    return false;
  }

  is_parsed_ = true;

  pb::Header header;
  HeadedRecordReader reader(filename_, &header);
  if (header.file_type() != pb::Header::SPECTRA) {
    carp(CARP_ERROR, "File '%s' is not a valid spectrum records file", filename_.c_str());
    return false;
  }
  pb::Spectrum pb_spectrum;
  while (!reader.Done()) {
    reader.Read(&pb_spectrum);
    vector<int> charges;
    for (int i = 0; i < pb_spectrum.charge_state_size(); i++) {
      charges.push_back(pb_spectrum.charge_state(i));
    }
    Crux::Spectrum* spectrum = new Crux::Spectrum(
      pb_spectrum.spectrum_number(),
      pb_spectrum.spectrum_number(),
      pb_spectrum.precursor_m_z(),
      charges,
      filename_);
    double mzDenom = pb_spectrum.peak_m_z_denominator();
    double intensityDenom = pb_spectrum.peak_intensity_denominator();
    uint64_t total = 0;
    for (int i = 0; i < pb_spectrum.peak_m_z_size(); i++) {
      total += pb_spectrum.peak_m_z(i);
      spectrum->addPeak(
        pb_spectrum.peak_intensity(i) / intensityDenom,
        total / mzDenom);
    }
    addSpectrum(spectrum);
  }
  if (!reader.OK()) {
    carp(CARP_ERROR, "Error reading spectrum records file '%s'", filename_.c_str());
    for (deque<Crux::Spectrum*>::const_iterator i = spectra_.begin(); i != spectra_.end(); i++) {
      delete *i;
    }
    spectra_.clear();
    return false;
  }
  return true;
}

Crux::Spectrum* SpectrumRecordSpectrumCollection::getSpectrum(int first_scan) {
  parse();
  return SpectrumCollection::getSpectrum(first_scan);
}

bool SpectrumRecordSpectrumCollection::getSpectrum(int first_scan, Crux::Spectrum* spectrum) {
  parse();
  return SpectrumCollection::getSpectrum(first_scan, spectrum);
}

