#include <cmath>
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "tide/records.h"
#include "tide/mass_constants.h"

#include "Peak.h"
#include "SpectrumRecordWriter.h"
#include "carp.h"
#include "crux-utils.h"

// For printing uint64_t values
#define __STDC_FORMAT_MACROS
#ifndef _MSC_VER
#include <inttypes.h>
#endif

int SpectrumRecordWriter::scanCounter_ = 0;
int SpectrumRecordWriter::removePrecursorPeak_ = 0;
FLOAT_T SpectrumRecordWriter::removePrecursorTolerance_ = 0;

/**
 * Converts a spectra file to spectrumrecords format for use with tide-search.
 * Spectra file is read by pwiz. Returns true on successful conversion.
 */
bool SpectrumRecordWriter::convert(
  const string& infile, ///< spectra file to convert
  string outfile  ///< spectrumrecords file to output
) {
  // Check options to remove peaks around precursor
  removePrecursorPeak_ = get_boolean_parameter("remove-precursor-peak");
  removePrecursorTolerance_ = get_double_parameter("remove-precursor-tolerance");
  /*if (removePrecursorPeak_ < 0 || removePrecursorPeak_ > 2) {
    carp(CARP_FATAL, "remove_precursor_peak must be 0, 1, or 2.");
  }*/

  // Open infile
  pwiz::msdata::MSDataFile* msd;
  try {
    msd = new pwiz::msdata::MSDataFile(infile);
  } catch (...) {
    return false;
  }

  // Write outfile
  pb::Header header;
  header.set_file_type(pb::Header::SPECTRA);

  const vector<pwiz::msdata::SourceFilePtr>& sourceFilePtrs =
    msd->fileDescription.sourceFilePtrs;
  for (vector<pwiz::msdata::SourceFilePtr>::const_iterator s = sourceFilePtrs.begin();
       s != sourceFilePtrs.end();
       ++s) {
    pb::Header_Source* source = header.add_source();
    source->set_filename((*s)->location + "/" + (*s)->name);
    size_t pos = (*s)->name.rfind('.');
    string extension = (pos == string::npos) ? "UNKNOWN" : (*s)->name.substr(pos + 1);
    source->set_filetype(extension);
  }

  header.mutable_spectra_header()->set_sorted(false);

  HeadedRecordWriter writer(outfile, header);
  if (!writer.OK()) {
    return false;
  }

  scanCounter_ = 0;

  // Go through the spectrum list and write each spectrum
  pwiz::msdata::SpectrumList& sl = *(msd->run.spectrumListPtr);
  for (size_t i = 0; i < sl.size(); ++i) {
    pwiz::msdata::SpectrumPtr s;
    try {
      s = sl.spectrum(i, true);
    } catch (...) {
      carp(CARP_FATAL, "Could not parse %s. Error occurred on spectrum %d "
                       "of the spectrum list.", infile.c_str(), i);
    }

    pb::Spectrum pb_spectrum = getPbSpectrum(s);
    if (pb_spectrum.charge_state_size() > 0) {
      writer.Write(&pb_spectrum);
    }
  }

  delete msd;
  return true;
}

/**
 * Return a pb::Spectrum from a pwiz SpectrumPtr
 * If spectrum is ms1, or has no precursors/peaks then return empty pb::Spectrum
 */
pb::Spectrum SpectrumRecordWriter::getPbSpectrum(
  const pwiz::msdata::SpectrumPtr& s
) {
  if (s->cvParam(pwiz::cv::MS_ms_level).valueAs<int>() == 1 ||
      s->precursors.empty() || s->precursors[0].selectedIons.empty() ||
      s->defaultArrayLength == 0) {
    return pb::Spectrum();
  }

  pb::Spectrum pb_spectrum;

  setScanNumber(&pb_spectrum, s);
  setPrecursorMz(&pb_spectrum, s);
  addChargeStates(&pb_spectrum, s);
  addPeaks(&pb_spectrum, s);

  return pb_spectrum;
}

/**
 * Set scan number for a pb::Spectrum
 */
void SpectrumRecordWriter::setScanNumber(
  pb::Spectrum* spectrum,
  const pwiz::msdata::SpectrumPtr& s
) {
  string scan_num = pwiz::msdata::id::translateNativeIDToScanNumber(
    pwiz::cv::MS_scan_number_only_nativeID_format, s->id);
  if (scanCounter_ > 0 || scan_num.empty()) {
    carp_once(CARP_INFO, "Parser could not determine scan numbers for this "
                         "file, using ordinal numbers as scan numbers.");
    spectrum->set_spectrum_number(++scanCounter_);
  } else {
    spectrum->set_spectrum_number(atoi(scan_num.c_str()));
  }
}

/**
 * Set precursor m/z for a pb::Spectrum
 */
void SpectrumRecordWriter::setPrecursorMz(
  pb::Spectrum* spectrum,
  const pwiz::msdata::SpectrumPtr& s
) {
  pwiz::msdata::SelectedIon& si = s->precursors[0].selectedIons[0];
  double mz = si.cvParam(pwiz::cv::MS_selected_ion_m_z).valueAs<double>();
  spectrum->set_precursor_m_z(mz);
}

/**
 * Add charge state information to a pb::Spectrum
 */
void SpectrumRecordWriter::addChargeStates(
  pb::Spectrum* spectrum,
  const pwiz::msdata::SpectrumPtr& s
) {
  pwiz::msdata::SelectedIon& si = s->precursors[0].selectedIons[0];
  pwiz::msdata::CVParam chargeParam = si.cvParam(pwiz::cv::MS_charge_state);
  if (!chargeParam.empty()) {
    spectrum->mutable_charge_state()->Add(chargeParam.valueAs<int>());
  } else {
    for (vector<pwiz::msdata::CVParam>::const_iterator param = si.cvParams.begin();
         param != si.cvParams.end();
         ++param) {
      if (param->cvid == pwiz::cv::MS_possible_charge_state) {
        spectrum->mutable_charge_state()->Add(param->valueAs<int>());
      }
    }
    if (spectrum->charge_state_size() == 0) {
      // Need to calculate charge states
      addCalculatedChargeStates(spectrum, s);
    }
  }
}

/**
 * Add calculated charge state information to a pb::Spectrum
 */
void SpectrumRecordWriter::addCalculatedChargeStates(
  pb::Spectrum* spectrum,
  const pwiz::msdata::SpectrumPtr& s
) {
  sortSpectrumPtrPeaks(s);
  const pwiz::msdata::BinaryDataArray& mzs = *(s->getMZArray());
  const pwiz::msdata::BinaryDataArray& intensities = *(s->getIntensityArray());
  vector<Peak*> peaks;
  for (size_t i = 0; i < s->defaultArrayLength; ++i) {
    peaks.push_back(new Peak(intensities.data[i], mzs.data[i]));
  }
  switch (choose_charge(spectrum->precursor_m_z(), peaks)) {
  case SINGLE_CHARGE_STATE:
    spectrum->mutable_charge_state()->Add(1);
    break;
  case MULTIPLE_CHARGE_STATE:
    spectrum->mutable_charge_state()->Add(2);
    spectrum->mutable_charge_state()->Add(3);
    break;
  default:
    carp(CARP_FATAL, "Could not determine charge state for scan %d",
         spectrum->spectrum_number());
  }
  for (vector<Peak*>::iterator i = peaks.begin(); i != peaks.end(); ++i) {
    delete *i;
  }
}

/**
 * Add peaks to a pb::Spectrum
 */
void SpectrumRecordWriter::addPeaks(
  pb::Spectrum* spectrum,
  const pwiz::msdata::SpectrumPtr& s
) {
  // Get each m/z, intensity pair
  const pwiz::msdata::BinaryDataArray& mzs = *(s->getMZArray());
  const pwiz::msdata::BinaryDataArray& intensities = *(s->getIntensityArray());
  int mz_denom = getDenom(mzs.data);
  int intensity_denom = getDenom(intensities.data);
  spectrum->set_peak_m_z_denominator(mz_denom);
  spectrum->set_peak_intensity_denominator(intensity_denom);
  uint64_t last = 0;
  int last_index = -1;
  uint64_t intensity_sum = 0;

  for (size_t i = 0; i < s->defaultArrayLength; ++i) {
    const double& peakMz = mzs.data[i];
    if (removePrecursorPeak(*spectrum, peakMz)) {
      continue;
    }
    uint64_t mz = peakMz * mz_denom + 0.5;
    uint64_t intensity = intensities.data[i] * intensity_denom + 0.5;
    if (mz < last) {
      // Unsorted peaks, sort and restart loop
      sortSpectrumPtrPeaks(s);

      spectrum->clear_peak_m_z();
      spectrum->clear_peak_intensity();
      last = 0;
      last_index = -1;
      intensity_sum = 0;
      i = -1;
      continue;
    } else if (mz == last) {
      intensity_sum += intensity;
      spectrum->set_peak_intensity(last_index, intensity_sum);
    } else {
      spectrum->add_peak_m_z(mz - last);
      spectrum->add_peak_intensity(intensity);
      last = mz;
      intensity_sum = intensity;
      ++last_index;
    }
  }
}

/**
 * Check if this peak should be excluded
 * Charge states must be set for pb_spectrum
 */
bool SpectrumRecordWriter::removePrecursorPeak(
  const pb::Spectrum& pb_spectrum,
  double peakMz
) {
  switch (removePrecursorPeak_) {
  case 1:
    return fabs(pb_spectrum.precursor_m_z() - peakMz) <= removePrecursorTolerance_;
  /*case 2: {
    // all charge reduced precursor peaks
    for (int i = 0; i < pb_spectrum.charge_state_size(); ++i) {
      int charge = pb_spectrum.charge_state(i);
      double mass = pb_spectrum.precursor_m_z() * charge -
        (charge - 1) - MassConstants::proton;
      for (int j = 1; j < charge; ++j) {
        double mz = (mass + (j - 1) * MassConstants::proton) / j;
        if (fabs(mz - peakMz) <= removePrecursorTolerance_) {
          return true;
        }
      }
    }
    return false;
  }*/
  default:
    return false;
  }
}

/**
 * Sort peaks in a pwiz SpectrumPtr
 */
void SpectrumRecordWriter::sortSpectrumPtrPeaks(
  const pwiz::msdata::SpectrumPtr& s
) {
  vector<pwiz::msdata::MZIntensityPair> sortedPeaks;
  s->getMZIntensityPairs(sortedPeaks);
  sort(sortedPeaks.begin(), sortedPeaks.end(), comparePeaks);
  s->setMZIntensityPairs(sortedPeaks, pwiz::cv::MS_number_of_detector_counts);
}

/**
 * See whether all vals can be accomodated by denom when rendered as a fraction.
 */
bool SpectrumRecordWriter::checkDenom(
  const vector<double>& vals, ///< values to check
  int denom ///< denominator to check against
) {
  double d_denom = denom;
  for (int i = 0; i < vals.size(); ++i) {
    double x = vals[i] * d_denom;
    if (fabs(x - google::protobuf::uint64(x+0.5)) >= 0.001)
      return false;
  }
  return true;
}

/**
 * See how much precision is given in the vals array.
 */
int SpectrumRecordWriter::getDenom(
  const vector<double>& vals  ///< values to check
) {
  const int kMaxPrecision = 10000; // store at most 4 digits of precision
  for (int precision = 1; precision < kMaxPrecision; precision *= 10)
    if (checkDenom(vals, precision))
      return precision;
  return kMaxPrecision;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
