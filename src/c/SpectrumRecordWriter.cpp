#include <cmath>
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "spectrum.pb.h"
#include "tide/records.h"

#include "Peak.h"
#include "SpectrumRecordWriter.h"
#include "carp.h"
#include "crux-utils.h"

// For printing uint64_t values
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

/**
 * Converts a spectra file to spectrumrecords format for use with tide-search.
 * Spectra file is read by pwiz. Returns true on successful conversion.
 */
bool SpectrumRecordWriter::convert(
  const string& infile, ///< spectra file to convert
  string outfile  ///< spectrumrecords file to output
) {

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

  vector<Peak*> peaks;
  int scan_counter = 0;

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
    if (s->cvParam(pwiz::cv::MS_ms_level).valueAs<int>() == 1 ||
        s->precursors.empty() || s->precursors[0].selectedIons.empty() ||
        s->defaultArrayLength == 0) {
      continue;
    }
    // Get scan number
    pb::Spectrum pb_spectrum;
    string scan_num = pwiz::msdata::id::translateNativeIDToScanNumber(
      pwiz::cv::MS_scan_number_only_nativeID_format, s->id);
    if (scan_counter > 0 || scan_num.empty()) {
      carp_once(CARP_INFO, "Parser could not determine scan numbers for this "
                           "file, using ordinal numbers as scan numbers.");
      pb_spectrum.set_spectrum_number(++scan_counter);
    } else {
      pb_spectrum.set_spectrum_number(atoi(scan_num.c_str()));
    }
    // Get precursor m/z
    pwiz::msdata::Precursor& precur = s->precursors[0];
    pwiz::msdata::SelectedIon& si = precur.selectedIons[0];
    double mz = si.cvParam(pwiz::cv::MS_selected_ion_m_z).valueAs<double>();
    pb_spectrum.set_precursor_m_z(mz);
    // Get charge states
    pwiz::msdata::CVParam chargeParam = si.cvParam(pwiz::cv::MS_charge_state);
    bool calcCharge = false;
    if (!chargeParam.empty()) {
      pb_spectrum.mutable_charge_state()->Add(chargeParam.valueAs<int>());
    } else {
      bool possibleCharge = false;
      for (vector<pwiz::msdata::CVParam>::const_iterator param = si.cvParams.begin();
           param != si.cvParams.end();
           ++param) {
        if (param->cvid == pwiz::cv::MS_possible_charge_state) {
          possibleCharge = true;
          pb_spectrum.mutable_charge_state()->Add(param->valueAs<int>());
        }
      }
      if (!possibleCharge) {
        carp(CARP_WARNING, "Scan %s has no charge state, it will be assigned",
             scan_num.c_str());
        calcCharge = true;
      }
    }
    // Get each m/z, intensity pair
    const pwiz::msdata::BinaryDataArray& mzs = *(s->getMZArray());
    const pwiz::msdata::BinaryDataArray& intensities = *(s->getIntensityArray());
    int mz_denom = getDenom(mzs.data);
    int intensity_denom = getDenom(intensities.data);
    pb_spectrum.set_peak_m_z_denominator(mz_denom);
    pb_spectrum.set_peak_intensity_denominator(intensity_denom);
    uint64_t last = 0;
    int last_index = -1;
    uint64_t intensity_sum = 0;
    for (size_t j = 0; j < s->defaultArrayLength; ++j) {
      if (calcCharge) {
        peaks.push_back(new Peak(intensities.data[j], mzs.data[j]));
      }
      uint64_t mz = mzs.data[j] * mz_denom + 0.5;
      uint64_t intensity = intensities.data[j] * intensity_denom + 0.5;
      if (mz < last) {
        // Unsorted peaks, sort and restart loop
        carp_once(CARP_WARNING, "Spectrum with unsorted peaks found. "
                                "All unsorted peaks will be sorted.");
        vector<pwiz::msdata::MZIntensityPair> sortedPeaks;
        for (size_t k = 0; k < s->defaultArrayLength; ++k) {
          sortedPeaks.push_back(pwiz::msdata::MZIntensityPair(
            mzs.data[k], intensities.data[k]));
        }
        sort(sortedPeaks.begin(), sortedPeaks.end(), comparePeaks);
        s->setMZIntensityPairs(sortedPeaks, pwiz::cv::MS_number_of_detector_counts);

        peaks.clear();
        pb_spectrum.clear_peak_m_z();
        pb_spectrum.clear_peak_intensity();
        last = 0;
        last_index = -1;
        intensity_sum = 0;
        j = -1;
        continue;
      }
      if (mz == last) {
        intensity_sum += intensity;
        pb_spectrum.set_peak_intensity(last_index, intensity_sum);
      } else {
        pb_spectrum.add_peak_m_z(mz - last);
        pb_spectrum.add_peak_intensity(intensity);
        last = mz;
        intensity_sum = intensity;
        ++last_index;
      }
    }
    if (!peaks.empty()) {
      switch (choose_charge(mz, peaks)) {
      case SINGLE_CHARGE_STATE:
        pb_spectrum.mutable_charge_state()->Add(1);
        break;
      case MULTIPLE_CHARGE_STATE:
        pb_spectrum.mutable_charge_state()->Add(2);
        pb_spectrum.mutable_charge_state()->Add(3);
        break;
      default:
        carp(CARP_FATAL, "Could not determine charge state for scan %d",
             scan_num.c_str());
      }
      for (vector<Peak*>::iterator j = peaks.begin(); j != peaks.end(); ++j) {
        delete *j;
      }
      peaks.clear();
    }
    // Write spectrum
    writer.Write(&pb_spectrum);
  }

  delete msd;
  return true;
}

/**
 * See whether all vals can be accomodated by denom when rendered as a fraction.
 */
bool SpectrumRecordWriter::checkDenom(
  const vector<double>& vals, ///< values to check
  int denom ///< denominator to check against
) {
  double d_denom = denom;
  for (int i = 0; i < vals.size(); ++ i) {
    double x = vals[i] * d_denom;
    if (fabs(x - google::protobuf::uint64(x+0.5)) < 0.001)
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
  const int kMaxPrecision = 10000; // store at most 3 digits of precision
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
