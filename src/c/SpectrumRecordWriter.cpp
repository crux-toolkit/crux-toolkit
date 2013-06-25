#include <cmath>
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumInfo.hpp"
#include "spectrum.pb.h"
#include "tide/records.h"

#include "SpectrumRecordWriter.h"

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

  // Go through the spectrum list and write each spectrum
  pwiz::msdata::SpectrumList& sl = *(msd->run.spectrumListPtr);
  for (size_t i = 0; i < sl.size(); ++i) {
    pwiz::msdata::SpectrumPtr s = sl.spectrum(i, true);
    if (s->cvParam(pwiz::cv::MS_ms_level).valueAs<int>() > 1 &&
        !s->precursors.empty() &&
        !s->precursors[0].selectedIons.empty()) {
      // Get scan number
      pb::Spectrum pb_spectrum;
      string scan_num = pwiz::msdata::id::translateNativeIDToScanNumber(
        pwiz::cv::MS_scan_number_only_nativeID_format, s->id);
      if (scan_num.empty()) {
        scan_num = "0";
      }
      pb_spectrum.set_spectrum_number(atoi(scan_num.c_str()));
      // Get precursor m/z
      pwiz::msdata::Precursor& precur = s->precursors[0];
      pwiz::msdata::SelectedIon& si = precur.selectedIons[0];
      double mz = si.cvParam(pwiz::cv::MS_selected_ion_m_z).valueAs<double>();
      pb_spectrum.set_precursor_m_z(mz);
      // Get charge states
      pwiz::msdata::CVParam chargeParam = si.cvParam(pwiz::cv::MS_charge_state);
      if (!chargeParam.empty()) {
        pb_spectrum.mutable_charge_state()->Add(chargeParam.valueAs<int>());
      } else {
        for (vector<pwiz::msdata::CVParam>::const_iterator param = si.cvParams.begin();
             param != si.cvParams.end();
             ++param) {
          if (param->cvid == pwiz::cv::MS_possible_charge_state) {
            pb_spectrum.mutable_charge_state()->Add(param->valueAs<int>());
          }
        }
      }
      // Get each m/z, intensity pair
      const pwiz::msdata::BinaryDataArray& mzs = *(s->getMZArray());
      const pwiz::msdata::BinaryDataArray& intensities = *(s->getIntensityArray());
      int mz_denom = getDenom(mzs.data);
      int intensity_denom = getDenom(intensities.data);
      pb_spectrum.set_peak_m_z_denominator(mz_denom);
      pb_spectrum.set_peak_intensity_denominator(intensity_denom);
      google::protobuf::uint64 last = 0;
      int last_index = -1;
      google::protobuf::uint64 intensity_sum = 0;
      for (size_t i = 0; i < s->defaultArrayLength; ++i) {
        google::protobuf::uint64 mz =
          google::protobuf::uint64(mzs.data[i] * mz_denom + 0.5);
        google::protobuf::uint64 intensity =
          google::protobuf::uint64(intensities.data[i] * intensity_denom + 0.5);
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
      // Write spectrum
      writer.Write(&pb_spectrum);
    }
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
