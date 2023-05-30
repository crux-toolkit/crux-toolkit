#include <cmath>
#include <memory>
#include "app/tide/records.h"
#include "app/tide/mass_constants.h"

#include "model/Peak.h"
#include "SpectrumCollectionFactory.h"
#include "SpectrumRecordWriter.h"
#include "io/carp.h"
#include "util/crux-utils.h"

// For printing uint64_t values
#define __STDC_FORMAT_MACROS
#ifndef _MSC_VER
#include <inttypes.h>
#endif

int SpectrumRecordWriter::scanCounter_ = 0;

/**
 * Converts a spectra file to spectrumrecords format for use with tide-search.
 * Spectra file is read by pwiz. Returns true on successful conversion.
 */
bool cmp_spectra(pb::Spectrum a1, pb::Spectrum a2)  {
  return a1.precursor_m_z() < a2. precursor_m_z();
}
bool SpectrumRecordWriter::convert(
  const string& infile, ///< spectra file to convert
  string outfile,  ///< spectrumrecords file to output
  int ms_level,   /// MS level to extract (1 or 2)
  bool dia_mode  /// whether it's used in DIAmeter
) {
  carp(CARP_INFO, "Converting ms_level %d ... ", ms_level);
  auto_ptr<Crux::SpectrumCollection> spectra(SpectrumCollectionFactory::create(infile.c_str()));

  // added by Yang
  if ( ms_level < 1 || ms_level > 2 ) { carp(CARP_FATAL, "ms_level must be 1 or 2 instead of %d.", ms_level); }


  // Open infile
  try {
    if (!spectra->parse(ms_level, dia_mode)) {
      return false;
    }
  } catch (const std::exception& e) {
    carp(CARP_ERROR, "%s", e.what());
    return false;
  } catch (...) {
    return false;
  }

  // Write outfile
  pb::Header header;
  header.set_file_type(pb::Header::SPECTRA);

  pb::Header_Source* source = header.add_source();
  source->set_filename(infile);
  size_t pos = infile.rfind('.');
  string extension = (pos == string::npos) ? "UNKNOWN" : infile.substr(pos + 1);
  source->set_filetype(extension);

  header.mutable_spectra_header()->set_sorted(false);

  HeadedRecordWriter writer(outfile, header);
  if (!writer.OK()) {
    return false;
  }

  scanCounter_ = 0;
   carp(CARP_DETAILED_DEBUG, "starting to convert spectrum to pb..." );
  // go through the spectrum list and write each spectrum

  vector<pb::Spectrum> all_spectra;  //Added by AKF

  for (SpectrumIterator i = spectra->begin(); i != spectra->end(); ++i) {
	(*i)->sortPeaks(_PEAK_LOCATION); // Sort by m/z

    vector<pb::Spectrum> pb_spectra = getPbSpectra(*i);
    for (vector<pb::Spectrum>::const_iterator j = pb_spectra.begin();
         j != pb_spectra.end();
         ++j) { 
//      writer.Write(&*j);  //Added by AKF
        assert(j->has_precursor_m_z());
        all_spectra.push_back(*j);
    }
  }
  // Below, added by AKF
  // all_spectra contains all spectra.
  // sort all_spectra
  std::sort(all_spectra.begin(), all_spectra.end(), cmp_spectra);
  for (vector<pb::Spectrum>::const_iterator j = all_spectra.begin();
         j != all_spectra.end();
         ++j) { 
      writer.Write(&*j);  //Added by AKF
  }
  return true;
}

/**
 * Return a pb::Spectrum from a pwiz SpectrumPtr
 * If spectrum is ms1, or has no precursors/peaks then return empty pb::Spectrum
 */
vector<pb::Spectrum> SpectrumRecordWriter::getPbSpectra(
  const Crux::Spectrum* s
) {
  vector<pb::Spectrum> spectra;

  if (s->getNumZStates() == 0 || s->getNumPeaks() == 0) {
	carp(CARP_DETAILED_DEBUG, "numZStates: %d \t numPeaks: %d", s->getNumZStates(), s->getNumPeaks() );
    return spectra;
  }

  // Get scan number
  int scan_num = s->getFirstScan();
  if (scanCounter_ > 0 || scan_num <= 0) {
    carp_once(CARP_INFO, "Parser could not determine scan numbers for this "
                         "file, using ordinal numbers as scan numbers.");
    scan_num = ++scanCounter_;
  }

  const vector<SpectrumZState>& zStates = s->getZStates();
  for (vector<SpectrumZState>::const_iterator i = zStates.begin(); i != zStates.end(); ++i) {
    spectra.push_back(pb::Spectrum());
    pb::Spectrum& newSpectrum = spectra.back();

    // added by Yang
    newSpectrum.set_ms1_spectrum_number(s->getMS1Scan());
    newSpectrum.set_iso_window_lower_mz(s->getIsoWindowLowerMZ());
    newSpectrum.set_iso_window_upper_mz(s->getIsoWindowUpperMZ());

    newSpectrum.set_spectrum_number(scan_num);
    newSpectrum.set_precursor_m_z(i->getMZ());
    newSpectrum.mutable_charge_state()->Add(i->getCharge());
    addPeaks(&newSpectrum, s);
    if (newSpectrum.peak_m_z_size() == 0) {
      spectra.pop_back();
    }
  }

  return spectra;
}

/**
 * Add peaks to a pb::Spectrum
 */
void SpectrumRecordWriter::addPeaks(
  pb::Spectrum* spectrum,
  const Crux::Spectrum* s
) {
  int mz_denom, intensity_denom;
  getDenoms(s, &mz_denom, &intensity_denom);
  spectrum->set_peak_m_z_denominator(mz_denom);
  spectrum->set_peak_intensity_denominator(intensity_denom);
  uint64_t last = 0;
  int last_index = -1;
  uint64_t intensity_sum = 0;

  for (PeakIterator i = s->begin(); i != s->end(); ++i) {
    FLOAT_T peakMz = (*i)->getLocation();
    uint64_t mz = peakMz * mz_denom + 0.5;
    uint64_t intensity = (*i)->getIntensity() * intensity_denom + 0.5;
    if (mz < last) {
      // Unsorted peaks, this should never happen since peaks get sorted earlier
      carp(CARP_FATAL, "Peaks are not sorted");
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
 * See how much precision is given in the peak data
 */
void SpectrumRecordWriter::getDenoms(
  const Crux::Spectrum* s,  ///< spectra with peaks to check
  int* mzDenom, ///< out parameter for m/z denom
  int* intensityDenom ///< out parameter for intensity denom
) {
  const int kMaxPrecision = 10000; // store at most 4 digits of precision
  *mzDenom = kMaxPrecision;
  *intensityDenom = kMaxPrecision;
  for (int precision = 1; precision < kMaxPrecision; precision *= 10) {
    bool mzDenomOk = (*mzDenom == kMaxPrecision);
    bool intensityDenomOk = (*intensityDenom == kMaxPrecision);
    if (!mzDenomOk && !intensityDenomOk) {
      return;
    }
    for (PeakIterator i = s->begin(); i != s->end(); ++i) {
      if (mzDenomOk) {
        double mzX = (*i)->getLocation() * precision;
        if (fabs(mzX - google::protobuf::uint64(mzX + 0.5)) >= 0.001) {
          mzDenomOk = false;
        }
      }
      if (intensityDenomOk) {
        double intensityX = (*i)->getIntensity() * precision;
        if (fabs(intensityX - google::protobuf::uint64(intensityX + 0.5)) >= 0.001) {
          intensityDenomOk = false;
        }
      }
      if (!mzDenomOk && !intensityDenomOk) {
        break;
      }
    }
    if (mzDenomOk) {
      *mzDenom = precision;
    }
    if (intensityDenomOk) {
      *intensityDenom = precision;
    }
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
