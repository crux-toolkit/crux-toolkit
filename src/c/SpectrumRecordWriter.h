#ifndef SPECTRUM_RECORD_WRITER_H
#define SPECTRUM_RECORD_WRITER_H

#include "pwiz/data/msdata/MSData.hpp"

#include "spectrum.pb.h"

using namespace std;

/**
 * A class for converting spectra file to the spectrumrecords format for use
 * with tide-search.
 */
class SpectrumRecordWriter {

public:

  /**
   * Converts a spectra file to spectrumrecords format for use with tide-search.
   * Spectra file is read by pwiz. Returns true on successful conversion.
   */
  static bool convert(
    const string& infile, ///< spectra file to convert
    string outfile  ///< spectrumrecords file to output
  );

protected:

  static int scanCounter_;
  static int removePrecursorPeak_;
  static FLOAT_T removePrecursorTolerance_;

  /**
   * Return a pb::Spectrum from a pwiz SpectrumPtr
   * If spectrum is ms1, or has no precursors/peaks then return empty pb::Spectrum
   */
  static pb::Spectrum getPbSpectrum(
    const pwiz::msdata::SpectrumPtr& s
  );

  /**
   * Set scan number for a pb::Spectrum
   */
  static void setScanNumber(
    pb::Spectrum* spectrum,
    const pwiz::msdata::SpectrumPtr& s
  );

  /**
   * Set precursor m/z for a pb::Spectrum
   */
  static void setPrecursorMz(
    pb::Spectrum* spectrum,
    const pwiz::msdata::SpectrumPtr& s
  );

  /**
   * Add charge state information to a pb::Spectrum
   */
  static void addChargeStates(
    pb::Spectrum* spectrum,
    const pwiz::msdata::SpectrumPtr& s
  );

  /**
   * Add calculated charge state information to a pb::Spectrum
   */
  static void addCalculatedChargeStates(
    pb::Spectrum* spectrum,
    const pwiz::msdata::SpectrumPtr& s
  );

  /**
   * Add peaks to a pb::Spectrum
   */
  static void addPeaks(
    pb::Spectrum* spectrum,
    const pwiz::msdata::SpectrumPtr& s
  );

  /**
   * Check if this peak should be excluded
   * Charge states must be set for pb_spectrum
   */
  static bool removePrecursorPeak(
    const pb::Spectrum& pb_spectrum,
    double peakMz
  );

  /**
   * Sort peaks in a pwiz SpectrumPtr
   */
  static void sortSpectrumPtrPeaks(
    const pwiz::msdata::SpectrumPtr& s
  );

  /**
   * See whether all vals can be accomodated by denom when rendered as a fraction.
   */
  static bool checkDenom(
    const vector<double>& vals, ///< values to check
    int denom ///< denominator to check against
  );

  /**
   * See how much precision is given in the vals array.
   */
  static int getDenom(
    const vector<double>& vals  ///< values to check
  );

  static bool comparePeaks(
    const pwiz::msdata::MZIntensityPair& x,
    const pwiz::msdata::MZIntensityPair& y
  ) {
    return x.mz < y.mz;
  }

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
