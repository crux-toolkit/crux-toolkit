#ifndef SPECTRUM_RECORD_WRITER_H
#define SPECTRUM_RECORD_WRITER_H

#include "Spectrum.h"
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
   * Return a pb::Spectrum from a Crux::Spectrum
   * Returns a default instance if there is a problem
   */
  static std::vector<pb::Spectrum> getPbSpectra(
    const Crux::Spectrum* s
  );

  /**
   * Add peaks to a pb::Spectrum
   */
  static void addPeaks(
    pb::Spectrum* spectrum,
    const Crux::Spectrum* s
  );

  /**
   * Check if this peak should be excluded
   * Precursor m/z and charges must be set for pb_spectrum
   */
  static bool removePrecursorPeak(
    const pb::Spectrum& pb_spectrum,
    double peakMz
  );

  /**
   * See how much precision is given in the vals array.
   */
  static void getDenoms(
    const Crux::Spectrum* s,  ///< values to check
    int* mzDenom, ///< out parameter for m/z denom
    int* intensityDenom ///< out parameter for intensity denom
  );

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
