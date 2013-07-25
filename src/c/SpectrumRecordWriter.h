#ifndef SPECTRUM_RECORD_WRITER_H
#define SPECTRUM_RECORD_WRITER_H

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

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
