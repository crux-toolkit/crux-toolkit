#ifndef XHHC_SCORER_H
#define XHHC_SCORER_H

#include "scorer.h"

#include "xhhc_ion_series.h"

#include <map>
/* All of these are slightly modified versions of scorer.c,
 * made to passed LinkedIonSeries objects rather than crux's 
 * ION_SERIES_T */
#define NORMALIZE 0
//#define MAX_MZ 1200;
//#define MIN_MZ 400;
//#define SCALE_SPECTRUM 0
#define NO_FLANKS 1

class Scorer {

  public:
    Scorer() : print_spectrums_(false), max_mz(0.0) {
      scorer = NULL;
	current_spectrum = NULL;
    }

    Scorer(FLOAT_T a_max_mz);

    FLOAT_T get_max_mz() { return max_mz; }
    
    // main scoring function
    // creates theoretical spectrum from ion series and
    // scores against observed spectrum
    // returns xcorr score
    float score_spectrum_vs_series(Spectrum* spectrum, 
                                   LinkedIonSeries& ion_series);

    static int get_matched_by_ions(Spectrum* spectrum,
				   LinkedIonSeries& ion_series);

    static void add_intensity_map(std::map<int, FLOAT_T>& theoretical, 
                                  int idx, 
                                  FLOAT_T intensity);
/* 
    void add_intensity(
      FLOAT_T* intensity_array, ///< the intensity array to add intensity at index add_idx -out
      int add_idx,            ///< the idex to add the intensity -in
      FLOAT_T intensity         ///< the intensity to add -in
    );
*/
    static bool xlink_create_map_theoretical(
      LinkedIonSeries& ion_series,
      std::map<int, FLOAT_T>& theoretical);

    // modified version of create_intensity_array_theoretical from crux
    bool hhc_create_intensity_array_theoretical(
      LinkedIonSeries& ion_series,
      FLOAT_T* theoretical       ///< the empty theoretical spectrum -out
    );
    void set_print(bool print) { print_spectrums_ = print;}
    FLOAT_T hhc_gen_score_xcorr(
      Spectrum* spectrum,    ///< the spectrum to score -in
      LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
    );

    FLOAT_T getIonCurrentExplained(LinkedIonSeries& ion_series, 
      Spectrum* spectrum, 
      FLOAT_T& explained, 
      int& by_observed);


  private:
    void print_spectrums(FLOAT_T* theoretical, Spectrum* spectrum);
    bool print_spectrums_;
    SCORER_T* scorer;
    Spectrum* current_spectrum;
    FLOAT_T max_mz;
    static int ion_counter;
};

// Everything below here is in a crux .c file but not included in headers, not sure
// how to include them
////////////////////////////////////
/*

#define MAX_I_LINES 2 // number of 'I' lines albe to parse for one spectrum object
#define MAX_D_LINES 2 // number of 'D' lines albe to parse for one spectrum object

struct spectrum{
  int              first_scan;    ///< The number of the first scan
  int              last_scan;     ///< The number of the last scan
  int              id;            ///< A unique identifier
                                  // FIXME, this field is not set when parsing
  SPECTRUM_TYPE_T  spectrum_type; ///< The type of spectrum. 
  FLOAT_T            precursor_mz;  ///< The m/z of precursor (MS-MS spectra)
  int*             possible_z;    ///< The possible charge states of this spectrum
  int              num_possible_z;///< The number of possible charge states of this spectrum
  PEAK_T*          peaks;         ///< The spectrum peaks
  FLOAT_T            min_peak_mz;   ///< The minimum m/z of all peaks
  FLOAT_T            max_peak_mz;   ///< The maximum m/z of all peaks
  int              num_peaks;     ///< The number of peaks
  double           total_energy;  ///< The sum of intensities in all peaks
  char*            filename;      ///< Optional filename
  char*            i_lines[MAX_I_LINES]; ///< store i lines, upto MAX_I_LINES
  char*            d_lines[MAX_D_LINES]; ///< store d lines, upto MAX_D_LINES 
  BOOLEAN_T        has_peaks;  ///< Does the spectrum contain peak information
  BOOLEAN_T        sorted_by_mz; ///< Are the spectrum peaks sorted by m/z...
  BOOLEAN_T        sorted_by_intensity; ///< ... or by intensity?
  BOOLEAN_T        has_mz_peak_array; ///< Is the mz_peak_array populated.
  PEAK_T**         mz_peak_array;  ///< Allows rapid peak retrieval by mz.
};    


// copied from scorer.c in crux
struct scorer {
  SCORER_TYPE_T type; ///< The type of scorer
  FLOAT_T sp_beta; ///< used for Sp: the beta variable 
  FLOAT_T sp_max_mz; ///< used for Sp: the max mz for the intensity array
  int sp_b_y_ion_matched; ///< The most recent ion_collection number of the b, y ion matched while scoring for SP
  int sp_b_y_ion_possible; ///< The most recent ion_collection number of the b, y ion possible while scoring for SP
  FLOAT_T sp_b_y_ion_fraction_matched; ///< The ratio of matched and possible.

  FLOAT_T* intensity_array; ///< used for Sp: the intensity array, which can be indexed using the m/z
  FLOAT_T max_intensity; ///< the max intensity in the intensity array
  BOOLEAN_T initialized; ///< has the scorer been initialized?
  int last_idx; ///< the last index in the array, the data size of the array

  /// used for xcorr
  FLOAT_T* observed; ///< used for Xcorr: observed spectrum intensity array
  FLOAT_T* theoretical; ///< used for Xcorr: theoretical spectrum intensity array
};

FLOAT_T cross_correlation(
  SCORER_T* scorer,  ///< the scorer object that contains observed spectrum -in
  FLOAT_T* theoretical ///< the theoretical spectrum to score against the observed spectrum -in
  );


void add_intensity(
  FLOAT_T* intensity_array, ///< the intensity array to add intensity at index add_idx -out
  int add_idx,            ///< the idex to add the intensity -in
  FLOAT_T intensity         ///< the intensity to add -in
  );


BOOLEAN_T create_intensity_array_xcorr(
  SPECTRUM_T* spectrum,    ///< the spectrum to score(observed) -in
  SCORER_T* scorer,        ///< the scorer object -in/out
  int charge               ///< the peptide charge -in 
  );

*/

#endif
