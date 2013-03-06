/**
 * \file xhhc_scorer.h 
 * \brief object to score spectrum vs. xlink peptides
 */
/*
 * AUTHOR: Sean McIlwain
 * $Revision: 1.0 $
 *****************************************************************************/
#ifndef XHHC_SCORER_H
#define XHHC_SCORER_H

#include "Scorer.h"
#include "XHHC_Peptide.h"

#include "LinkedIonSeries.h"

#include <map>

class XHHC_Scorer {
 private:
  bool print_spectrums_; ///< Indicator of whether to print the spectrums
  Scorer* scorer_; ///< The crux scorer
  Crux::Spectrum* current_spectrum_; ///< Spectrum that this scorer works on
  FLOAT_T max_mz_; ///< Maximum mz allocated for scorer
  int max_bin_; ///< Maximum bin idx for scorer.
  FLOAT_T bin_width_; ///< bin width for scorer
  FLOAT_T bin_offset_; ///<bin offset for scorer

  static int ion_counter_; ///< Ion counter

  /**
   * Prints out an annotated spectrum
   * creates three files for spectacle.pl: spectrums.out, 
   * theoretical.out, observed.out
   */
  void printSpectrums(
    FLOAT_T* theoretical, ///< The theoretical spectrum array
    Crux::Spectrum* spectrum ///< The spectrum to print
    );
  
  /**
   * Initializes an empty XHHC_Scorer object
   */
  void init();

 public:

  /**
   * \returns an empty XHHC_Scorer object
   */
  XHHC_Scorer();

  /**
   * \returns an XHHC_Scorer object with the max_mz initialized
   */
  XHHC_Scorer(
    FLOAT_T max_mz ///< max_mz for the scorer.
    );

  /**
   * Destructor
   */
  virtual ~XHHC_Scorer();

  /**
   * \returns the max_mz
   */
  FLOAT_T getMaxMz();

  /**
   * \returns the xcorr score for the spectrum and the ion_series
   */
  FLOAT_T scoreSpectrumVsSeries(
    Crux::Spectrum* spectrum, ///< Spectrum to score
    LinkedIonSeries& ion_series ///< LinkedIonSeries to score
    );

  /**
   * \returns the number of BY ions in the Ionseries that are
   * matched to peaks in the spectrum
   */
  static int getMatchedBYIons(
    Crux::Spectrum* spectrum, ///< Spectrum to match peaks with
    LinkedIonSeries& ion_series ///< IonSeries to use
    );

  /**
   * adds an intensity in the theoretical map at peak idx
   */
  static void addIntensityMap(
    std::map<int, FLOAT_T>& theoretical, ///< the theoretical map 
    int idx, ///< the idx to add
    FLOAT_T intensity ///< the corresponding intensity
    );

  /**
   * \returns whether create a theoretical map of peak_idx,intensity pairs
   * is sucessful
   */
  static bool xlinkCreateMapTheoretical(
    LinkedIonSeries& ion_series, ///< LinkedIonSeries to create the map from -in
    std::map<int, FLOAT_T>& theoretical ///< The theoretical map -out
  );

  /**
   * Fills in the theoretical array using the LinkedIonSeries
   */
  bool hhcCreateIntensityArrayTheoretical(
    LinkedIonSeries& ion_series, /// The linked ion series to use -in
    FLOAT_T* theoretical       ///< the empty theoretical spectrum -out
    );
  
  /*
   * Sets the print_spectrums_ variable
   */
  void setPrint(
    bool print_spectrums ///< indicator of whether we are printing spectrums
    ); 

  /**
   * Generates the theoretical/observed arrays and scores the spectrum for xcorr
   */
  FLOAT_T hhcGenScoreXcorr(
    Crux::Spectrum* spectrum,    ///< the spectrum to score -in
    LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
    );

  /**
   * \returns the sum of the spectrum peak intensities of the by-ions that are matched
   * to the LinkedIonSeries
   */ 
  FLOAT_T getIonCurrentExplained(
    LinkedIonSeries& ion_series, ///<The LinkedIonSeries to match 
    Crux::Spectrum* spectrum, ///<The spectrum to match
    FLOAT_T& explained, ///<The ion current explained -out
    int& by_observed ///<The number of by ions matched
    );
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
