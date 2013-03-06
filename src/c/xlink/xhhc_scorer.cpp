#include "xhhc_scorer.h"
#include "xhhc.h"
#include "LinkedPeptide.h"

#include "Spectrum.h"
#include "Scorer.h"

#include <fstream>

// constants for print_spectrums
#define NORMALIZE 0
#define MAX_MZ 1200;
#define MIN_MZ 400;
#define NO_FLANKS 1

using namespace Crux;

/**
 * Initializes an empty XHHC_Scorer object
 */
void XHHC_Scorer::init() {
  print_spectrums_ = false;
  scorer_ = NULL;
  max_mz_ = 0.0;
  max_bin_ = 0;
  current_spectrum_ = NULL;

  bin_width_ = get_double_parameter("mz-bin-width");
  bin_offset_ = get_double_parameter("mz-bin-offset");

}

/**
 * \returns an empty XHHC_Scorer object
 */
XHHC_Scorer::XHHC_Scorer() {
  init();
}

/**
 * Destructor
 */
XHHC_Scorer::~XHHC_Scorer() {

  if (scorer_ != NULL) {
    delete scorer_;
  }

}
  
/**
 * \returns an XHHC_Scorer object with the max_mz initialized
 */
XHHC_Scorer::XHHC_Scorer(
  FLOAT_T max_mz ///< max_mz for the scorer.
  ) {

  init();
  max_mz_ = max_mz;
}

/**
 * \returns the max_mz
 */
FLOAT_T XHHC_Scorer::getMaxMz() {
  
  return max_mz_;
}

/*
 * Sets the print_spectrums_ variable
 */
void XHHC_Scorer::setPrint(
  bool print_spectrums ///< indicator of whether we are printing spectrums
  ) { 

  print_spectrums_ = print_spectrums;
}

/**
 * \returns the number of BY ions in the Ionseries that are
 * matched to peaks in the spectrum
 */
int XHHC_Scorer::getMatchedBYIons(
  Spectrum* spectrum,
  LinkedIonSeries& ion_series) {

  FLOAT_T bin_width = bin_width_mono;
  vector<LinkedPeptide>& ions = ion_series.getIons();

  int ans = 0;

  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    if (ion -> getMZ(MONO) >= 400 && ion -> getMZ(MONO) <= 1200) {
    if (ion -> getIonType() == B_ION || ion -> getIonType() == Y_ION) {
      Peak * peak = spectrum->getNearestPeak(ion->getMZ(AVERAGE), 
                                              bin_width);
      if (peak != NULL) {
	ans++;
      }
    }
  }
  }
  return ans;
}

/**
 * \returns the xcorr score for the spectrum and the ion_series
 */
FLOAT_T XHHC_Scorer::scoreSpectrumVsSeries(
  Spectrum* spectrum, ///< Spectrum to score
  LinkedIonSeries& ion_series ///< LinkedIonSeries to score
  ) {

  FLOAT_T score = 0.0;
  score =  hhcGenScoreXcorr(spectrum, ion_series);
  return score; 
} 
 
/**
 * Generates the theoretical/observed arrays and scores the spectrum for xcorr
 */
FLOAT_T XHHC_Scorer::hhcGenScoreXcorr(
  Spectrum* spectrum,    ///< the spectrum to score -in
  LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
  ) {

  FLOAT_T final_score = 0;
  FLOAT_T* theoretical = NULL;

  if (current_spectrum_ != spectrum) {
    current_spectrum_ = spectrum;
    if (scorer_ != NULL) { 
      delete scorer_; 
      scorer_=NULL;
    }
    scorer_ = new Scorer(XCORR);
    if (!scorer_->createIntensityArrayXcorr(spectrum, ion_series.getCharge())) {
      carp(CARP_FATAL, "failed to produce XCORR");
    }
  }

  max_mz_ = scorer_->getSpMaxMz();
  max_bin_ = scorer_->getMaxBin();
  // create theoretical array
  theoretical = (FLOAT_T*)mycalloc((size_t)max_bin_, sizeof(FLOAT_T));
  
 if (!hhcCreateIntensityArrayTheoretical(ion_series, theoretical)) {
    carp(CARP_ERROR, "failed to create theoretical spectrum for Xcorr");
    return false;
  }
   
  // do cross correlation between observed spectrum(in scorer) and theoretical spectrum.
  // use the two intensity arrays that were created
  final_score = scorer_->crossCorrelation(theoretical);
  if (print_spectrums_) {
    printSpectrums(theoretical, spectrum);
  }
  // free theoretical spectrum
  free(theoretical);
  return final_score;
}

/**
 * adds an intensity in the theoretical map at peak idx
 */
void XHHC_Scorer::addIntensityMap(
  map<int, FLOAT_T>& theoretical, ///< the theoretical map 
  int idx, ///< the idx to add   
  FLOAT_T intensity ///< the corresponding intensity
  ) {
  
  map<int, FLOAT_T>::iterator iter = theoretical.find(idx);
  if (iter == theoretical.end())
    theoretical[idx] = intensity;
  else
    iter -> second = max(intensity, iter -> second);
}


/**
 * \returns whether create a theoretical map of peak_idx,intensity pairs
 * is sucessful
 */
bool XHHC_Scorer::xlinkCreateMapTheoretical(
  LinkedIonSeries& ion_series, ///< LinkedIonSeries to create the map from -in
  map<int, FLOAT_T>& theoretical ///< The theoretical map -out
  ) {

  theoretical.clear();
  int intensity_array_idx = 0;
  FLOAT_T bin_width = bin_width_mono;
  vector<LinkedPeptide>& ions = ion_series.getIons();
  // while there are ion's in ion iterator, add matched observed peak intensity
  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    intensity_array_idx = (int)(ion->getMZ(MONO) / bin_width + 0.5);
    //ion_type = ion->getIonType();
    //ion_charge = ion->getCharge();

    // is it B, Y ion?
    
    // neutral loss peak?
    // Add peaks of intensity 50.0 for B, Y type ions. 
    // In addition, add peaks of intensity of 25.0 to +/- 1 m/z flanking each B, Y ion.
    // Skip ions that are located beyond max mz limit
    addIntensityMap(theoretical, intensity_array_idx, 50);
    if (get_boolean_parameter("use-flanking-peaks")) {
      if (intensity_array_idx > 0) {
	addIntensityMap(theoretical, intensity_array_idx - 1, 25);
      }
      addIntensityMap(theoretical, intensity_array_idx + 1, 25);
     
    }

    // add neutral loss of water and NH3
    //  mass_z + (modification_masses[(int)ion_modification]/(FLOAT_T)charge) * modification_count;  
    
    //TODO put this back in
    //if(ion_type == B_ION){
      int h2o_array_idx = (int)((ion->getMZ(MONO) - (MASS_H2O_MONO/ion->getCharge()) ) / bin_width + 0.5);
      addIntensityMap(theoretical, h2o_array_idx, 10);
    //}
    
    int co_array_idx = (int)((ion -> getMZ(MONO) - (MASS_CO_MONO/ion->getCharge())) / bin_width + 0.5);
    addIntensityMap(theoretical, co_array_idx, 10);

    int nh3_array_idx = (int)((ion->getMZ(MONO) -  (MASS_NH3_MONO/ion->getCharge())) / bin_width + 0.5);
    addIntensityMap(theoretical, nh3_array_idx, 10);        
  }
  return true;
}

/**
 * Fills in the theoretical array using the LinkedIonSeries
 */
bool XHHC_Scorer::hhcCreateIntensityArrayTheoretical(
  LinkedIonSeries& ion_series, ///< The linked ion series to use -in
  FLOAT_T* theoretical       ///< the empty theoretical spectrum -out
  ) {

  ION_TYPE_T ion_type;
  int intensity_array_idx = 0;
  vector<LinkedPeptide>& ions = ion_series.getIons();
  // while there are ion's in ion iterator, add matched observed peak intensity
  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    intensity_array_idx = INTEGERIZE(ion->getMZ(MONO), bin_width_, bin_offset_);
    ion_type = ion->getIonType();
    //ion_charge = ion->getCharge();
    // skip ions that are located beyond max mz limit
    // is it B, Y ion?
    // neutral loss peak?
    // Add peaks of intensity 50.0 for B, Y type ions. 
    // In addition, add peaks of intensity of 25.0 to +/- 1 m/z flanking each B, Y ion.
    // Skip ions that are located beyond max mz limit
    if((intensity_array_idx)< max_bin_){
      //if (ion->type() == Y_ION)
      //add_intensity(theoretical, intensity_array_idx, 51);
      Scorer::addIntensity(theoretical, intensity_array_idx, 50);
      if (get_boolean_parameter("use-flanking-peaks") &&
        intensity_array_idx > 0) {
        Scorer::addIntensity(theoretical, intensity_array_idx - 1, 25);
      }
      if(get_boolean_parameter("use-flanking-peaks") &&
	 ((intensity_array_idx + 1)< max_bin_)) {
        Scorer::addIntensity(theoretical, intensity_array_idx + 1, 25);
      }
    }

    // add neutral loss of water and NH3
    //  mass_z + (modification_masses[(int)ion_modification]/(FLOAT_T)charge) * modification_count;  


    if(ion_type == B_ION){
      int h2o_array_idx = 
        INTEGERIZE((ion->getMZ(MONO) - (MASS_H2O_MONO/ion->getCharge())), bin_width_, bin_offset_);

      if (h2o_array_idx < max_bin_)
	Scorer::addIntensity(theoretical, h2o_array_idx, 10);
      }

      int nh3_array_idx = 
        INTEGERIZE((ion->getMZ(MONO) -  (MASS_NH3_MONO/ion->getCharge())),
		   bin_width_, bin_offset_);
      if (nh3_array_idx < max_bin_) {
        Scorer::addIntensity(theoretical, nh3_array_idx, 10);
      }
  }
  return true;
}

/**
 * \returns the sum of the spectrum peak intensities of the by-ions that are matched
 * to the LinkedIonSeries
 */ 
FLOAT_T XHHC_Scorer::getIonCurrentExplained(
  LinkedIonSeries& ion_series, ///<The LinkedIonSeries to match 
  Spectrum* spectrum, ///<The spectrum to match
  FLOAT_T& explained, ///<The ion current explained -out
  int& by_observed ///<The number of by ions matched
  ) {

  max_mz_ = spectrum->getMaxPeakMz();
  FLOAT_T* theoretical = (FLOAT_T*)mycalloc((size_t)max_mz_+1, sizeof(FLOAT_T));
  hhcCreateIntensityArrayTheoretical(ion_series, theoretical);

  explained = 0.0;
  by_observed = 0;

  FLOAT_T bin_width = bin_width_mono;

  FLOAT_T ans = 0.0;

  map<int, bool> by_found;

  for (PeakIterator peak_iter = spectrum->begin();
    peak_iter != spectrum->end();
    ++peak_iter) {

    Peak * peak = *peak_iter;
    FLOAT_T spec_mz = peak->getLocation();
    FLOAT_T spec_intensity = peak->getIntensity();
    //for each peak in the spectrum, find the array index for the theoretical.
    int intensity_array_idx = (int)(spec_mz / bin_width + 0.5);

    if (theoretical[intensity_array_idx] >= 45) {
    //it is a b-y ions that is observed.
      if (by_found.find(intensity_array_idx) == by_found.end()) {
        by_found[intensity_array_idx] = true;
        by_observed ++;
      }
      ans += spec_intensity;
    }
  }
  free(theoretical);


  explained = ans;
  return ans;
}

/**
 * Prints out an annotated spectrum
 * creates three files for spectacle.pl: spectrums.out, 
 * theoretical.out, observed.out
 */
void XHHC_Scorer::printSpectrums(
  FLOAT_T* theoretical, ///< The theoretical spectrum array
  Spectrum* spectrum ///< The spectrum to print
  ){
   
  ofstream theoretical_file;
  ofstream observed_file;
  ofstream spectrums_file;

  theoretical_file.open("theoretical.out");
  observed_file.open("observed.out");
  spectrums_file.open("spectrums.out");

  theoretical_file << "> theoretical" << endl;
  observed_file << "> observed" << endl;
  spectrums_file << "> spectrums" << endl;
  bool noflanks = get_boolean_parameter("use-flanking-peaks");
  int normalize = NORMALIZE;
  int max_mz = MAX_MZ;
  int min_mz = MIN_MZ;
  // keep track of colors
  map<Peak *, string> peak_colors;
  carp(CARP_DEBUG, "min mz: %d, max mz: %d\n", max_mz);
  FLOAT_T average = 0;

  for (PeakIterator peak_iter = spectrum->begin();
    peak_iter != spectrum->end();
    ++peak_iter) {
    
    average += (*peak_iter)->getIntensity();
  }

  average = average / spectrum->getNumPeaks();
  // make spectacle file for observed peaks
 
  for (PeakIterator peak_iter = spectrum->begin();
    peak_iter != spectrum->end();
    ++peak_iter) {

    Peak * peak = *peak_iter;
    FLOAT_T location = peak->getLocation();
    FLOAT_T intensity = peak->getIntensity(); 
    if (location > min_mz && location < max_mz) {
      if (normalize) {
        peak_colors[peak] = "blue";
      } else {
        observed_file << location<< "\t" << intensity << "\tnolabel\tred" << endl;
        spectrums_file << location<< "\t" << intensity  << "\tnolabel\tblue" << endl;
      }
    }
  }
  

  observed_file.close();
  // make spectacle file for theoretical peaks
  FLOAT_T* index = theoretical;
  int i = 0;
  int match_count = 0;
  int mismatch_count = 0;
  while (i <= max_mz)  {
    if (((*index > 1 && !noflanks) || *index > 26) && i >= min_mz) {
        theoretical_file << i << "\t" << *index << "\tnolabel\tred" << endl;
      Peak * peak = spectrum->getNearestPeak(i, 1);
      if (peak != NULL) {
	++match_count;
	peak_colors[peak] = "green";
	spectrums_file << i << "\t" << -10000 << "\tnolabel\tgreen" << endl;
	//spectrums_file << get_peak_location(peak) << "\t" << pow (get_peak_intensity(peak) * average * normalize, 0.2) << "\tnolabel\tgreen" << endl;
      } else {
	++mismatch_count;
        spectrums_file << i << "\t" << -10000 << "\tnolabel\tred" << endl;
      }
    }
    ++i;
    ++index;
  }
  FLOAT_T location;
  FLOAT_T intensity;
  for (map<Peak *, string>::iterator it = peak_colors.begin(); it != peak_colors.end(); ++it) {
    location = it->first->getLocation();
    intensity = it->first->getIntensity();;
    //spectrums_file << location << "\t" << pow(intensity * average * normalize, 0.2) << "\tnolabel\t" << it->second << endl;
    spectrums_file << location << "\t" << intensity << "\tnolabel\t" << it->second << endl;
  }

  theoretical_file.close();
  spectrums_file.close();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
