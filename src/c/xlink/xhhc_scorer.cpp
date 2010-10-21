#include "xhhc_scorer.h"
#include "xhhc.h"

#include "Spectrum.h"

#include <fstream>

// constants for print_spectrums
#define NORMALIZE 0
#define MAX_MZ 1200;
#define MIN_MZ 400;
#define NO_FLANKS 1


Scorer::Scorer(FLOAT_T a_max_mz) {
  max_mz = a_max_mz;
}


int Scorer::get_matched_by_ions(Spectrum* spectrum,
				LinkedIonSeries& ion_series) {
  FLOAT_T bin_width = bin_width_mono;
  vector<LinkedPeptide>& ions = ion_series.ions();

  int ans = 0;

  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    if (ion -> get_mz(MONO) >= 400 && ion -> get_mz(MONO) <= 1200) {
    if (ion -> type() == B_ION || ion -> type() == Y_ION) {
      PEAK_T* peak = spectrum->getNearestPeak(ion->get_mz(AVERAGE), 
                                              bin_width);
      if (peak != NULL) {
	ans++;
      }
    }
  }
  }
  return ans;
}

float Scorer::score_spectrum_vs_series(
      Spectrum* spectrum,
      LinkedIonSeries& ion_series
  ) {
    //SCORER_T* scorer = new_scorer(XCORR);
    float score = 0.0;
    score =  hhc_gen_score_xcorr(spectrum, ion_series);
    return score; 
  } 
 
FLOAT_T Scorer::hhc_gen_score_xcorr(
    Spectrum* spectrum,    ///< the spectrum to score -in
    LinkedIonSeries& ion_series ///< the ion series to score against the spectrum -in
  )
  
  {

    //cout <<"hhc_gen_score_xcorr() - start."<<endl;

    FLOAT_T final_score = 0;
    FLOAT_T* theoretical = NULL;

    if (current_spectrum != spectrum) {
      //cout <<"Creating scorer"<<endl;
      current_spectrum = spectrum;
      if (scorer != NULL) { free_scorer(scorer); scorer=NULL;}
      scorer = new_scorer(XCORR);
      if (!create_intensity_array_xcorr(spectrum, scorer, ion_series.charge())) {
        carp(CARP_FATAL, "failed to produce XCORR");
      }
    } else {
      //cout <<"Using same scorer"<<endl;
    }

    max_mz = get_scorer_sp_max_mz(scorer);

    // create theoretical array
    //cout <<"Creating theoretical array"<<endl;
    theoretical = (FLOAT_T*)mycalloc((size_t)max_mz, sizeof(FLOAT_T));
  
    // create intensity array for theoretical spectrum 
    //if(!create_intensity_array_theoretical(scorer, ion_series, theoretical)){
   
   if (!hhc_create_intensity_array_theoretical(ion_series, theoretical)) {
      carp(CARP_ERROR, "failed to create theoretical spectrum for Xcorr");
      return FALSE;
    }
   
   //cout <<"Doing cross correlation"<<endl;
    // do cross correlation between observed spectrum(in scorer) and theoretical spectrum.
    // use the two intensity arrays that were created
    final_score = cross_correlation(scorer, theoretical);
    //cout <<"Done cross correlation"<<endl;
    if (print_spectrums_) {
      //cout <<"printing spectrums"<<endl;
      print_spectrums(theoretical, spectrum);
    }
    //print_spectrums(theoretical, spectrum, 300, 750, 1);

    // free theoretical spectrum
    //cout <<"freeing theoretical"<<endl;
    free(theoretical);
    //cout <<"hhc_gen_score_xcorr() - done."<<endl;
    // return score
    return final_score;
  }

void Scorer::add_intensity_map(std::map<int, FLOAT_T>& theoretical, int idx, FLOAT_T intensity) {
  
  std::map<int, FLOAT_T>::iterator iter = theoretical.find(idx);
  if (iter == theoretical.end())
    theoretical[idx] = intensity;
  else
    iter -> second = max(intensity, iter -> second);
}

bool Scorer::xlink_create_map_theoretical(
					  LinkedIonSeries& ion_series,
					  std::map<int, FLOAT_T>& theoretical) {
  theoretical.clear();
  //ION_T* ion = NULL;
  int ion_charge = 0;
  ION_TYPE_T ion_type;
  int intensity_array_idx = 0;
  FLOAT_T bin_width = bin_width_mono;
  vector<LinkedPeptide>& ions = ion_series.ions();
  // while there are ion's in ion iterator, add matched observed peak intensity
  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    intensity_array_idx = (int)(ion->get_mz(MONO) / bin_width + 0.5);
    ion_type = ion->type();
    ion_charge = ion->charge();

    // is it B, Y ion?
    
    // neutral loss peak?
    // Add peaks of intensity 50.0 for B, Y type ions. 
    // In addition, add peaks of intensity of 25.0 to +/- 1 m/z flanking each B, Y ion.
    // Skip ions that are located beyond max mz limit
    add_intensity_map(theoretical, intensity_array_idx, 50);
    if (get_boolean_parameter("xcorr-use-flanks")) {
      if (intensity_array_idx > 0) {
	add_intensity_map(theoretical, intensity_array_idx - 1, 25);
      }
      add_intensity_map(theoretical, intensity_array_idx + 1, 25);
     
    }

    // add neutral loss of water and NH3
    //  mass_z + (modification_masses[(int)ion_modification]/(FLOAT_T)charge) * modification_count;  
    
    //TODO put this back in
    //if(ion_type == B_ION){
      int h2o_array_idx = (int)((ion->get_mz(MONO) - (MASS_H2O_MONO/ion->charge()) ) / bin_width + 0.5);
      add_intensity_map(theoretical, h2o_array_idx, 10);
    //}
    
    int co_array_idx = (int)((ion -> get_mz(MONO) - (MASS_CO_MONO/ion->charge())) / bin_width + 0.5);
    add_intensity_map(theoretical, co_array_idx, 10);

    int nh3_array_idx = (int)((ion->get_mz(MONO) -  (MASS_NH3_MONO/ion->charge())) / bin_width + 0.5);
    add_intensity_map(theoretical, nh3_array_idx, 10);        
  }
  return true;
}


bool Scorer::hhc_create_intensity_array_theoretical(
    LinkedIonSeries& ion_series,
    FLOAT_T* theoretical       ///< the empty theoretical spectrum -out
    )
  {
    //ION_T* ion = NULL;
    int ion_charge = 0;
    ION_TYPE_T ion_type;
    int intensity_array_idx = 0;
    FLOAT_T bin_width = bin_width_mono;
    vector<LinkedPeptide>& ions = ion_series.ions();
    // while there are ion's in ion iterator, add matched observed peak intensity
    for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    //while(ion_iterator_has_next(ion_iterator)){
      //cout << "ion " << *ion << "\tmz " << ion->get_mz() << endl;
      intensity_array_idx = (int)(ion->get_mz(MONO) / bin_width + 0.5);
      //if (intensity_array_idx <= 1) continue;
      //cout << "index " << intensity_array_idx << endl;
      ion_type = ion->type();
      ion_charge = ion->charge();
      //cout << "m/z: " << ion->get_mz() << " charge: " << ion->charge() << endl;
      // skip ions that are located beyond max mz limit

      // is it B, Y ion?

      // neutral loss peak?
      // Add peaks of intensity 50.0 for B, Y type ions. 
      // In addition, add peaks of intensity of 25.0 to +/- 1 m/z flanking each B, Y ion.
      // Skip ions that are located beyond max mz limit
      if((intensity_array_idx)< max_mz){
        //if (ion->type() == Y_ION)
        //add_intensity(theoretical, intensity_array_idx, 51);
        add_intensity(theoretical, intensity_array_idx, 50);
        if (get_boolean_parameter("xcorr-use-flanks") &&
	    intensity_array_idx > 0) {
	      add_intensity(theoretical, intensity_array_idx - 1, 25);
      }
	if(get_boolean_parameter("xcorr-use-flanks") &&
	   ((intensity_array_idx + 1)< max_mz)) {
        add_intensity(theoretical, intensity_array_idx + 1, 25);
	}
      }

      // add neutral loss of water and NH3
      //  mass_z + (modification_masses[(int)ion_modification]/(FLOAT_T)charge) * modification_count;  


      if(ion_type == B_ION){
        int h2o_array_idx = (int)((ion->get_mz(MONO) - (MASS_H2O_MONO/ion->charge()) ) / bin_width + 0.5);
	if (h2o_array_idx < max_mz)
	  add_intensity(theoretical, h2o_array_idx, 10);
      }

      int nh3_array_idx = (int)((ion->get_mz(MONO) -  (MASS_NH3_MONO/ion->charge())) / bin_width + 0.5);
      if (nh3_array_idx < max_mz)
	add_intensity(theoretical, nh3_array_idx, 10);        
    }
    return true;
  }


FLOAT_T Scorer::getIonCurrentExplained(LinkedIonSeries& ion_series, 
  Spectrum* spectrum, 
  FLOAT_T& explained, 
  int& by_observed) {

  max_mz = spectrum->getMaxPeakMz();
  FLOAT_T* theoretical = (FLOAT_T*)mycalloc((size_t)max_mz+1, sizeof(FLOAT_T));
  hhc_create_intensity_array_theoretical(ion_series, theoretical);

  explained = 0.0;
  by_observed = 0;

  FLOAT_T bin_width = bin_width_mono;

  FLOAT_T ans = 0.0;

  map<int, bool> by_found;

  for (PeakIterator peak_iter = spectrum->begin();
    peak_iter != spectrum->end();
    ++peak_iter) {

    PEAK_T* peak = *peak_iter;
    FLOAT_T spec_mz = get_peak_location(peak);
    FLOAT_T spec_intensity = get_peak_intensity(peak);
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



//creates three files for spectacle.pl: spectrums.out, theoretical.out, observed.out
void Scorer::print_spectrums(FLOAT_T* theoretical, Spectrum* spectrum) {
   
  ofstream theoretical_file;
  ofstream observed_file;
  ofstream spectrums_file;

  theoretical_file.open("theoretical.out");
  observed_file.open("observed.out");
  spectrums_file.open("spectrums.out");

  theoretical_file << "> theoretical" << endl;
  observed_file << "> observed" << endl;
  spectrums_file << "> spectrums" << endl;
  bool noflanks = get_boolean_parameter("xcorr-use-flanks");
  int normalize = NORMALIZE;
  int max_mz = MAX_MZ;
  int min_mz = MIN_MZ;
  // keep track of colors
  map<PEAK_T*, string> peak_colors;
  carp(CARP_DEBUG, "min mz: %d, max mz: %d\n", max_mz);
  FLOAT_T average = 0;

  for (PeakIterator peak_iter = spectrum->begin();
    peak_iter != spectrum->end();
    ++peak_iter) {

    average += get_peak_intensity(*peak_iter);
  }

  average = average / spectrum->getNumPeaks();
  //cout << "AVERAGE " << average << endl;
  // make spectacle file for observed peaks
 
  for (PeakIterator peak_iter = spectrum->begin();
    peak_iter != spectrum->end();
    ++peak_iter) {

    PEAK_T* peak = *peak_iter;
    FLOAT_T location = get_peak_location(peak);
    FLOAT_T intensity = get_peak_intensity(peak); 
    if (location > min_mz && location < max_mz) {
    if (normalize) {
      peak_colors[peak] = "blue";
      //observed_file << location<< "\t" << pow(intensity * average * normalize, 0.2) << "\tnolabel\tred" << endl;
      //spectrums_file << location<< "\t" << pow(intensity * average * normalize, 0.2) << "\tnolabel\tblue" << endl;
      //spectrums_file << location<< "\t" << pow(intensity * average * normalize, 0.25) << "\tnolabel" << endl;
    } else {
      observed_file << location<< "\t" << intensity << "\tnolabel\tred" << endl;
      spectrums_file << location<< "\t" << intensity  << "\tnolabel\tblue" << endl;
      //spectrums_file << location<< "\t" << intensity  << "\tnolabel" << endl;
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
      PEAK_T* peak = spectrum->getNearestPeak(i, 1);
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
  for (map<PEAK_T*, string>::iterator it = peak_colors.begin(); it != peak_colors.end(); ++it) {
    location = get_peak_location(it->first);
    intensity = get_peak_intensity(it->first);
    //spectrums_file << location << "\t" << pow(intensity * average * normalize, 0.2) << "\tnolabel\t" << it->second << endl;
    spectrums_file << location << "\t" << intensity << "\tnolabel\t" << it->second << endl;
  }

  //cout << "match: " << match_count << " mismatch: " << mismatch_count << endl;
  theoretical_file.close();
  spectrums_file.close();
}

