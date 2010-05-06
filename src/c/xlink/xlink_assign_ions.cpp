#include "xhhc.h"
#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"

#include "objects.h"
#include "scorer.h"
#include "spectrum_collection.h"
#include "DelimitedFile.h"

#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
//#define bin_width_mono 1.0005079


#define NUM_ARGUMENTS 8
#define NUM_OPTIONS 5


void print_spectrum(SPECTRUM_T* spectrum, LinkedIonSeries& ion_series);
int main(int argc, char** argv){

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  int num_options = NUM_OPTIONS;
  const char* option_list[NUM_OPTIONS] = {
    "verbosity",
    "version",
    "use-mgf",
    "ion-tolerance",
    "precision"
  };

  

  /* Define required command line arguments */
  int num_arguments = NUM_ARGUMENTS;
  const char* argument_list[NUM_ARGUMENTS] = {"peptide A",
                                              "peptide B",
					      "pos A",
					      "pos B",
					      "link mass",
					      "charge state",
					      "scan number",
					      "ms2 file"};

  
  /* for debugging of parameter processing */
  //set_verbosity_level( CARP_DETAILED_DEBUG );
  set_verbosity_level( CARP_ERROR );
  
  /* Set default values for parameters in parameter.c */
  initialize_parameters();

  /* Define optional and required command line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  parse_cmd_line_into_params_hash(argc, argv, "xlink-assign-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
  char* peptideA = get_string_parameter("peptide A");
  char* peptideB = get_string_parameter("peptide B");
  
  int posA     = get_int_parameter("pos A");
  int posB     = get_int_parameter("pos B");
  int charge   = get_int_parameter("charge state"); 
  int scan_num = get_int_parameter("scan number"); 

  char* ms2_file = get_string_parameter("ms2 file");

  LinkedPeptide::linker_mass = get_double_parameter("link mass");
 
  // create new ion series
  
  // a single peptide linked to itself
  if (strcmp(peptideB, "NULL") == 0) {
    cout << "B is null" << endl; 
    peptideB = NULL;
  }

  // read ms2 file
  SPECTRUM_COLLECTION_T* collection = new_spectrum_collection(ms2_file);
  parse_spectrum_collection(collection);

  
  //cout << "lp " << lp << endl; 

  SPECTRUM_ITERATOR_T* spectrum_iterator = new_spectrum_iterator(collection);
  SPECTRUM_T* spectrum = NULL;
  SPECTRUM_T* current_spectrum = NULL;

  //TODO allow a binary search on both mgf and ms2 files.

  while (spectrum_iterator_has_next(spectrum_iterator)) {
    current_spectrum = spectrum_iterator_next(spectrum_iterator);
    if (get_spectrum_first_scan(current_spectrum) == scan_num) {
      spectrum = current_spectrum;
      break;
    }
  }
  free(spectrum_iterator);



  // search for spectrum with correct scan number
  if (spectrum == NULL) {
    carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
    free_spectrum_collection(collection);
    exit(1);
  }

  //created linked peptide.
  LinkedPeptide lp = LinkedPeptide(peptideA, peptideB, posA-1, posB-1, charge);

  LinkedIonSeries ion_series(charge);
  ion_series.add_linked_ions(lp);
  print_spectrum(spectrum, ion_series);

  // free heap
  free_spectrum_collection(collection);
  //free_spectrum(spectrum);
}

void print_spectrum(SPECTRUM_T* spectrum, LinkedIonSeries& ion_series) {

      int total_by_ions = ion_series.get_total_by_ions();
      int matched_by_ions = Scorer::get_matched_by_ions(spectrum, ion_series);
      FLOAT_T frac_by_ions = (double)matched_by_ions / (double) total_by_ions;

      carp(CARP_INFO, "total theoretical ions:%d",total_by_ions);
      carp(CARP_INFO,"theoretical ions matched:%d",matched_by_ions);
      carp(CARP_INFO,"frac theoretical ions matched:%f",frac_by_ions);
      carp(CARP_INFO,"npeaks:%d",get_spectrum_num_peaks(spectrum));

      FLOAT_T bin_width = get_double_parameter("ion-tolerance");
      vector<LinkedPeptide>& ions = ion_series.ions();
      
      map<PEAK_T*, LinkedPeptide> matched;
      double matched_intensity = 0;
      for (vector<LinkedPeptide>::iterator ion = ions.begin();
	   ion != ions.end(); 
	   ++ion) {
        
	//if (ion -> get_mz(MONO) >= 400 && ion -> get_mz(MONO) <= 1400) {
	  if (ion -> type() == B_ION || ion -> type() == Y_ION) {
	    PEAK_T* peak = get_nearest_peak(spectrum, ion -> get_mz(MONO), bin_width);
	    if (peak != NULL) {
              if (matched.find(peak) == matched.end()) {
                matched_intensity += get_peak_intensity(peak);
              }
	      matched[peak] = *ion;
	    } else {
              carp(CARP_DETAILED_DEBUG,"Ion clash!");
            }
	  }
	//}
      }

      double total_intensity = get_spectrum_total_energy(spectrum);
      double frac_intensity = matched_intensity / total_intensity;

      carp(CARP_INFO,"matched intensity:%lf",matched_intensity);
      carp(CARP_INFO,"total intensity:%lf",total_intensity);
      carp(CARP_INFO,"frac intensity:%lf",frac_intensity);
      

      //now print out the spectrum

      DelimitedFile result_file;

      unsigned int mz_obs_col = result_file.addColumn("m/z obs");
      unsigned int int_col = result_file.addColumn("intensity");
      unsigned int mz_calc_col = result_file.addColumn("m/z calc");
      unsigned int mass_obs_col = result_file.addColumn("mass obs");
      unsigned int mass_calc_col = result_file.addColumn("mass calc");
      unsigned int ppm_col = result_file.addColumn("ppm(Mass)");
      unsigned int ion_col = result_file.addColumn("ion");
      unsigned int seq_col = result_file.addColumn("sequence");
      

      PEAK_ITERATOR_T* peak_iter = new_peak_iterator(spectrum);
      while (peak_iterator_has_next(peak_iter)) {
	PEAK_T* peak = peak_iterator_next(peak_iter);
        
	//if (get_peak_location(peak) >= 400 && get_peak_location(peak) <= 1400) {
          unsigned int row_idx = result_file.addRow();
          result_file.setValue(mz_obs_col, row_idx, get_peak_location(peak));
          result_file.setValue(int_col, row_idx, get_peak_intensity(peak));

          

	  if (matched.find(peak) != matched.end()) {
              LinkedPeptide& ion = matched[peak];
              double mz_calc = ion.get_mz(MONO);
              double mz_obs = get_peak_location(peak);
              int charge = ion.charge();
              double mass_calc = ion.mass(MONO);
              double mass_obs = (mz_obs - MASS_PROTON) * (double)charge;
              double ppm =  fabs(mass_calc - mass_obs) / mass_calc * 1e6;
              
              result_file.setValue(mz_calc_col, row_idx, mz_calc);
              result_file.setValue(mass_obs_col, row_idx, mass_obs);
              result_file.setValue(mass_calc_col, row_idx, mass_calc);
              result_file.setValue(ppm_col, row_idx, ppm);
              
              ostringstream ion_string_oss;
           


              if (ion.type() == B_ION) {
                ion_string_oss << "b";
              } else if (ion.type() == Y_ION) {
                ion_string_oss << "y";
              } else {
                ion_string_oss << "u";
              }
              ion_string_oss << "(+" << charge << ")";

              result_file.setValue(ion_col, row_idx, ion_string_oss.str());
              result_file.setValue(seq_col, row_idx, ion);
          }
        //}
      }
      free_peak_iterator(peak_iter);

      cout << result_file;

}
