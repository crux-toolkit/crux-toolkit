#include "xhhc_scorer.h"
#include "xhhc_ion_series.h"
//#include "xhhc_search.h"
#include "xlink_compute_qvalues.h"


//CRUX INCLUDES
#include "objects.h"
#include "Spectrum.h"
#include "FilteredSpectrumChargeIterator.h"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <ctime>

using namespace std;

//typedef map<char, set<char> > BondMap;

double bonf_correct(double nlp_value, int nt);



void get_ions_from_window(vector<LinkedPeptide>& filtered_ions,
  vector<LinkedPeptide>& all_ions,
  FLOAT_T precursor_mass,
  FLOAT_T window,
  WINDOW_TYPE_T window_type
);

void get_ions_from_mass_range(
  vector<LinkedPeptide>& filtered_ions,
  vector<LinkedPeptide>& all_ions,
  double min_mass,
  double max_mass
);

void get_ions_from_mz_range(vector<LinkedPeptide>& filtered_ions,
	vector<LinkedPeptide>& all_ions,
	FLOAT_T precursor_mass,
	int charge,
	FLOAT_T mass_window,
	int decoy_iterations);

void plot_weibull(vector<pair<FLOAT_T, LinkedPeptide> >& scores, 
                  Spectrum* spectrum, int charge); 


string get_protein_ids_locations(vector<PEPTIDE_T*>& peptides);

int xlink_search_main(int argc, char** argv) {

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);

  /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
    "precursor-window",
    "precursor-window-type",
    "precursor-window-decoy",
    "precursor-window-type-decoy",
    "max-ion-charge",
    "min-weibull-points",
    /* TODO: Implement or remove
    "missed-link-cleavage",
    */
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "top-match",
    "xlink-include-linears",
    "xlink-include-deadends",
    "xlink-include-selfloops",
    "xcorr-use-flanks",
    "use-mgf"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {
    "ms2 file", 
    "protein database", 
    "link sites", 
    "link mass"
  };

  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize_run(XLINK_SEARCH_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);
  
  carp(CARP_INFO, "Beginning crux xlink-search");


  //int num_missed_cleavages = 0;
  char* ms2_file = get_string_parameter("ms2 file");

  FLOAT_T precursor_window = get_double_parameter("precursor-window");
  FLOAT_T precursor_window_decoy = get_double_parameter("precursor-window-decoy");
  WINDOW_TYPE_T precursor_window_type = 
    get_window_type_parameter("precursor-window-type");
  WINDOW_TYPE_T window_type_decoy = 
    get_window_type_parameter("precursor-window-type-decoy");



  char* database = get_string_parameter("protein database");
  char* links = get_string_parameter("link sites");

  unsigned int min_weibull_points = 
    (unsigned int)get_int_parameter("min-weibull-points");

  int scan_num = 0;
  //int charge = 1;
  SpectrumZState zstate;
  int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  int top_match = get_int_parameter("top-match");

  FLOAT_T linker_mass = get_double_parameter("link mass");

  //MASS_Type_T 

  LinkedPeptide::linker_mass = linker_mass;
  vector<LinkedPeptide> all_ions;
  carp(CARP_DETAILED_DEBUG,"Calling find all precursor ions");
  carp(CARP_INFO, "Building xlink database");
  find_all_precursor_ions(all_ions, links, "K", database,1);
  carp(CARP_DETAILED_DEBUG,"Sort");
  // sort filtered ions and decoy ions by mass
  //sort(all_ions.begin(), all_ions.end());


  carp(CARP_INFO,"Loading Spectra");
  Spectrum* spectrum = new Spectrum();
  SpectrumCollection* spectra = new SpectrumCollection(ms2_file);
  spectra->parse();

  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);
 
  FLOAT_T score;
 // best pvalues

  char* output_directory = get_string_parameter("output-dir");
  const char* target_filename = "search.target.txt";
  
  string target_path = string(output_directory) + "/" + string(target_filename);
  ofstream search_target_file(target_path.c_str());
  //set precision
  search_target_file << setprecision(get_int_parameter("precision"));
  //print header
  search_target_file << "scan\t";
  search_target_file << "charge\t";
  search_target_file << "spectrum precursor m/z\t";
  search_target_file << "spectrum neutral mass\t";
  search_target_file << "peptide mass mono\t";
  search_target_file << "peptide mass average\t";
  search_target_file << "mass error(ppm)\t";
  search_target_file << "xcorr score\t";
  search_target_file << "xcorr rank\t";
  search_target_file << "p-value\t";
  search_target_file << "matches/spectrum\t";
  search_target_file << "sequence\t";
  search_target_file << "protein id(loc) 1\t";
  search_target_file << "protein id(loc) 2\t";
  search_target_file << "by total\t";
  search_target_file << "by observable (0-1200)\t";
  search_target_file << "by observable bin (0-1200)\t";
  search_target_file << "by observable (0-max)\t";
  search_target_file << "by observable bin (0-max)\t";
  search_target_file << "by observed bin\t";
  search_target_file << "ion current total\t";
  search_target_file << "ion current observed"<<"\t";
  search_target_file << "ions observable bin (0-1200)"<<endl;

  const char *decoy_filename = "search.decoy.txt";
  string decoy_path = string(output_directory) + "/" + string(decoy_filename);

  ofstream search_decoy_file (decoy_path.c_str());
  //set precision
  search_decoy_file << setprecision(get_int_parameter("precision"));
  //print header
  search_decoy_file << "scan\t";
  search_decoy_file << "charge\t";
  search_decoy_file << "spectrum precursor m/z\t";
  search_decoy_file << "spectrum neutral mass\t";
  search_decoy_file << "peptide mass mono\t";
  search_decoy_file << "peptide mass average\t";
  search_decoy_file << "mass error(ppm)\t";
  search_decoy_file << "xcorr score\t";
  search_decoy_file << "xcorr rank\t";
  search_decoy_file << "p-value\t";
  search_decoy_file << "matches/spectrum\t";
  search_decoy_file << "sequence"<<endl;

  Scorer hhc_scorer;
  // main loop over spectra in ms2 file
 
  int search_count = 0;

  // for every observed spectrum 
  while (spectrum_iterator->hasNext()) {
    int charge;
    spectrum = spectrum_iterator->next(zstate);

    charge = zstate.getCharge();

    //SCORER_T* scorer = new_scorer(XCORR);
    scan_num = spectrum->getFirstScan();

    if (search_count % 100 == 0)
      carp(CARP_INFO,"count %d scan %d charge %d", search_count, scan_num, charge);
    search_count++;

    //vector<pair<FLOAT_T, LinkedPeptide> > linked_scores;
    //vector<pair<FLOAT_T, LinkedPeptide> > single_scores;
    vector<pair<FLOAT_T, LinkedPeptide> > scores;

    vector<LinkedPeptide> target_xpeptides;
    vector<LinkedPeptide> target_decoy_xpeptides;
    vector<LinkedPeptide> decoy_train_xpeptides;
    vector<LinkedPeptide> decoy_xpeptides;

    FLOAT_T precursor_mz = spectrum->getPrecursorMz();
    FLOAT_T precursor_mass = zstate.getNeutralMass(); 
 


    clock_t start_clock = clock();

    
    carp(CARP_DEBUG, "finding target xpeptides in mass window...%g", precursor_window);
    get_ions_from_window(
      target_xpeptides,
      all_ions,
      precursor_mass,
      precursor_window,
      precursor_window_type
      );

    if (target_xpeptides.size() < 1) {
      carp(CARP_INFO, "not enough precursors found in range, skipping scan %d charge %d", scan_num, charge);
      continue;
    }
    

    carp(CARP_DEBUG, "finding training xpeptides in decoy precursor window..%g", precursor_window_decoy);
    get_ions_from_window(
	target_decoy_xpeptides,
	all_ions,
	precursor_mass,
	precursor_window_decoy,
	window_type_decoy);
    
    carp(CARP_DETAILED_DEBUG, "Creating decoys for target window");
    //create the decoys from the target found in the target_mass_window.
    for (vector<LinkedPeptide>::iterator ion = target_xpeptides.begin();
	 ion != target_xpeptides.end(); ++ion) {
        add_decoys(decoy_xpeptides, *ion);
    }
    
    
    carp(CARP_DETAILED_DEBUG, "Creating decoys for decoy mass window");
    //create the decoys from the target found in the decoy_mass_window.
    while (decoy_train_xpeptides.size() < min_weibull_points) {
      for (vector<LinkedPeptide>::iterator ion = target_decoy_xpeptides.begin();
	   ion != target_decoy_xpeptides.end(); ++ion) {
	add_decoys(decoy_train_xpeptides, *ion);
      }
    }    

    carp(CARP_DEBUG, "num targets:%d",target_xpeptides.size());
    carp(CARP_DEBUG, "num decoys:%d", decoy_xpeptides.size());
    carp(CARP_DEBUG, "num training decoys:%d", decoy_train_xpeptides.size());

    clock_t candidate_clock = clock();

    int max_charge = min(max_ion_charge, charge);

    LinkedIonSeries ion_series = LinkedIonSeries(max_charge);

    // for every ion in the mass window
    carp(CARP_DEBUG, "Scoring targets");
    for (unsigned int idx=0;idx<target_xpeptides.size();idx++) {
      //LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.clear();
      ion_series.add_linked_ions(target_xpeptides[idx]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      scores.push_back(make_pair(score, target_xpeptides[idx]));
    }

    clock_t target_clock = clock();

    carp(CARP_DEBUG, "Scoring decoys.");
    for (unsigned int idx=0;idx<decoy_xpeptides.size();idx++) {
      //LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.clear();
      ion_series.add_linked_ions(decoy_xpeptides[idx]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      scores.push_back(make_pair(score, decoy_xpeptides[idx]));
    }


    //use the decoy scores to build the estimator.
    // create arrays to pass to crux's weibull methods
    FLOAT_T* linked_decoy_scores_array = new FLOAT_T[decoy_train_xpeptides.size()+target_xpeptides.size()];

    clock_t decoy_clock = clock();
    carp(CARP_DEBUG, "scoring training decoys...");
    // score all training decoys
    for (unsigned int idx=0;idx<decoy_train_xpeptides.size();idx++) {
      //LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.clear();
      ion_series.add_linked_ions(decoy_train_xpeptides[idx]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      linked_decoy_scores_array[idx] = score;
    }
  
    

    clock_t train_decoy_clock = clock();




    
    for (unsigned int idx=0;idx<scores.size();idx++) {
      if (!scores[idx].second.is_decoy())
	linked_decoy_scores_array[idx+decoy_train_xpeptides.size()] = scores[idx].first;
    }
    
    // sort scores
    sort(scores.begin(), scores.end(), greater<pair<FLOAT_T, LinkedPeptide> >());

    clock_t create_array_clock = clock();


   // weibull parameters for candidates
    FLOAT_T eta_linked = 0.0;
    FLOAT_T beta_linked  = 0.0;
    FLOAT_T shift_linked  = 0.0;
    FLOAT_T correlation_linked  = 0.0;

    // fit weibull to decoys

    hhc_estimate_weibull_parameters_from_xcorrs(linked_decoy_scores_array, 
						decoy_train_xpeptides.size(), 
						&eta_linked, &beta_linked, 
						&shift_linked, &correlation_linked, 
						spectrum, charge);
    
    clock_t weibull_clock = clock();
    double candidate_time = (double(candidate_clock) - double(start_clock)) / CLOCKS_PER_SEC;
    double target_time = (double(target_clock) - double(candidate_clock)) / CLOCKS_PER_SEC;
    double decoy_time = (double(decoy_clock) - double(target_clock)) / CLOCKS_PER_SEC;
    double train_decoy_time = (double(train_decoy_clock) - double(decoy_clock)) / CLOCKS_PER_SEC;
    double create_array_time =(double(create_array_clock) - double(train_decoy_clock)) / CLOCKS_PER_SEC;
    double weibull_time = (double(weibull_clock) - double(create_array_clock)) / CLOCKS_PER_SEC;

    carp(CARP_DEBUG, "candidate:%g", candidate_time);
    carp(CARP_DEBUG, "target:%g", target_time);
    carp(CARP_DEBUG, "decoy:%g", decoy_time);
    carp(CARP_DEBUG, "train decoy:%g", train_decoy_time);
    carp(CARP_DEBUG, "create array:%g", create_array_time);
    carp(CARP_DEBUG, "weibull:%g", weibull_time);
    carp(CARP_DEBUG, "========================");


    int ndecoys = 0;
    int ntargets = 0;
    unsigned int score_index = 0;

    while (score_index < scores.size() && (ndecoys < top_match || ntargets < top_match)) {
 
      
      double ppm_error = fabs(scores[score_index].second.mass(MONO) - precursor_mass) / 
          scores[score_index].second.mass(MONO) * 1e6;

      double pvalue = compute_weibull_pvalue(scores[score_index].first, eta_linked, beta_linked, shift_linked);
      double pvalue_bonf = pvalue;//bonf_correct(pvalue, decoy_xpeptides.size());
	
      if (pvalue != pvalue) {
        pvalue = 1;
	pvalue_bonf = 1;
      } else if (pvalue_bonf != pvalue_bonf) {
        pvalue_bonf = 1;
      }

      if (scores[score_index].second.is_decoy() && ndecoys < top_match) {
        ndecoys++;
	search_decoy_file << scan_num << "\t"; 
	search_decoy_file << charge << "\t"; 
	search_decoy_file << precursor_mz << "\t";
	search_decoy_file << precursor_mass << "\t";
	search_decoy_file << scores[score_index].second.mass(MONO) << "\t";
	search_decoy_file << scores[score_index].second.mass(AVERAGE) << "\t";
        search_decoy_file << ppm_error << "\t";
	search_decoy_file << scores[score_index].first <<"\t";
	search_decoy_file << ndecoys << "\t";
	search_decoy_file << pvalue << "\t";
	search_decoy_file << decoy_xpeptides.size() << "\t";
	search_decoy_file << scores[score_index].second<<endl;

      } else if (!scores[score_index].second.is_decoy() && ntargets < top_match) {
	ntargets++;
	search_target_file << scan_num << "\t"; 
	search_target_file << charge << "\t"; 
	search_target_file << precursor_mz << "\t";
	search_target_file << precursor_mass << "\t";
	search_target_file << scores[score_index].second.mass(MONO) << "\t";
	search_target_file << scores[score_index].second.mass(AVERAGE) << "\t";
        search_target_file << ppm_error << "\t";
	search_target_file << scores[score_index].first <<"\t";
	search_target_file << ntargets << "\t";
	search_target_file << pvalue << "\t";
	search_target_file << target_xpeptides.size() << "\t";
	search_target_file << scores[score_index].second<<"\t";

        //output protein ids/peptide locations.  If it is a linear, dead or self loop, only
        //use the 1st field.
        string sequence1  = scores[score_index].second.peptides()[0].sequence();
        vector<PEPTIDE_T*>& peptides1 = get_peptides_from_sequence(sequence1);
        string result_string = get_protein_ids_locations(peptides1);
        search_target_file << result_string << "\t";
        //if it is cross-linked peptide, use the second field
        if (scores[score_index].second.is_linked()) {
          string sequence2  = scores[score_index].second.peptides()[1].sequence();
          vector<PEPTIDE_T*>& peptides2 = get_peptides_from_sequence(sequence2);
          string result_string = get_protein_ids_locations(peptides2);
          search_target_file << result_string;
        }
        search_target_file <<"\t";
                //get theoretical ions count for (0-1200, with 1Da bins).
        Scorer scorer;
        LinkedIonSeries ion_series(charge);
        ion_series.add_linked_ions(scores[score_index].second);

        FLOAT_T ion_current_observed;
        FLOAT_T ion_current_total = spectrum->getTotalEnergy();
        int by_total = ion_series.get_total_by_ions();
        int by_observable;
        int by_observable2;
        int by_observable_bin;
        int by_observable_bin2;
        int by_observed_bin;
        int ions_observable;
        int ions_observable_bin;
        ion_series.get_observable_ions(0, 1200, bin_width_mono, ions_observable, ions_observable_bin);
        ion_series.get_observable_by_ions(0, 1200, bin_width_mono, by_observable, by_observable_bin);
        ion_series.get_observable_by_ions(0, spectrum->getMaxPeakMz(), bin_width_mono, by_observable2, by_observable_bin2);
        scorer.getIonCurrentExplained(ion_series, spectrum, ion_current_observed, by_observed_bin);
        
  

        search_target_file << by_total << "\t";
        search_target_file << by_observable << "\t";
        search_target_file << by_observable_bin << "\t";
        search_target_file << by_observable2 << "\t";
        search_target_file << by_observable_bin2 << "\t";
        search_target_file << by_observed_bin << "\t";
        search_target_file << ion_current_total << "\t";
        search_target_file << ion_current_observed << "\t";
        search_target_file << ions_observable_bin;
        
        search_target_file << endl;


        } 
      score_index++;
    }

    delete [] linked_decoy_scores_array;
    //free_spectrum(spectrum);

    carp(CARP_DETAILED_DEBUG,"Done with spectrum %d", scan_num);
  } // get next spectrum
  search_target_file.close();
  search_decoy_file.close();
  //free_spectrum_collection(spectra);
  //free_spectrum(spectrum);

  //Calculate q-values.
  carp(CARP_INFO,"Computing Q-Values");
  xlink_compute_qvalues();

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux search-for-xlinks.");

  return(0);
}


void get_ions_from_window(vector<LinkedPeptide>& filtered_ions,
  vector<LinkedPeptide>& all_ions,
			  FLOAT_T precursor_mass,
			  FLOAT_T window,
			  WINDOW_TYPE_T window_type) {

  double min_mass = 0;
  double max_mass = 0;
  
  if (window_type == WINDOW_MASS) {
    min_mass = precursor_mass - window;
    max_mass = precursor_mass + window;
  } else if (window_type == WINDOW_PPM) {
    min_mass = precursor_mass / (1.0 + window * 1e-6);
    max_mass = precursor_mass / (1.0 - window * 1e-6);
  } else {
    carp(CARP_FATAL,"Precursor m/z window type not supported!");
  }

  get_ions_from_mass_range(filtered_ions, all_ions, min_mass, max_mass);

}


void get_ions_from_mass_range(
  vector<LinkedPeptide>& filtered_ions,
  vector<LinkedPeptide>& all_ions,
  double min_mass,
  double max_mass) {

  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");

  filtered_ions.clear();
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin();
    ion != all_ions.end();
    ++ion) {
    double mass = ion -> mass(mass_type);
    if (mass >= min_mass && mass <= max_mass) {
      filtered_ions.push_back(*ion);
    }
  }

}

// get all precursor ions within given mass window
void get_ions_from_mz_range(vector<LinkedPeptide>& filtered_ions,
	vector<LinkedPeptide>& all_ions,
	FLOAT_T precursor_mass,
	int charge,
	FLOAT_T mass_window,
	int decoy_iterations) {
  FLOAT_T min_mass = precursor_mass - mass_window;
  FLOAT_T max_mass = precursor_mass + mass_window;
  carp(CARP_DETAILED_DEBUG,"get_ions_from_mz_range()");
  carp(CARP_DETAILED_DEBUG,"min_mass %g max_mass %g", min_mass, max_mass);

  FLOAT_T ion_mass;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin();
	ion != all_ions.end(); ++ion) {
    ion->set_charge(charge);
    ion->calculate_mass(get_mass_type_parameter("isotopic-mass"));
    ion_mass = ion->mass(get_mass_type_parameter("isotopic-mass"));
    if (ion_mass >= min_mass && ion_mass <= max_mass) {
      filtered_ions.push_back(*ion);
      for (int i = decoy_iterations; i > 0; --i)
        add_decoys(filtered_ions, *ion);
    }
  }
}


#define BONF_CUTOFF_P 1e-4
#define BONF_CUTOFF_NP 1e-2

double bonf_correct(double nlp_value, int n) {
  if (nlp_value != nlp_value) return 0;
  if (nlp_value == 0) return 0;

  double NL_BONF_CUTOFF_P = (-log(BONF_CUTOFF_P));
  double NL_BONF_CUTOFF_NP= (-log(BONF_CUTOFF_NP));


  double ans = nlp_value - log((double)n);
 
  if ((nlp_value <= NL_BONF_CUTOFF_P) || 
      (ans <= NL_BONF_CUTOFF_NP)) { 
    double p = exp(-nlp_value);
    ans = -log(1-pow((1-p), n));
  }
  return ans;
}


void get_protein_ids_locations(PEPTIDE_T *peptide, 
  set<string>& protein_ids_locations) {

  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator = 
    new_peptide_src_iterator(peptide);

  std::ostringstream protein_field_stream;

  if (peptide_src_iterator_has_next(peptide_src_iterator)) {
    while(peptide_src_iterator_has_next(peptide_src_iterator)){
      PEPTIDE_SRC_T* peptide_src = peptide_src_iterator_next(peptide_src_iterator);
      Protein* protein = get_peptide_src_parent_protein(peptide_src);
      char* protein_id = protein->getId();
      int peptide_loc = get_peptide_src_start_idx(peptide_src);
      std::ostringstream protein_loc_stream;
      protein_loc_stream << protein_id << "(" << peptide_loc << ")";
      free(protein_id);
      protein_ids_locations.insert(protein_loc_stream.str());
    }
  }
  free(peptide_src_iterator);
}

string get_protein_ids_locations(vector<PEPTIDE_T*>& peptides) {
  set<string> protein_ids_locations;

  for (unsigned int idx=0;idx<peptides.size();idx++) {
    get_protein_ids_locations(peptides[idx], protein_ids_locations);
  }

  set<string>::iterator result_iter = protein_ids_locations.begin();

  string protein_field_string = *result_iter;

  while(++result_iter != protein_ids_locations.end()) {
    protein_field_string += "," + *result_iter;
  }

  return protein_field_string;

}


