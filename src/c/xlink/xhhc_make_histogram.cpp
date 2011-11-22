#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"

#include "objects.h"
#include "IonConstraint.h"
#include "SpectrumCollectionFactory.h"

#include <fstream>
#include <math.h>
#include <iostream>

//#define PARAM_ESTIMATION_SAMPLE_COUNT 500
#define MIN_WEIBULL_MATCHES 40
#define MIN_XCORR_SHIFT -5.0
#define MAX_XCORR_SHIFT  5.0
#define XCORR_SHIFT 0.05
using namespace std;

//typedef map<char, set<char> > BondMap;

void plot_weibull(vector<pair<FLOAT_T, LinkedPeptide> >& scores, 
                  Spectrum* spectrum, 
                  int charge); 

int main(int argc, char** argv) {
  const char* missed_link_cleavage = "K";
  //int num_missed_cleavages = 0;
  char* ms2_file = NULL;
  char* min_mass_string = NULL;
  char* max_mass_string = NULL;
  char* database = NULL;
  char* links = NULL;
  char* linker_mass_string = NULL;
  int decoy_iterations = 5;
  int charge = 1;
  int scan_num = 0;
  bool open_modification = false;
  int open_modification_int = 0;
  parse_arguments_set_req(
    "protein database", 
    "database containing all proteins", 
    (void *) &database, 
    STRING_ARG);

  parse_arguments_set_req(
    "links", 
    "comma delimited pair of amino acid link sites, ex. A:K,A:D", 
    (void *) &links, 
    STRING_ARG);
/*
  parse_arguments_set_req(
    "max charge", 
    "maximum charge for ions", 
    (void *) &charge, 
    INT_ARG);
*/
  parse_arguments_set_req(
    "linker mass", 
    "combined mass of linker and linker modifications", 
    (void *) &linker_mass_string, 
    STRING_ARG);

  parse_arguments_set_req(
    "scan-number", 
    "The scan number for the MS-MS spectrum to extract from the ms2 file. This is an integer in the range [1, 100000], and uniquely identifies a particular MS-MS spectrum within an .ms2 file.",
    (void *) &scan_num, INT_ARG);

  parse_arguments_set_req(
    "ms2-filename", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);
 
  parse_arguments_set_opt(
    "open-modification",
    "",
    (void *) &open_modification_int, 
    INT_ARG);

  parse_arguments_set_opt(
    "decoy-iterations",
    "",
    (void *) &decoy_iterations,
    INT_ARG);

  // not implemented yet
  /*TODO  
  parse_arguments_set_opt(
    "missed-link-cleavage",
    "",
    (void *) &missed_link_cleavage, 
    STRING_ARG);
  */

  // not implemented yet
  /*
  parse_arguments_set_opt(
    "num-missed-cleavages", 
    "maximum number of missed cleavages (not including one at link site)", 
    (void *) &num_missed_cleavages, 
    INT_ARG);
  */
  parse_arguments_set_opt(
    "charge",
    "peptide charge", 
    (void *) &charge,
    INT_ARG); 

  parse_arguments_set_opt(
    "min-mass", 
    "", 
    (void *) &min_mass_string, 
    STRING_ARG);

  parse_arguments_set_opt(
    "max-mass", 
    "", 
    (void *) &max_mass_string, 
    STRING_ARG);

  initialize_parameters();
  if (!parse_arguments(argc, argv, 0)) {
   char* error_message;
   char* usage = parse_arguments_get_usage("xhhc-make-histogram");
   int result = parse_arguments_get_error(&error_message);
   fprintf(stderr, "Error in command line. Error # %d\n", result);
   fprintf(stderr, "%s\n", error_message);
   fprintf(stderr, "%s", usage);
   free(usage);
   exit(1);
 }
  // something wrong with DOUBLE_ARG
  FLOAT_T linker_mass = atof(linker_mass_string);
  // and boolean arg
  if (open_modification_int == 1)
    open_modification = true;
  cout << "ms2 " << ms2_file << " charge " << charge << " scan num " << scan_num << endl;

  vector<LinkedPeptide> all_ions;
  
  //find_all_precursor_ions(all_ions, links, linker_mass, charge, missed_link_cleavage, database);
  find_all_precursor_ions(all_ions, (char*)links, (char*)missed_link_cleavage, database, charge);

  FLOAT_T max_mass = all_ions.back().mass(AVERAGE);
  FLOAT_T min_mass = 0.0;
  if (min_mass_string != NULL) min_mass = atof(min_mass_string);
  if (max_mass_string != NULL) max_mass = atof(max_mass_string);
  if (max_mass < min_mass) {
    carp(CARP_FATAL, "max mass must be larger than min mass");
  }

  cout << "min " << min_mass << " max " << max_mass << endl;

  int num_ions = 0;
  vector<LinkedPeptide> filtered_ions;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {

      ion->calculate_mass(AVERAGE);
    // if the mass is in the range
    if (min_mass <= ion->mass(AVERAGE) && ion->mass(AVERAGE) <= max_mass) {
      //ion->set_charge(charge);
      ++num_ions;
      filtered_ions.push_back(*ion);
      // print out the ion
      //cout << ion->mass() << "\t" << *ion << endl;
      // iterate to add shuffled decoys
      //int i = 5;
      //if (charge == 4) i = 7;
      for (int i = decoy_iterations; i > 0; --i) {
        add_decoys(filtered_ions, *ion);
        //add_decoys(filtered_ions, *ion, links, charge, linker_mass);
      }
    } 
 }     
  
  // sort filtered ions and decoy ions by mass
  cout << "sorting ...";
  sort(filtered_ions.begin(), filtered_ions.end());
  cout << "done" << endl;
  cout << "scan " << scan_num << " +" << charge << "<br>" << endl;
  cout << "precursors  " << num_ions << "<br>" << endl;
  cout << "decoys      " << filtered_ions.size() - num_ions << "<br>" << endl;
  cout << "total       " << filtered_ions.size() << "<br>" << endl;

  Spectrum* spectrum = new Spectrum();
  SpectrumCollection* collection = SpectrumCollectionFactory::create(ms2_file);

  XHHC_Scorer xhhc_scorer;
  if(!collection->getSpectrum(scan_num, spectrum)){
    carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
    delete collection;
    delete spectrum;
    exit(1);
  }
  
  FLOAT_T score = 0;
  // Pragya's open modification method
  if (open_modification) {
  FLOAT_T mod_mass;
    Scorer* scorer = new Scorer(XCORR);
    IonSeries* ion_series = NULL; 
    IonConstraint* ion_constraint = 
	IonConstraint::newIonConstraintSequestXcorr(charge);
    set<pair<FLOAT_T, string> > scores;
    stringstream ss;
    // for every precursor in the mass window
    for (vector<LinkedPeptide>::iterator ion = filtered_ions.begin(); ion != filtered_ions.end(); ++ion) {
      if (ion->size() == 2) {
	vector<XHHC_Peptide> peptides = ion->peptides();
	// score the first peptide with modification of second peptide	
        mod_mass = linker_mass + peptides[1].mass(MONO);	
	ion_series = new IonSeries((char*)peptides[0].sequence().c_str(), ion->charge(), ion_constraint);
	hhc_predict_ions(ion_series, mod_mass, peptides[0].link_site());
	score = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);
	//score = xhhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
        ss.str("");
	ss << peptides[0].sequence() << " mod " << peptides[1].sequence() << ", " << peptides[0].link_site();
        scores.insert(make_pair(score, ss.str()));	
	// score second peptide with modification of first peptide
        mod_mass = linker_mass + peptides[0].mass(MONO);	
	ion_series = new IonSeries((char*)peptides[1].sequence().c_str(), ion->charge(), ion_constraint);
	hhc_predict_ions(ion_series, mod_mass, peptides[1].link_site());
        //score = xhhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
	score = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);
	ss.str("");
	ss << peptides[1].sequence() << " mod " << peptides[0].sequence() << ", " << peptides[1].link_site();
        scores.insert(make_pair(score, ss.str()));
      }
    }
    //sort(scores.begin(), scores.end());
    int i = 1;
    for (set<pair<FLOAT_T, string> >::reverse_iterator score_pair = scores.rbegin();
		score_pair != scores.rend();
		++score_pair) {

	cout << i << "\t" << score_pair->first << "\t" << score_pair->second << endl;
	++i;
    }

  } else { // linked peptide method
    vector<pair<FLOAT_T, LinkedPeptide> > scores;
    //LinkedIonSeries ion_series;
  // for every ion in the mass window
    LinkedIonSeries ion_series = LinkedIonSeries(charge);
    for (vector<LinkedPeptide>::iterator ion = filtered_ions.begin(); ion != filtered_ions.end(); ++ion) {
      ion_series.clear();
      ion_series.add_linked_ions(*ion);
      score = xhhc_scorer.scoreSpectrumVsSeries(spectrum, ion_series);
     //score = hhc_score_spectrum_v_ion_series(scorer, spectrum, ion_series);
      scores.push_back(make_pair(score, *ion));
    }
    sort(scores.begin(), scores.end());
    FLOAT_T range = scores.back().first - scores.front().first;
    cout << "xcorr range " << range << "<br>" << endl;
    plot_weibull(scores, spectrum, charge);
    reverse(scores.begin(), scores.end());
    for (int i = 0; i < 20; ++i) {
	cout << i+1 << "\t" << scores[i].first << "\t" << scores[i].second << endl;
    }
  }
  //free_scorer(scorer);
  delete collection;
  delete spectrum;
  return 0;
}

// for running experiments. plots fit and pvalues
void plot_weibull(vector<pair<FLOAT_T, LinkedPeptide> >& scores, 
                  Spectrum* spectrum, int charge) {
  
  ofstream target_fit_file ("fit.target");
  ofstream decoy_fit_file ("fit.decoy");
  ofstream target_score_file ("scores.target");
  ofstream decoy_score_file ("scores.decoy");
  ofstream target_pvalue_file ("pvalues.target");
  ofstream decoy_pvalue_file ("pvalues.decoy");

  int num_scores = scores.size();
  int num_targets = 0;
  int num_decoys = 0;

  FLOAT_T* decoy_scores_array = new FLOAT_T[num_scores];
  FLOAT_T* target_scores_array = new FLOAT_T[num_scores];

  for (vector<pair<FLOAT_T, LinkedPeptide> >::iterator score_pair = scores.begin();
	score_pair != scores.end(); ++score_pair) {
    if (score_pair->second.is_decoy()) {
      decoy_score_file << score_pair->first << endl;
      decoy_scores_array[num_decoys++] = score_pair->first;
    } else {
      target_score_file << score_pair->first << endl;
      target_scores_array[num_targets++] = score_pair->first;
    }
  }
  FLOAT_T eta_target = 0.0;
  FLOAT_T beta_target = 0.0;
  FLOAT_T shift_target = 0.0;
  FLOAT_T correlation_target = 0.0;

  FLOAT_T eta_decoy = 0.0;
  FLOAT_T beta_decoy = 0.0;
  FLOAT_T shift_decoy = 0.0;
  FLOAT_T correlation_decoy = 0.0;

  FLOAT_T y;

  // plot fit for targets
  hhc_estimate_weibull_parameters_from_xcorrs(target_scores_array, 
                                              num_targets, &eta_target, 
	&beta_target, &shift_target, &correlation_target, spectrum, charge);
  for (FLOAT_T x = scores.front().first; x <= scores.back().first; x = x + 0.01) {
      y = (beta_target / eta_target) * pow(((x+shift_target)/eta_target), beta_target - 1) * exp(- pow((x+shift_target)/eta_target, beta_target));
      target_fit_file << x << "\t" << y << endl;
  }

  // plot fit for decoys 
  hhc_estimate_weibull_parameters_from_xcorrs(decoy_scores_array, num_decoys, &eta_decoy, 
	&beta_decoy, &shift_decoy, &correlation_decoy, spectrum, charge);
  for (FLOAT_T x = scores.front().first; x <= scores.back().first; x = x + 0.01) {
      y = (beta_decoy / eta_decoy) * pow(((x+shift_decoy)/eta_decoy), beta_decoy - 1) 
	* exp(- pow((x+shift_decoy)/eta_decoy, beta_decoy));
      decoy_fit_file << x << "\t" << y << endl;
  }

  delete target_scores_array;
  delete decoy_scores_array;

  cout << "target correlation " << correlation_target << " <br>" << endl;
  cout << "decoy correlation " << correlation_decoy << " <br>" << endl;

  FLOAT_T pvalue_target;
  FLOAT_T pvalue_decoy;

  for (vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator score_pair = scores.rbegin();
       score_pair != scores.rend(); ++score_pair) {
    pvalue_target = compute_weibull_pvalue(score_pair->first, eta_target, beta_target, shift_target);
    pvalue_decoy = compute_weibull_pvalue(score_pair->first, eta_decoy, beta_decoy, shift_decoy);

    if (score_pair->second.is_decoy()) {
      target_pvalue_file << pvalue_target << endl;
      decoy_pvalue_file << pvalue_decoy << endl;
    }
  }
} 
