/**
 * \file xhhc_search.cpp
 * AUTHOR: Sean McIlwain
 * \brief Main file for crux-search-for-xlinks.
 *
 * Given an ms2 file and a fasta file or index, compare all spectra to
 * peptides and cross-linked products in the fasta file/index and 
 * return high scoring matches.
 * Products are determined by parameters for length, mass, mass
 * tolerance, cleavages. Score each spectrum with
 * respect to all candidates, and rank by score. 
 * tab-delimited file format.
 */

//XLINK INCLUDES
#include "xhhc_scorer.h"
#include "LinkedIonSeries.h"
#include "xlink_compute_qvalues.h"
#include "SearchForXLinks.h"
#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

//CRUX INCLUDES
#include "objects.h"
#include "Spectrum.h"
#include "SpectrumCollectionFactory.h"
#include "FilteredSpectrumChargeIterator.h"

//C++ includes
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace Crux;

//PRIVATE FUNCTIONS

/**
 * populate the filtered_ions vector with the products that
 * fit within the precursor mass window
 */
void get_ions_from_window(
  vector<LinkedPeptide>& filtered_ions, ///< filtered ions -out
  vector<LinkedPeptide>& all_ions, ///< all ions in the database -in
  FLOAT_T precursor_mass, ///<precursor mass -in
  FLOAT_T window, ///< window value -in
  WINDOW_TYPE_T window_type ///< window type -in
);

/**
 * populated the filtered_ions vector with the products that
 * have mass that is within the (min_mass,max_mass).
 */
void get_ions_from_mass_range(
  vector<LinkedPeptide>& filtered_ions, ///< filtered ions -out
  vector<LinkedPeptide>& all_ions, ///< all ions in the database -in
  double min_mass, ///< min mass of ions to select -in
  double max_mass ///< max mass of ions to select -in
);

/**
 * \returns a comma delimited string of protein_id(peptide_start)
 * locations for the peptides
 */ 
string get_protein_ids_locations(
  vector<Peptide*>& peptides ///< vector of peptides
  );

/**
 * main method for search-for-xlinks using the original code for the paper.
 */
int SearchForXLinks::xhhcSearchMain() {

  carp(CARP_INFO, "Beginning crux xlink-search (original)");

  //Get parameters
  const char* ms2_file = get_string_parameter_pointer("ms2 file");

  FLOAT_T precursor_window = get_double_parameter("precursor-window");
  FLOAT_T precursor_window_weibull = get_double_parameter("precursor-window-weibull");
  WINDOW_TYPE_T precursor_window_type = 
    get_window_type_parameter("precursor-window-type");
  WINDOW_TYPE_T window_type_weibull = 
    get_window_type_parameter("precursor-window-type-weibull");

  const char* output_directory = get_string_parameter_pointer("output-dir");

  unsigned int min_weibull_points = 
    (unsigned int)get_int_parameter("min-weibull-points");
  int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");
  int top_match = get_int_parameter("top-match");
  FLOAT_T linker_mass = get_double_parameter("link mass");
  int scan_num = 0;
  SpectrumZState zstate;
 
  //Set the static variable for the linker mass
  LinkedPeptide::setLinkerMass(linker_mass);

  vector<LinkedPeptide> all_ions;
  carp(CARP_DETAILED_DEBUG,"Calling find all precursor ions");
  find_all_precursor_ions(all_ions);
  carp(CARP_DETAILED_DEBUG,"Sort");
  // sort filtered ions and decoy ions by mass
  //sort(all_ions.begin(), all_ions.end());
  if (get_boolean_parameter("xlink-print-db")) {
    ostringstream oss;
    oss << output_directory << "/" << "xlink_peptides.txt";
    string temp = oss.str();
    ofstream peptides_file(temp.c_str());
    peptides_file << "mass\tsequence"<<endl;
    for (unsigned int idx=0;idx < all_ions.size();idx++) {
      peptides_file << all_ions.at(idx).getMass(MONO) << "\t";
      peptides_file << all_ions.at(idx) << endl;
    }
    peptides_file.flush();
    carp(CARP_INFO, "outputted database to xlink_peptides.txt");
    return 0;
  }

  carp(CARP_INFO, "Loading Spectra");
  Spectrum* spectrum = NULL;
  Crux::SpectrumCollection* spectra = SpectrumCollectionFactory::create(ms2_file);
  spectra->parse();

  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);
 
  FLOAT_T score;
 // best pvalues

  const char* target_filename = "search-for-xlinks.target.txt";
  
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

  const char *decoy_filename = "search-for-xlinks.decoy.txt";
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

  XHHC_Scorer hhc_scorer;
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
    

    carp(CARP_DEBUG, "finding training xpeptides in decoy precursor window..%g", precursor_window_weibull);
    get_ions_from_window(
  target_decoy_xpeptides,
  all_ions,
  precursor_mass,
  precursor_window_weibull,
  window_type_weibull);
    
    carp(CARP_DETAILED_DEBUG, "Creating decoys for target window");
    //create the decoys from the target found in the target_mass_window.
    for (vector<LinkedPeptide>::iterator ion = target_xpeptides.begin();
   ion != target_xpeptides.end(); ++ion) {
        add_decoy(decoy_xpeptides, *ion);
    }
    
    
    carp(CARP_DETAILED_DEBUG, "Creating decoys for decoy mass window");
    //create the decoys from the target found in the decoy_mass_window.
    while ((decoy_train_xpeptides.size() + target_xpeptides.size()) < min_weibull_points) {
      for (vector<LinkedPeptide>::iterator ion = target_decoy_xpeptides.begin();
     ion != target_decoy_xpeptides.end(); ++ion) {
  add_decoy(decoy_train_xpeptides, *ion);
      }
    }    

    size_t num_training_points = decoy_train_xpeptides.size() + target_xpeptides.size();

    carp(CARP_DEBUG, "num targets:%d",target_xpeptides.size());
    carp(CARP_DEBUG, "num decoys:%d", decoy_xpeptides.size());
    carp(CARP_DEBUG, "num training points:%d", num_training_points);

    clock_t candidate_clock = clock();

    int max_charge = min(max_ion_charge, charge);

    LinkedIonSeries ion_series = LinkedIonSeries(max_charge);

    // for every ion in the mass window
    carp(CARP_DEBUG, "Scoring targets");
    for (unsigned int idx=0;idx<target_xpeptides.size();idx++) {
      //LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.clear();
      ion_series.addLinkedIons(target_xpeptides[idx]);
      score = hhc_scorer.scoreSpectrumVsSeries(spectrum, ion_series);
      scores.push_back(make_pair(score, target_xpeptides[idx]));
    }

    carp(CARP_DEBUG, "Scoring decoys.");
    for (unsigned int idx=0;idx<decoy_xpeptides.size();idx++) {
      //LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.clear();
      ion_series.addLinkedIons(decoy_xpeptides[idx]);
      score = hhc_scorer.scoreSpectrumVsSeries(spectrum, ion_series);
      scores.push_back(make_pair(score, decoy_xpeptides[idx]));
    }


    //use the decoy scores to build the estimator.
    // create arrays to pass to crux's weibull methods
    FLOAT_T* linked_decoy_scores_array = new FLOAT_T[num_training_points];

    clock_t decoy_clock = clock();
    carp(CARP_DEBUG, "scoring training decoys...");
    // score all training decoys
    for (unsigned int idx=0;idx<decoy_train_xpeptides.size();idx++) {
      //LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.clear();
      ion_series.addLinkedIons(decoy_train_xpeptides[idx]);
      score = hhc_scorer.scoreSpectrumVsSeries(spectrum, ion_series);
      linked_decoy_scores_array[idx] = score;
    }
  
    clock_t train_decoy_clock = clock();
  
    for (unsigned int idx=0;idx<scores.size();idx++) {
      if (!scores[idx].second.isDecoy())
  linked_decoy_scores_array[idx+decoy_train_xpeptides.size()] = scores[idx].first;
    }


   // weibull parameters for candidates
    FLOAT_T eta_linked = 0.0;
    FLOAT_T beta_linked  = 0.0;
    FLOAT_T shift_linked  = 0.0;
    FLOAT_T correlation_linked  = 0.0;

    // fit weibull to decoys

    carp(CARP_DEBUG, "Fitting weibull to %d scores", num_training_points);
    hhc_estimate_weibull_parameters_from_xcorrs(linked_decoy_scores_array, 
            num_training_points, 
            &eta_linked, &beta_linked, 
            &shift_linked, &correlation_linked, 
            spectrum, charge);

    //okay we don't need the training peptides/scores anymore
    delete []linked_decoy_scores_array;
    vector<LinkedPeptide>().swap(decoy_train_xpeptides);
        
    // sort scores
    
    carp(CARP_DEBUG, "sorting %u scores", scores.size());
    sort(scores.begin(), scores.end(), greater<pair<FLOAT_T, LinkedPeptide> >());
    carp(CARP_INFO, "done sorting");
    int ndecoys = 0;
    int ntargets = 0;
    unsigned int score_index = 0;

    while (score_index < scores.size() && (ndecoys < top_match || ntargets < top_match)) {
 
      
      double ppm_error = fabs(scores[score_index].second.getMass(MONO) - precursor_mass) / 
          scores[score_index].second.getMass(MONO) * 1e6;

      double pvalue = compute_weibull_pvalue(scores[score_index].first, eta_linked, beta_linked, shift_linked);
  
      if (pvalue != pvalue) {
        pvalue = 1;
      }

      if (scores[score_index].second.isDecoy() && ndecoys < top_match) {
        ndecoys++;
        search_decoy_file << scan_num << "\t"; 
        search_decoy_file << charge << "\t"; 
        search_decoy_file << precursor_mz << "\t";
        search_decoy_file << precursor_mass << "\t";
        search_decoy_file << scores[score_index].second.getMass(MONO) << "\t";
        search_decoy_file << scores[score_index].second.getMass(AVERAGE) << "\t";
        search_decoy_file << ppm_error << "\t";
        search_decoy_file << scores[score_index].first <<"\t";
        search_decoy_file << ndecoys << "\t";
        search_decoy_file << pvalue << "\t";
        search_decoy_file << decoy_xpeptides.size() << "\t";
        search_decoy_file << scores[score_index].second<<endl;

      } else if (!scores[score_index].second.isDecoy() && ntargets < top_match) {
        ntargets++;
        search_target_file << scan_num << "\t"; 
        search_target_file << charge << "\t"; 
        search_target_file << precursor_mz << "\t";
        search_target_file << precursor_mass << "\t";
        search_target_file << scores[score_index].second.getMass(MONO) << "\t";
        search_target_file << scores[score_index].second.getMass(AVERAGE) << "\t";
        search_target_file << ppm_error << "\t";
        search_target_file << scores[score_index].first <<"\t";
        search_target_file << ntargets << "\t";
        search_target_file << pvalue << "\t";
        search_target_file << target_xpeptides.size() << "\t";
        search_target_file << scores[score_index].second<<"\t";

        //output protein ids/peptide locations.  If it is a linear, dead or self loop, only
        //use the 1st field.
        string sequence1  = scores[score_index].second.getPeptides()[0].getSequence();
        vector<Peptide*>& peptides1 = get_peptides_from_sequence(sequence1);
        string result_string = get_protein_ids_locations(peptides1);
        search_target_file << result_string << "\t";
        //if it is cross-linked peptide, use the second field
        if (scores[score_index].second.isCrossLinked()) {
          string sequence2  = scores[score_index].second.getPeptides()[1].getSequence();
          vector<Peptide*>& peptides2 = get_peptides_from_sequence(sequence2);
          string result_string = get_protein_ids_locations(peptides2);
          search_target_file << result_string;
        }
        search_target_file <<"\t";
                //get theoretical ions count for (0-1200, with 1Da bins).
        XHHC_Scorer scorer;
        LinkedIonSeries ion_series(charge);
        ion_series.addLinkedIons(scores[score_index].second);

        FLOAT_T ion_current_observed;
        FLOAT_T ion_current_total = spectrum->getTotalEnergy();
        int by_total = ion_series.getTotalBYIons();
        int by_observable;
        int by_observable2;
        int by_observable_bin;
        int by_observable_bin2;
        int by_observed_bin;
        int ions_observable;
        int ions_observable_bin;
        ion_series.getObservableIons(0, 1200, bin_width_mono, ions_observable, ions_observable_bin);
        ion_series.getObservableBYIons(0, 1200, bin_width_mono, by_observable, by_observable_bin);
        ion_series.getObservableBYIons(0, spectrum->getMaxPeakMz(), bin_width_mono, by_observable2, by_observable_bin2);
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

    //free_spectrum(spectrum);

    carp(CARP_DETAILED_DEBUG,"Done with spectrum %d", scan_num);
  } // get next spectrum
  search_target_file.close();
  search_decoy_file.close();

  //Calculate q-values.
  carp(CARP_INFO,"Computing Q-Values");

  //free memory.
  xlink_compute_qvalues();
  delete spectrum_iterator;
  delete spectra;
  free_peptides();

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );
  for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    free_peptide_mod(peptide_mods[mod_idx]);
  }
  free(peptide_mods);

  free_parameters();
  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux search-for-xlinks.");

  return(0);
}

/**
 * populate the filtered_ions vector with the ions that
 * fit within the precursor mass window
 */
void get_ions_from_window(
  vector<LinkedPeptide>& filtered_ions, ///< filtered ions -out
  vector<LinkedPeptide>& all_ions, ///< all ions in the database -in
  FLOAT_T precursor_mass, ///<precursor mass -in
  FLOAT_T window, ///< window value -in
  WINDOW_TYPE_T window_type ///< window type -in
  ) {

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

/**
 * populated the filtered_ions vector with the products that
 * have mass that is within the (min_mass,max_mass).
 */
void get_ions_from_mass_range(
  vector<LinkedPeptide>& filtered_ions, ///< filtered ions -out
  vector<LinkedPeptide>& all_ions, ///< all ions in the database -in
  double min_mass, ///< min mass of ions to select -in
  double max_mass ///< max mass of ions to select -in
  ) {

  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");

  filtered_ions.clear();
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin();
    ion != all_ions.end();
    ++ion) {
    double mass = ion -> getMass(mass_type);
    if (mass >= min_mass && mass <= max_mass) {
      filtered_ions.push_back(*ion);
    }
  }
}

/**
 * given a peptide, populate the protein_ids_locations
 * with a set of strings that are protein_id(peptide start index).
 */
void get_protein_ids_locations(
  Peptide *peptide, ///< peptide to generate locations from -in
  set<string>& protein_ids_locations ///< set of protein_id(peptide start index). -out
  ) {

  std::ostringstream protein_field_stream;

  for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
       iter != peptide->getPeptideSrcEnd();
       ++iter) {

    PeptideSrc* peptide_src = *iter;
    Protein* protein = peptide_src->getParentProtein();
    char* protein_id = protein->getIdPointer();
    int peptide_loc = peptide_src->getStartIdx();
    std::ostringstream protein_loc_stream;
    protein_loc_stream << protein_id << "(" << peptide_loc << ")";
    protein_ids_locations.insert(protein_loc_stream.str());
    
  }

}

/**
 * \returns a comma delimited string of protein_id(peptide_start)
 * locations for the peptides
 */ 
string get_protein_ids_locations(
  vector<Peptide*>& peptides ///< vector of peptides
  ) {
  
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


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
