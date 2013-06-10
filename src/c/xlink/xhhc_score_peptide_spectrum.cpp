//TODO - Change cout to carps

#include "xhhc.h"
#include "LinkedIonSeries.h"
#include "xhhc_scorer.h"
#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

#include "objects.h"
#include "IonConstraint.h"
#include "Scorer.h"
#include "SpectrumCollectionFactory.h"


#include <math.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <iostream>
#include <fstream>

using namespace Crux;

#define bin_width_mono 1.0005079


#define NUM_ARGUMENTS 8
#define NUM_OPTIONS 4


double get_concat_score(char* peptideA, char* peptideB, int link_site, 
                        int charge, Spectrum* spectrum);
void print_spectrum(Spectrum* spectrum, LinkedIonSeries& ion_series);
int main(int argc, char** argv){

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  int num_options = NUM_OPTIONS;
  const char* option_list[NUM_OPTIONS] = {
    "verbosity",
    "version",
    "use-flanking-peaks",
    "xlink-score-method"
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
  set_verbosity_level( CARP_ERROR );
  
  /* Set default values for parameters in parameter.c */
  initialize_parameters();

  /* Define optional and required command line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  parse_cmd_line_into_params_hash(argc, argv, "xlink-score-peptide-spectrum");

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

  LinkedPeptide::setLinkerMass(get_double_parameter("link mass"));
 
  // create new ion series
  
  // a single peptide linked to itself
  if (strcmp(peptideB, "NULL") == 0) {
    cout << "B is null" << endl; 
    peptideB = NULL;
  }

  // read ms2 file
  Crux::SpectrumCollection* collection = SpectrumCollectionFactory::create(ms2_file);

  // search for spectrum with correct scan number
  Spectrum* spectrum = collection->getSpectrum(scan_num);
  if( spectrum == NULL ){
    carp(CARP_ERROR, "Failed to find spectrum with scan_num: %d", scan_num);
    delete collection;
    exit(1);
  }
  
  //created linked peptide.
  LinkedPeptide lp = LinkedPeptide(peptideA, peptideB, posA, posB, charge);

  cout <<"LinkedPeptide:"<<lp<<" mass:"<<lp.getMass(MONO)<<endl;
  
  XHHC_Scorer xhhc_scorer;
  xhhc_scorer.setPrint(false);



  string scoremethod(get_string_parameter("xlink-score-method"));

  if (scoremethod=="composite") {

    LinkedIonSeries ion_series(charge);

    //cout << lp << endl;
    
    ion_series.addLinkedIons(lp);
       
    double score = xhhc_scorer.scoreSpectrumVsSeries(spectrum, ion_series);

    cout <<score<<endl;

    bool do_print_spectra = true;
    if (do_print_spectra) {
      print_spectrum(spectrum, ion_series);
    }
  } else if (scoremethod=="modification") {
    
    LinkedIonSeries ion_seriesA;
    ion_seriesA.addLinkedIons(lp, SPLITTYPE_A);
    double scoreA = xhhc_scorer.scoreSpectrumVsSeries(spectrum, ion_seriesA);
    
    LinkedIonSeries ion_seriesB;
    ion_seriesB.addLinkedIons(lp, SPLITTYPE_B);

    

    double scoreB = xhhc_scorer.scoreSpectrumVsSeries(spectrum, ion_seriesB);

    if (scoreA > scoreB)
      cout << scoreA << "\t" << scoreB << endl;
    else
      cout << scoreB << "\t" << scoreA << endl;

  } else if (scoremethod=="concatenation") {


    vector<double> scores;
    double score1 = get_concat_score(peptideA, peptideB, posA, charge, spectrum);
    scores.push_back(score1);

    double score2 = get_concat_score(peptideB, peptideA, posB, charge, spectrum);
    scores.push_back(score2);


    int lengthA = string(peptideA).length();
    int lengthB = string(peptideB).length();

    double score3 = get_concat_score(peptideA, peptideB, lengthA + posB, charge, spectrum);
    scores.push_back(score3);

    double score4 = get_concat_score(peptideB, peptideA, lengthB + posA, charge, spectrum);
    scores.push_back(score4);

    sort(scores.begin(), scores.end(), less<double>());
    cout <<scores[0];
    for (int i=1;i<4;i++)
      {
	cout <<"\t"<<scores[i];
      }

    cout << endl;
  }
  else {
    carp(CARP_ERROR,"Unknown method");
  }
  // free heap
  delete collection;
  delete spectrum;
}


double get_concat_score(char* peptideA, char* peptideB, int link_site, int charge, Spectrum* spectrum) {
  string lpeptide = string(peptideA) + string(peptideB); 
  
  IonConstraint* ion_constraint = IonConstraint::newIonConstraintSmart(XCORR, charge);
  
  IonSeries* ion_series = new IonSeries(lpeptide.c_str(), charge, ion_constraint);
  
  ion_series->predictIons();
  
  //modify ions.
  
  //int pepA_begin = 0;
  int pepB_begin = string(peptideA).length();
  int llength = lpeptide.length();
  
  for (IonIterator ion_iterator = ion_series->begin();
    ion_iterator != ion_series->end();
    ++ion_iterator) {

    Ion* ion = *ion_iterator;
    //check to see if if is the cterm of 1st peptide.
    int ion_charge = ion->getCharge();
    int cleavage_idx = ion->getCleavageIdx();
    ION_TYPE_T ion_type = ion->getType();

    //if contains cterm of 1st peptide, modify by -OH 
      
    carp(CARP_DEBUG,"====================");
    if (ion_type == B_ION) {
	carp(CARP_DEBUG,"B-ion");
	carp(CARP_DEBUG,"%s",lpeptide.substr(0,cleavage_idx).c_str());
      } else if (ion_type == Y_ION) {
	carp(CARP_DEBUG,"Y-ion");
	carp(CARP_DEBUG,"%s",lpeptide.substr(llength-cleavage_idx,llength).c_str());
      }
      else continue;

      carp(CARP_DEBUG,"cleavage idx:%d",cleavage_idx);
      //print_ion(ion, stdout);

      bool cterm_1st = false;
      if (ion_type == B_ION) {
	carp(CARP_DEBUG,"B-Ion");
	if (cleavage_idx >= pepB_begin) {
	  cterm_1st = true;
	}
      } else if (ion_type == Y_ION) {
	carp(CARP_DEBUG,"Y-ion");
	if (cleavage_idx > (llength - pepB_begin)) {
	  cterm_1st = true;
	}
      }

      bool nterm_2nd = false;
      if (ion_type == B_ION) {
	if (cleavage_idx > pepB_begin) {
	  nterm_2nd = true;
	}
      } else if (ion_type == Y_ION) {
	if (cleavage_idx >= (llength-pepB_begin)) {
	  nterm_2nd = true;
	}
      }

      bool has_link_site = false;
     if (ion_type == B_ION) {
	if (cleavage_idx > link_site) {
	  has_link_site = true;
	}
      } else if (ion_type == Y_ION) {
       if (cleavage_idx >= (llength- link_site)) {
	  has_link_site = true;
	}
      }
      
     carp(CARP_DEBUG,"cterm:%d",cterm_1st);
     carp(CARP_DEBUG,"nterm:%d",nterm_2nd);
     carp(CARP_DEBUG,"has site:%d",has_link_site);
      

      //if it contains the cterm of the 1st peptide, modify by -OH
      if (cterm_1st) {
	FLOAT_T old_mass = (ion->getMassZ() - MASS_H_MONO) * (FLOAT_T)ion_charge;
	FLOAT_T new_mass = old_mass + MASS_H2O_MONO - MASS_H_MONO;
	FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
	ion->setMassZ(new_mz);
      }
      //if contains the nterm of 2nd peptide, modify by -H
      if (nterm_2nd) {
	FLOAT_T old_mass = (ion->getMassZ() - MASS_H_MONO) * (FLOAT_T)ion_charge;
	FLOAT_T new_mass = old_mass + MASS_H_MONO;
	FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
	ion->setMassZ(new_mz);
      }
      //if contains the link site, modify by link mass.
      if (has_link_site) {
	FLOAT_T old_mass = (ion->getMassZ() - MASS_H_MONO) * (FLOAT_T)ion_charge;
	FLOAT_T new_mass = old_mass + LinkedPeptide::getLinkerMass();
	FLOAT_T new_mz = (new_mass + (FLOAT_T)ion_charge) / (FLOAT_T)ion_charge;
	ion->setMassZ(new_mz);
      }
    

    

      //print_ion(ion, stdout);
    
    }

    Scorer* scorer = new Scorer(XCORR); 

    // calculate the score
    FLOAT_T score = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);
    return score;



}

FLOAT_T* get_observed_raw(Spectrum* spectrum, int charge) {
  FLOAT_T peak_location = 0;
  int mz = 0;
  FLOAT_T intensity = 0;
  FLOAT_T bin_width = bin_width_mono;
  FLOAT_T precursor_mz = spectrum->getPrecursorMz();
  FLOAT_T experimental_mass_cut_off = precursor_mz*charge + 50;

  // set max_mz and malloc space for the observed intensity array
  FLOAT_T sp_max_mz = 512;

  if(experimental_mass_cut_off > 512){
    int x = (int)experimental_mass_cut_off / 1024;
    FLOAT_T y = experimental_mass_cut_off - (1024 * x);
    sp_max_mz = x * 1024;

    if(y > 0){
      sp_max_mz += 1024;
    }
  }

  // DEBUG
  // carp(CARP_INFO, "experimental_mass_cut_off: %.2f sp_max_mz: %.3f", experimental_mass_cut_off, scorer->sp_max_mz);
  FLOAT_T* observed = (FLOAT_T*)mycalloc((int)sp_max_mz, sizeof(FLOAT_T));
  
  // DEBUG
  // carp(CARP_INFO, "max_peak_mz: %.2f, region size: %d",get_spectrum_max_peak_mz(spectrum), region_selector);
  
  for (PeakIterator peak_iterator = spectrum->begin();
    peak_iterator != spectrum->end();
    ++peak_iterator) {

    Peak *peak = *peak_iterator;
    peak_location = peak->getLocation();
    
    // skip all peaks larger than experimental mass
    if(peak_location > experimental_mass_cut_off){
      continue;
    }
    
    // skip all peaks within precursor ion mz +/- 15
    if(peak_location < precursor_mz + 15 &&  peak_location > precursor_mz - 15) {
      continue;
    }
    
    // map peak location to bin
    mz = (int)(peak_location / bin_width + 0.5);

    // get intensity
    // sqrt the original intensity
    intensity = peak->getIntensity();

    // set intensity in array with correct mz, only if max peak in the bin
    if(observed[mz] < intensity){
      observed[mz] = intensity;
      }
    }    
  

  return observed;

}



void print_spectrum(Spectrum* spectrum, LinkedIonSeries& ion_series) {


      Scorer* scorer = new Scorer(XCORR);
      scorer->createIntensityArrayObserved(spectrum, ion_series.getCharge());

      FLOAT_T* observed_raw = get_observed_raw(spectrum, ion_series.getCharge());
      FLOAT_T* observed_processed = scorer->getIntensityArrayObserved();


      FLOAT_T max_mz = scorer->getSpMaxMz();

        


      XHHC_Scorer xhhc_scorer(max_mz);

      FLOAT_T* theoretical = (FLOAT_T*)mycalloc((size_t)max_mz, sizeof(FLOAT_T));
      xhhc_scorer.hhcCreateIntensityArrayTheoretical(ion_series, theoretical);


      
      ofstream fout("ion_match.out");
      for (int i=0;i<max_mz;i++) {
        fout << i << "\t" 
             << observed_raw[i] << "\t"
             << observed_processed[i] << "\t" 
             << theoretical[i] << "\t" 
             << (theoretical[i] != 0) << "\t" 
             << (theoretical[i] > 25) << endl;
      }

      fout.close();
}
