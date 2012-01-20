#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"

#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

#include <iostream>

FLOAT_T* mass_matrix = NULL;
int mass_matrix_size = -1;
// main constructor


LinkedIonSeries::LinkedIonSeries() {
  charge_ = 0;
  fragment_mass_type = get_mass_type_parameter("fragment-mass");
}

LinkedIonSeries::LinkedIonSeries(int charge) {
  charge_ = charge;
  fragment_mass_type = get_mass_type_parameter("fragment-mass");
}

// prints out tab delimited information about the ion series
void LinkedIonSeries::print() {

  cout <<"Sorting ions"<<endl;

  LinkedPeptide::sortByMass(all_ions, fragment_mass_type);

  cout << "m/z\ttype\tion" << endl;
  string ion_type;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    if (ion->getIonType() == B_ION) 
      ion_type = "B_ION";
    else 
      ion_type = "Y_ION";
    cout << ion->getMZ(fragment_mass_type) << "\t" << ion_type << "\t" << *ion << endl;
  }
}

#define SPLIT_BOTH 0
#define SPLIT_A 1
#define SPLIT_B 2

// cleaves linked_peptide at all positions, adding b and y ions
void LinkedIonSeries::add_linked_ions(LinkedPeptide& linked_peptide,int split_type) {

  //if (charge_ == 0) charge_ = linked_peptide.charge();
  linked_peptide.setCharge(charge_);

  fragments.clear();
  // split the precursor at every cleavage site

  switch(split_type) {
  case SPLIT_A:
    linked_peptide.splitA(fragments);
    break;
  case SPLIT_B:
    linked_peptide.splitB(fragments);
    break;
  case SPLIT_BOTH:
  default:
    linked_peptide.split(fragments);
  }
  
  for (vector<pair<LinkedPeptide, LinkedPeptide> >::iterator ion_pair = 
    fragments.begin(); 
    ion_pair != fragments.end(); 
    ++ion_pair) {

    // if b-ion and not a neutral loss
    if (ion_pair->first.getCharge() != 0) {
      ion_pair->first.setIonType(B_ION); 
      ion_pair->first.calculateMass(fragment_mass_type);
      //cout << ion_pair->first.get_mz() << " B " << ion_pair->first << endl;
      all_ions.push_back(ion_pair->first);
    }
    // if y-ion and not a neutral loss
    if (ion_pair->second.getCharge() != 0) {
      ion_pair->second.setIonType(Y_ION); 
      ion_pair->second.calculateMass(fragment_mass_type);
      //cout << ion_pair->second.get_mz() << " Y " << ion_pair->second << endl;
      all_ions.push_back(ion_pair->second);
    }
  }
  
}

int LinkedIonSeries::get_total_by_ions() {
 
  
  int ans = 0;
  vector<LinkedPeptide>::iterator ion_iter;

  for (ion_iter = all_ions.begin(); 
       ion_iter != all_ions.end(); 
       ++ion_iter) {
        //if (ion_iter -> get_mz(MONO) >= 400 && ion_iter -> get_mz(MONO) <= 1200) {
	  if (ion_iter -> getIonType() == B_ION || ion_iter -> getIonType() == Y_ION) {
	    ans++;
	  }
	//}
  }
  return ans;
}

void add_ion_bin(map<int, bool>& observed, 
                int& ions, int& ions_bin, 
                FLOAT_T mz, FLOAT_T bin_width,
                FLOAT_T min_mz, FLOAT_T max_mz, bool add_flanks) {

  if ((mz >= min_mz) && (mz <= max_mz)) {
    ions++;
    int bin_idx = (int)(mz / bin_width + 0.5);
    if (observed.find(bin_idx) == observed.end()) {
      observed[bin_idx] = true;
      ions_bin++;
    }
    if (add_flanks) {
      
      bin_idx = bin_idx - 1;
      FLOAT_T flank_mz = bin_idx * bin_width;
      if (flank_mz >= min_mz) {
        ions++;
        if (observed.find(bin_idx) == observed.end()) {
          observed[bin_idx] = true;
          ions_bin++;
        }
      }
      bin_idx = bin_idx + 2;
      flank_mz = bin_idx * bin_width;
      if (flank_mz <= max_mz) {
        ions++;
        if (observed.find(bin_idx) == observed.end()) {
          observed[bin_idx] = true;
          ions_bin++;
        }
      }
    }
  }
}

int LinkedIonSeries::get_observable_ions(
  FLOAT_T min_mz,
  FLOAT_T max_mz,
  FLOAT_T bin_width,
  int& ions_observable,
  int& ions_observable_bin) {
  
  ions_observable = 0;
  ions_observable_bin = 0;

  map<int, bool> observed;

  for (vector<LinkedPeptide>::iterator ion_iter = all_ions.begin();
    ion_iter != all_ions.end();
    ++ion_iter) {

    double mz = ion_iter -> getMZ(MONO);

    if (ion_iter -> getIonType() == B_ION || 
      ion_iter -> getIonType() == Y_ION) {
      add_ion_bin(observed, ions_observable, ions_observable_bin, mz, bin_width, min_mz, max_mz, true);
      double h2O_mz = ion_iter -> getMZ(MONO) - (MASS_H2O_MONO/ ion_iter -> getCharge());
      add_ion_bin(observed, ions_observable, ions_observable_bin, h2O_mz, bin_width, min_mz, max_mz, false);
      double nh3_mz = ion_iter -> getMZ(MONO) - (MASS_NH3_MONO/ ion_iter -> getCharge());
      add_ion_bin(observed, ions_observable, ions_observable_bin, nh3_mz, bin_width, min_mz, max_mz, false);
      double co_mz = ion_iter -> getMZ(MONO) - (MASS_CO_MONO / ion_iter -> getCharge());
      add_ion_bin(observed, ions_observable, ions_observable_bin, co_mz, bin_width, min_mz, max_mz, false);
    }
  }
  return 0;
}


int LinkedIonSeries::get_observable_by_ions(
  FLOAT_T min_mz, 
  FLOAT_T max_mz, 
  FLOAT_T bin_width,
  int &by_observable,
  int &by_observable_bin) {

  by_observable = 0;
  by_observable_bin = 0;

  map<int, bool> observed;

  vector<LinkedPeptide>::iterator ion_iter;

  for (ion_iter = all_ions.begin(); 
       ion_iter != all_ions.end(); 
       ++ion_iter) {
    if (ion_iter -> getMZ(MONO) >= min_mz && ion_iter -> getMZ(MONO) <= max_mz) {
      if (ion_iter ->getIonType() == B_ION || ion_iter -> getIonType() == Y_ION) {
        by_observable++;
        int bin_idx = 
          (int)(ion_iter->getMZ(MONO) / bin_width + 0.5);
        if (observed.find(bin_idx) == observed.end()) {
          observed[bin_idx] = true;
          by_observable_bin++;
        }
      }
    }
  }
  return by_observable_bin;
}



/***************************************
 *CRUX OVERRIDES
 ***************************************/


/**
 * \brief Creates an array in which element i is the sum of the masses
 * of amino acids 0 to (i-1).  At i=0 is stored the length of the
 * peptide.  
 * \returns an array of ion masses for all sub sequences
 */
void hhc_create_ion_mass_matrix(
  //char* peptide, ///< The peptide for this ion series. -in
  MODIFIED_AA_T* modified_seq, ///< the sequence
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int peptide_length, ///< the length of the peptide
  FLOAT_T linker_mass,
  int linker_site
  )
{
  if( modified_seq == NULL ){
  //if( peptide == NULL ){
    carp(CARP_ERROR, "Cannot create mass matrix from NULL seqence");
  }

  if (mass_matrix_size < peptide_length+1) {
    if (mass_matrix != NULL) {
      free(mass_matrix);
    }
    mass_matrix = (FLOAT_T*)mymalloc(sizeof(FLOAT_T)*(peptide_length+1));
    mass_matrix_size = peptide_length+1;
  }

  // at index 0, the length of the peptide is stored
  mass_matrix[0] = peptide_length;

  // add up AA masses
  int ion_idx = 1;
  // initialize first to be mass of c-term amino acid
  // mass_matrix[ion_idx] = get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  mass_matrix[ion_idx] = get_mass_mod_amino_acid(modified_seq[ion_idx-1], mass_type);
  // for open modification cross linking
  if (linker_site == 0) mass_matrix[ion_idx] += linker_mass; 
  //++ion_idx;
  //for(; ion_idx <= peptide_length; ++ion_idx){
  for(ion_idx = 2; ion_idx <= peptide_length; ++ion_idx){
    mass_matrix[ion_idx] = mass_matrix[ion_idx-1] + 
      get_mass_mod_amino_acid(modified_seq[ion_idx-1], mass_type);
    if (linker_site == ion_idx-1)
      mass_matrix[ion_idx] += linker_mass; 
      //get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  }
  // DEBUGGING
  /*
  fprintf(stderr, "seq:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%s", modified_aa_to_string(modified_seq[ion_idx]));
  }
  fprintf(stderr, "\nmas:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%.2f", get_mass_mod_amino_acid(modified_seq[ion_idx], MONO));
  }
  fprintf(stderr, "\nsum:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%.2f", mass_matrix[ion_idx+1]);
  }
  fprintf(stderr, "\n");
  */
}

/**
 * \brief The engine of ion series. Predicts all the ions from the
 * peptide that meet the ion constraint. All predicted ions are stored
 * in the ion_series as ion objects. 
 */
void hhc_predict_ions(
  IonSeries* ion_series, ///< the ion series to predict ions for -in
  FLOAT_T linker_mass,
  int linker_site
  )
{


  if(ion_series->getIsPredicted()){
    carp(CARP_WARNING, "The ion series has already been predicted and added");
    return;
  }


  IonConstraint* constraint = ion_series->getIonConstraint();
  
  // create a mass matrix
 
    hhc_create_ion_mass_matrix(ion_series->getModifiedAASeq(), 
			       constraint->getMassType(), 
			       ion_series->getPeptideLength(), 
			       linker_mass, linker_site);  
  /*
  printf("cumulative mass sum is:\n");
  int idx = 0;
  for(idx = 0; idx < mass_matrix[0]; idx++){
    printf("%i\t%f\n", idx, mass_matrix[idx]);
  }
  */
  // scan for the first and last  (S, T, E, D) and (R, K, Q, N), 
  // initialize to determine modification is ok.
  // the first, last of STED, RKQN are stored in ion_series.
  ion_series->scanForAAForNeutralLoss();
  
  // generate ions without any modifications
  if(!ion_series->generateIonsNoModification(mass_matrix)){
    carp(CARP_FATAL, "failed to generate ions, no modifications");
  }

  // create modification ions?
  if(constraint->getUseNeutralLosses()){
    
    // generate ions with nh3 modification
    if(abs(constraint->getModification(NH3)) > 0){
      if(!ion_series->generateIons(NH3)){
        carp(CARP_FATAL, "failed to generate ions, NH3 modifications");
      }
    }
    
    // generate ions with h2o modification
    if(abs(constraint->getModification(H2O)) > 0){
      if(!ion_series->generateIons(H2O)){
        carp(CARP_FATAL, "failed to generate ions, H2O modifications");
      }
    }

    // generate ions with isotope modification
    if(constraint->getModification(ISOTOPE) > 0){
      if(!ion_series->generateIons(ISOTOPE)){
        carp(CARP_FATAL, "failed to generate ions, ISOTOPE modifications");
      }
    }

    // generate ions with flank modification
    if(constraint->getModification(FLANK) > 0){
      if(!ion_series->generateIonsFlank()){
        carp(CARP_FATAL, "failed to generate ions, FLANK modifications");
      }
    }
    
    // add more modifications here

  }
  
  // ion series now been predicted
  ion_series->setIsPredicted(true);

}
