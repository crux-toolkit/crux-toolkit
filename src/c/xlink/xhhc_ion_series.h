#ifndef XHHC_ION_SERIES_H
#define XHHC_ION_SERIES_H

#include "xhhc.h"

//CRUX includes
#include "ion_series.h"
#include "scorer.h"
#include "spectrum.h"
#include "spectrum_collection.h"

// should this be somewhere else?
#define bin_width_mono 1.0005079

using namespace std;

class LinkedIonSeries {
  public:
    // constructors
    LinkedIonSeries();

    LinkedIonSeries(int charge);
    //LinkedIonSeries(char* sequenceA, char* sequenceB, int posA, int posB, int charge);
    //LinkedIonSeries(char* sequenceA, char* sequenceB, char* links, int charge);

    // getters
    int charge()                  { return charge_; }
    vector<LinkedPeptide>& ions() { return all_ions; }
    int size()                    { return all_ions.size(); }

    int get_total_by_ions();

    int get_observable_ions(
      FLOAT_T min_mz,
      FLOAT_T max_mz,
      FLOAT_T bin_width,
      int& ions_observable,
      int& ions_observable_bin);

    
    int get_observable_by_ions(
      FLOAT_T min_mz, 
      FLOAT_T max_mz, 
      FLOAT_T bin_width,
      int &by_observable,
      int &by_observable_bin);


    // other
    void set_charge(int charge) { charge_ = charge; }
    void clear()                { all_ions.clear(); }
    // splits
    void add_linked_ions(LinkedPeptide& linked_peptide, int split_type=0);
    // print tab-delimited list of all ions
    void print();



  private:
    
    int charge_; 
    // a list of all the ions 
    std::vector<LinkedPeptide> all_ions;

    //for add linked ions.
    std::vector<pair<LinkedPeptide, LinkedPeptide> > fragments;

    MASS_TYPE_T fragment_mass_type;

    //for ion mass matrix.
    //std::vector<FLOAT_T> mass_matrix;
};


void hhc_predict_ions(
  ION_SERIES_T* ion_series, ///< the ion series to predict ions for -in
  FLOAT_T linker_mass,
  int linker_site);

#endif
