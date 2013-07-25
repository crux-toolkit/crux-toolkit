/**
 * \file LinkedIonSeries.h 
 * AUTHOR: Sean McIlwain and Paul Draghicescu
 * CREATE DATE: 18 January 2012
 * \brief Object for keeping track of a (non-)crosslinked peptide ions.
 *****************************************************************************/
#ifndef LINKEDIONSERIES_H
#define LINKEDIONSERIES_H

#include "xhhc.h"

//CRUX includes
#include "IonSeries.h"
#include "Scorer.h"
#include "Spectrum.h"
#include "SpectrumCollection.h"

// TODO - Get rid of this dependency. (xhhc_scorer and xhhc_search).
#define bin_width_mono 1.0005079

class LinkedIonSeries {
 protected:

  int charge_; ///< maximum charge of the ion series
  std::vector<LinkedPeptide> all_ions_; ///< list of all ions
  std::vector<pair<LinkedPeptide, LinkedPeptide> > fragments_; ///< for add linked ions
  MASS_TYPE_T fragment_mass_type_; ///< mass type of the fragments.

  /**
   * adds an ion to the observed list, used in conjunction with 
   * getObservableIons
   */
  void addIonBin(
    map<int, bool>& observed, ///< the observed bin vector -in/out
    int& ions, ///< the number of ions -in/out
    int& ions_bin, ///< the number of binned ions -in/out
    FLOAT_T mz, ///< m/z of the ion
    FLOAT_T bin_width, ///< bin-width
    FLOAT_T min_mz, ///< max m/z for the range
    FLOAT_T max_mz, ///< min m/z for the range
    bool add_flanks ///< add flanks?
    );

 public:

  /**
   * \returns a blank LinkedIonSeries object
   */
  LinkedIonSeries();

  /**
   * \returns an ion series object assigning the charge.
   */
  LinkedIonSeries(
    int charge ///< The maximum charge of the ion series.
  );

  /**
   * Default destructor
   */
  virtual ~LinkedIonSeries();

  /**
   * \returns charge for this linked ion series
   */
  int getCharge();

  /**
   * \returns all the ions in this series
   */
  std::vector<LinkedPeptide>& getIons();

  /**
   * \returns the number of ions within this series
   */
  int getSize();

  /**
   * \returns the total number of by ions for this series
   */
  int getTotalBYIons();

  /**
   * \returns the total number of ions observable within a binned
   * range
   */
  int getObservableIons(
    FLOAT_T min_mz, ///< the minimum m/z -in
    FLOAT_T max_mz, ///< the maximum m/z -in
    FLOAT_T bin_width, ///< the width of the bins -in
    int& ions_observable, ///<number of ions observable -out
    int& ions_observable_bin ///<number of bins observable -out
  );

  /**
   * \returns the total number of b-y ions observable within a binned
   * range
   */  
  int getObservableBYIons(
    FLOAT_T min_mz, ///< the minimum m/z -in
    FLOAT_T max_mz, ///< the maximum m/z -in
    FLOAT_T bin_width, ///< the width of the bins -in
    int &by_observable, ///<number of ions observable -out
    int &by_observable_bin ///<number of bins observable -out
  );

  /**
   * sets the charge for this ion series
   */
  void setCharge(
    int charge ///< the charge
    );
  
  /**
   * removes all ions from this series
   */
  void clear();

  /**
   * Adds the linked ions for a linked peptide
   */
  void addLinkedIons(
    LinkedPeptide& linked_peptide, ///< The linked peptide 
    SPLITTYPE_T split_type=SPLITTYPE_BOTH ///< Which peptide to split
  );

  /**
   * prints out tab delimited information about the ion series
   */
  void print();


};
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

