#include "XLinkIonSeriesCache.h"

using namespace std;

vector<vector<IonSeries*> > XLinkIonSeriesCache::target_xlinkable_ion_series_;

vector<vector<IonSeries*> > XLinkIonSeriesCache::decoy_xlinkable_ion_series_;
vector<IonConstraint*> XLinkIonSeriesCache::xcorr_ion_constraint_;


IonSeries* XLinkIonSeriesCache::getXLinkablePeptideIonSeries(
  XLinkablePeptide& xpep,
  int charge
  ) {

  IonSeries* ans = NULL;
  int xpep_idx = xpep.getIndex();

  if (xpep_idx == -1) {
    //carp(CARP_DEBUG, "Unindexed xlinkable peptide. Returning NULL");
    return NULL;
  } else {

    bool decoy = xpep.isDecoy();
    //carp(CARP_INFO, "decoy %i pep_idx %i charge %i", decoy, xpep_idx, charge);
    vector<vector<IonSeries*> >* ion_cache = &target_xlinkable_ion_series_;
    if (decoy) {
      ion_cache = &decoy_xlinkable_ion_series_;
      //carp(CARP_INFO, "Getting decoy cache");
    }

    int charge_idx = charge-1;
    //carp(CARP_INFO, "getting pepidx:%i",xpep_idx);
    while(ion_cache->size() <= xpep_idx) {
      //carp(CARP_INFO, "Adding vector<IonSeries*>():%d", ion_cache->size());
      ion_cache->push_back(vector<IonSeries*>());
    }
  
    vector<IonSeries*>& level1 = (*ion_cache)[xpep_idx];

    if (level1.size() == 0) {
      IonSeries* ion_series1 = new IonSeries(getXCorrIonConstraint(1), 1);
      ion_series1->update(xpep.getSequence(), xpep.getModifiedSequencePtr());
      ion_series1->predictIons();
      //carp(CARP_INFO, "charge1:%d",ion_series1->getNumIons());
      level1.push_back(ion_series1);
      IonSeries* ion_series2 = new IonSeries();
      IonSeries::copy(ion_series1, ion_series2, false);
      ion_series2->setCharge(2);
      ion_series2->setIonConstraint(getXCorrIonConstraint(2));
      level1.push_back(ion_series2);

    }
    while(level1.size() <= charge_idx) {
      IonSeries* next_series = new IonSeries();
      IonSeries* current_series = level1.back();
      int next_charge = current_series->getCharge()+1;
      
      IonSeries::copy(current_series, next_series, false);
      //carp(CARP_INFO, "copy src:%d dest:%d", current_series->getNumIons(), next_series->getNumIons());
      next_series->setCharge(next_charge);
      next_series->setIonConstraint(getXCorrIonConstraint(next_charge));

      for (IonIterator ion_iter = level1[0]->begin();
	   ion_iter != level1[0]->end();
	   ++ion_iter) {
	Ion* newIon = Ion::newIon();
	Ion::copy(*ion_iter, newIon, (*ion_iter)->getPeptideSequence());
	FLOAT_T ion_mass = newIon->getMassFromMassZ();
	newIon->setCharge(next_charge-1);
	newIon->setMassZFromMass(ion_mass);
	next_series->addIon(newIon);
      }
      /*
      IonSeries* test_series = new IonSeries(getXCorrIonConstraint(next_charge), next_charge);
      test_series->update(xpep.getSequence(), xpep.getModifiedSequencePtr());
      test_series->predictIons();
      carp(CARP_INFO, "charge:%d next:%d test:%d", next_charge, next_series->getNumIons(), test_series->getNumIons());
      delete test_series;
      
      for (IonIterator iter = test_series->begin();
	   iter != test_series->end();
	   ++iter) {
	//carp(CARP_INFO, "charge:%d", (*iter)->getCharge());
      }
      */
      level1.push_back(next_series);
    }

    ans = level1[charge_idx];
    if (ans == NULL) {
      carp(CARP_FATAL, "Null ion series?");
    } else {
      //carp(CARP_INFO, "using cached ions");
    }
  }
  return ans;

}

IonConstraint* XLinkIonSeriesCache::getXCorrIonConstraint(
  int charge
  ) {

  int charge_idx = charge - 1;

  while(xcorr_ion_constraint_.size() <= charge_idx) {
    xcorr_ion_constraint_.push_back(IonConstraint::newIonConstraintSmart(XCORR, (xcorr_ion_constraint_.size()+1)));
  }
  //carp(CARP_INFO, "returning ion_constraint");
  return(xcorr_ion_constraint_[charge_idx]);

}


 
void XLinkIonSeriesCache::finalize() {
  for (size_t xpep_idx = 0;xpep_idx <  target_xlinkable_ion_series_.size(); xpep_idx++) {

    vector<IonSeries*> &level1 = target_xlinkable_ion_series_[xpep_idx];
    for (size_t charge_idx=0;charge_idx < level1.size(); charge_idx++) {
      if (level1[charge_idx]) {
	IonSeries::freeIonSeries(level1[charge_idx]);
      }
    }
  }

  for (size_t xpep_idx = 0;xpep_idx < decoy_xlinkable_ion_series_.size(); xpep_idx++) {
    vector<IonSeries*>& level1 = decoy_xlinkable_ion_series_[xpep_idx];
    for (size_t charge_idx=0;charge_idx < level1.size(); charge_idx++) {
      if (level1[charge_idx]) {
	IonSeries::freeIonSeries(level1[charge_idx]);
      }
    }
  }

  for (size_t charge_idx=0;charge_idx < xcorr_ion_constraint_.size();charge_idx++) {
    IonConstraint::free(xcorr_ion_constraint_[charge_idx]);
  }
  
}
