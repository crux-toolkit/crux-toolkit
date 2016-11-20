
#include "model/IonSeries.h"
#include "XLinkablePeptide.h"
#include "model/IonConstraint.h"

#include <vector>

class XLinkIonSeriesCache {

 protected:

  //key is: xpep.getIndex(), charge
  static std::vector<std::vector<IonSeries*> > target_xlinkable_ion_series_; 
  static std::vector<std::vector<IonSeries*> > decoy_xlinkable_ion_series_; 

  static std::vector<IonConstraint*> xcorr_ion_constraint_;

 public:

  static IonSeries* getXLinkablePeptideIonSeries(
    XLinkablePeptide& xpep,
    int charge
    );

  static IonConstraint* getXCorrIonConstraint(int charge);

  static void finalize();

};
