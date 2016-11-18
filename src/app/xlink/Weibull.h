
#include <set>
#include <vector>
#include <string>
#include "model/objects.h"

class Weibull {
    
 protected:
  FLOAT_T correlation_;
  FLOAT_T eta_;
  FLOAT_T beta_;
  FLOAT_T shift_;
  
  std::set<std::string> sequences_;
  std::vector<FLOAT_T> scores_;
  
  bool fit_called_;
  bool fit_success_;

  int duplicates_;
  
 public:
  Weibull();
  ~Weibull();
  
  FLOAT_T getEta();
  FLOAT_T getBeta();
  FLOAT_T getShift();
  FLOAT_T getCorrelation();
  
  void reset();
  void addPoint(const std::string& sequence, FLOAT_T score);
  bool fit();

  FLOAT_T getWeibullPValue(FLOAT_T score, bool logp=false);
  FLOAT_T getECDFPValue(FLOAT_T score, bool logp=false);
  FLOAT_T getPValue(FLOAT_T score, bool logp=false);


    
};
