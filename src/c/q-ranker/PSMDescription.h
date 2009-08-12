#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <set>
#include <string>
using namespace std;

namespace qranker {

class PSMDescription
{
public:
  PSMDescription();
  virtual ~PSMDescription();
  void clear() {proteinIds.clear();}
  void calcRegressionFeature();        

  double q,pep,sc;
  double * features;
  double * retentionFeatures;
  double retentionTime,massDiff,pI;
  unsigned int scan; 
  set<string> proteinIds;
  string peptide;  
};

} // qranker namspace

#endif /*PSMDESCRIPTION_H_*/
