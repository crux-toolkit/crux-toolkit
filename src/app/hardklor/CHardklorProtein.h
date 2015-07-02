#ifndef _CHARDKLORPROTEIN_H
#define _CHARDKLORPROTEIN_H

#include <vector>

#include "HardklorTypes.h"

using namespace std;

class CHardklorProtein {
 public:
  //Constructors & Destructors:
  CHardklorProtein();
  CHardklorProtein(const CHardklorProtein&);
  ~CHardklorProtein();

  //Functions:
  CHardklorProtein& operator=(const CHardklorProtein&);
  void add(int, int, double);
  sEnrichMercury& at(int);
  void clear();
  int size();

  //Data Members:
  sInt ID;
  double mz;
  double monoMass;
  double shft;
  double abun;
  double rTime;
  int charge;
  int C;
  char seq[64];
  
 protected:
 private:
  //Data Members:
  vector<sEnrichMercury> *enrich;

};

#endif
