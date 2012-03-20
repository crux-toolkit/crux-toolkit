#ifndef _CHARDKLORVARIANT_H
#define _CHARDKLORVARIANT_H

#include <vector>

#include "HardklorTypes.h"

using namespace std;

class CHardklorVariant {
 public:
  //Constructors & Destructors:
  CHardklorVariant();
  CHardklorVariant(const CHardklorVariant&);
  ~CHardklorVariant();

  //Operators:
  CHardklorVariant& operator=(const CHardklorVariant&);

  //Methods:
  void addAtom(const sInt&);
  void addAtom(const int&, const int&);
  void addEnrich(const sEnrichMercury&);
  void addEnrich(const int&, const int&, const double&);
  sInt& atAtom(int);
  sEnrichMercury& atEnrich(int);
  void clear();
  int sizeAtom();
  int sizeEnrich();

  //Data Members
  int ID;

 protected:

 private:
  //Data Members:
  vector<sInt> *atoms;
  vector<sEnrichMercury> *enrich;

};

#endif
