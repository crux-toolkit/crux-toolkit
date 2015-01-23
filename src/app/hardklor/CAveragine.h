#ifndef _CAVERAGINE_H
#define _CAVERAGINE_H

#include "CHardklorVariant.h"
#include "CPeriodicTable.h"

#include <fstream>
#include <vector>

using namespace std;

/*
const double AVE_MASS = 111.1254;
const double AVE_C = 4.9384;
const double AVE_H = 7.7583;
const double AVE_N = 1.3577;
const double AVE_O = 1.4773;
const double AVE_S = 0.0417;
*/

const double AVE_MASS = 111.2137;
const double AVE_C = 4.9558;
const double AVE_H = 7.8241;
const double AVE_N = 1.3571;
const double AVE_O = 1.4716;
const double AVE_S = 0.0390;

typedef struct {
  char symbol[3];
  int numIsotopes;
  vector<double> *mass;
  vector<double> *abundance;
} atomInfo;

class CAveragine {
  
 public:
  //Constructors & Destructors
  //CAveragine();
  CAveragine(char* fn=NULL, char* fn2=NULL);
  ~CAveragine();

  //Methods:
  void calcAveragine(double,CHardklorVariant);
  void clear();
  void getAveragine(char*);
  void setAveragine(int,int,int,int,int);
  int getElement(int);
  double getMonoMass();
  void loadTable(char*);
  void loadTableHardcoded();

 protected:

 private:
  //Data Members:
  //double monoMass;
  int *atoms;
  CPeriodicTable *PT;
  vector<atomInfo> *enrich;

};

#endif
