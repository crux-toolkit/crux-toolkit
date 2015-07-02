#ifndef _CAVERAGINE_H
#define _CAVERAGINE_H

#include "CHardklorVariant.h"
#include "CPeriodicTable.h"

#include <cstdlib>
#include <fstream>
#include <vector>
#include <cstring>

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

typedef struct atomInfo {
  char symbol[3];
  int numIsotopes;
  vector<double> *mass;
  vector<double> *abundance;
  atomInfo(){
    strcpy(symbol,"X");
    numIsotopes=0;
    mass = new vector<double>;
    abundance = new vector<double>;
  }
  atomInfo(const atomInfo& a){
    strcpy(symbol,a.symbol);
    numIsotopes=a.numIsotopes;
    mass = new vector<double>;
    abundance = new vector<double>;
    unsigned int i;
    for(i=0;i<a.mass->size();i++) mass->push_back(a.mass->at(i));
    for(i=0;i<a.abundance->size();i++) abundance->push_back(a.abundance->at(i));
  }
  ~atomInfo(){
    delete mass;
    delete abundance;
  }
  atomInfo& operator=(const atomInfo& a){
    if(&a!=this){
      strcpy(symbol,a.symbol);
      numIsotopes=a.numIsotopes;
      delete mass;
      delete abundance;
      mass = new vector<double>;
      abundance = new vector<double>;
      unsigned int i;
      for(i=0;i<a.mass->size();i++) mass->push_back(a.mass->at(i));
      for(i=0;i<a.abundance->size();i++) abundance->push_back(a.abundance->at(i));
    }
    return *this;
  }

} atomInfo;

class CAveragine {
  
 public:
  //Constructors & Destructors
  //CAveragine();
  CAveragine(char* fn="ISOTOPE.DAT", char* fn2="Hardklor.dat");
  ~CAveragine();

  //Methods:
  void calcAveragine(double,CHardklorVariant);
  void clear();
  void defaultValues();
  void getAveragine(char*);
  void setAveragine(int,int,int,int,int);
  int getElement(int);
  double getMonoMass();
  CPeriodicTable* getPT();
  void loadTable(char*);

 protected:

 private:
  //Data Members:
  //double monoMass;
  int *atoms;
  CPeriodicTable *PT;
  vector<atomInfo> *enrich;

};

#endif
