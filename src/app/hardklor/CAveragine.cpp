#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "CAveragine.h"

using namespace std;

CAveragine::CAveragine(char* fn, char* fn2){
	PT = new CPeriodicTable(fn2);
  atoms = new int[PT->size()];
  for(int i=0;i<PT->size();i++) atoms[i]=0;
  enrich = new vector<atomInfo>;
  loadTable(fn);
}

CAveragine::~CAveragine(){
  delete [] atoms;
  for(unsigned int i=0;i<enrich->size();i++){
    delete enrich->at(i).abundance;
    delete enrich->at(i).mass;
  }
  delete enrich;
}

void CAveragine::calcAveragine(double dMass, CHardklorVariant hv) {
  int i,j;
  int iso;
  int atom;
  double fixedD = dMass;
  double units;
  double *enr;
  double *abun;
  double maxAbun;
  double totAbun;
  enr = new double[PT->size()];
  for(i=0;i<PT->size();i++) enr[i]=0;

  for(i=0;i<hv.sizeEnrich();i++){
   
    atom = hv.atEnrich(i).atomNum;
    iso = enrich->at(atom).numIsotopes;
    abun = new double[iso];

    //Set our abundance array and find maximum
    maxAbun=0;
    for(j=0;j<iso;j++){
      abun[j]=enrich->at(atom).abundance->at(j);
      if(abun[j]>maxAbun) maxAbun=abun[j];
    }

    //Normalize to 1
    for(j=0;j<iso;j++) abun[j]/=maxAbun;
    
    //Normalize, then Calculate enrichment
    totAbun=0;
    enr[atom]=0;
    for(j=0;j<iso;j++){
      abun[j]/=maxAbun;
      if(j==hv.atEnrich(i).isotope) abun[j]=(1-hv.atEnrich(i).ape)*abun[j]+hv.atEnrich(i).ape;
      else abun[j]=(1-hv.atEnrich(i).ape)*abun[j];
      totAbun+=abun[j];
      enr[atom]+=abun[j]*enrich->at(atom).mass->at(j);
    }

    //Calculate average mass
    enr[atom] /= totAbun;
    enr[atom] -= PT->at(atom).mass;

    delete [] abun;
  }

  for(i=0;i<hv.sizeAtom();i++) {
    
    atoms[hv.atAtom(i).iLower]+=hv.atAtom(i).iUpper;
    fixedD -= hv.atAtom(i).iUpper*(PT->at(hv.atAtom(i).iLower).mass + enr[hv.atAtom(i).iLower]);
  }

  units = (fixedD)/(AVE_MASS + enr[1]*AVE_H + enr[6]*AVE_C + enr[7]*AVE_N + enr[8]*AVE_O + enr[16]*AVE_S);

  //Quick fix; assumes complete periodic table
  atoms[6] += (int)(AVE_C*units+0.5);
  atoms[8] += (int)(AVE_O*units+0.5);
  atoms[7] += (int)(AVE_N*units+0.5);
  atoms[16] += (int)(AVE_S*units+0.5);

  atoms[1] += (int)(fixedD-( ((int)(AVE_C*units+0.5)) * (PT->at(6).mass + enr[6]) + 
			     ((int)(AVE_N*units+0.5)) * (PT->at(7).mass + enr[7]) + 
			     ((int)(AVE_O*units+0.5)) * (PT->at(8).mass + enr[8]) +  
			     ((int)(AVE_S*units+0.5)) * (PT->at(16).mass + enr[16]) )+0.5);

  delete [] enr;

}

void CAveragine::clear(){
  delete [] atoms;
  atoms = new int[PT->size()];
  for(int i=0;i<PT->size();i++) atoms[i]=0;
}

int CAveragine::getElement(int i){
  return atoms[i];
}

void CAveragine::getAveragine(char *str){
  char tstr[10];
  int i;
  str[0]=0;

  //Create formula string
  for(i=0;i<PT->size();i++){
    if(atoms[i]>0){
      sprintf(tstr,"%s%i",PT->at(i).symbol,atoms[i]);
      strcat(str,tstr);
    }
  }

}

double CAveragine::getMonoMass(){
  double d=0;
  int i;

  for(i=0;i<PT->size();i++) d+=atoms[i]*PT->at(i).mass;
  return d;

}


void CAveragine::loadTableHardcoded() {
  cerr << "CAveragine::loadTableHardcoded"<<endl;

}

void CAveragine::loadTable(char* c){


  if (c == NULL) {
    loadTableHardcoded();
  } else {

    FILE *f;
    int  i,j;
    atomInfo a;
          double d=0;

    f = fopen(c,"rt");
    if(f==NULL) {
      cout << "Cannot read " << c << endl;
      loadTableHardcoded();
      return;
    }
  
    i=0;
    while(!feof(f)){
  
      enrich->push_back(a);
                  fscanf(f,"%2s\t%d\n",enrich->at(i).symbol,&enrich->at(i).numIsotopes);
      enrich->at(i).mass = new vector<double>;
      enrich->at(i).abundance = new vector<double>;
  
      for(j=0;j<enrich->at(i).numIsotopes;j++){
                          enrich->at(i).mass->push_back(d);
                          enrich->at(i).abundance->push_back(d);
                          fscanf(f,"%lf\t%lf\n",&enrich->at(i).mass->at(j),&enrich->at(i).abundance->at(j));
      }
  
                  fscanf(f," \n");
      i++;
  
    }
  
    fclose(f);
  
  }
}

//This is not to be used normally - just here to make my life easier.
void CAveragine::setAveragine(int C, int H, int N, int O, int S){
  atoms[1]=H;
  atoms[6]=C;
  atoms[7]=N;
  atoms[8]=O;
  atoms[16]=S;
}
