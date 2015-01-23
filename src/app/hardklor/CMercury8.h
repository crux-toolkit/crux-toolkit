#ifndef _CMercury8_H
#define _CMercury8_H


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "ctype.h"
#include <ctime>
#include "mercury.h"
#include "FFT.h"
using namespace std;

typedef struct
{
   char  Symbol[3];	/* Elemental symbol */
   int	 NumIsotopes;	/* Number of stable isotopes */
   float *IsoMass;	/* Array of isotopic masses */
   int   *IntMass;	/* Array of integer isotopic masses */
   float *IsoProb;	/* Array of isotopic probabilities */
   int   NumAtoms;	/* Number of occurances of element in molecular formula */
 
} Atomic5;

class CMercury8 {
 private:
  //Data Members:
  Atomic5 Element[MAXAtomNo+1];	/* 104 elements allows for Z=103 or Lr */
  Atomic5 Orig[MAXAtomNo+1];
  int AtomicNum[MAXIsotopes];	/* Atomic numbers of elements parsed from molecular formula */
  bool showOutput;
  bool bAccMass;
  bool bRelAbun;
  vector<int> EnrichAtoms;
  double monoMass;
  double zeroMass;

  //Functions:
  void AccurateMass(int,int);
  void AddElement(char[],int,int);
    void CalcFreq(Hardklor::complex*, int, int, int, int);
  void CalcMassRange(int*, double, int, int);
  void CalcVariances(double*, double*, int);
  void CalcWeights(double&,double&,double&,int&,int&,int&,int&,int);
    void ConvertMass(Hardklor::complex*, int, int, double, double, int, int, int, double, double);
  void GetPeaks(Hardklor::complex*, int, vector<Result>&, int, int);
  void InitializeData(const char* fn=NULL);
  void InitializeDataHardcoded();
  void MassToInt(Hardklor::complex*, int);
  void Mercury(int,int);
  int ParseMF(char[], int*);
  void RelativeAbundance(vector<Result>&);
 
 public:
  //Data Members:
  vector<Result> FixedData;
  vector<Result> FracAbunData;

  //Constructors & Destructors:
  CMercury8();
  CMercury8(char *fn);
  ~CMercury8();

  //Functions:
  void AccMass(bool);
  void Echo(bool);
  void Enrich(int,int,double d=0.99);
  double getMonoMass();
  double getZeroMass();
  int GoMercury(char*, int=1, const char* filename="\0");
  void Intro();
  void RelAbun(bool);
  void Reset();

};

#endif
