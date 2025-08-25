#ifndef _SPECANALYZESUPPORT_H
#define _SPECANALYZESUPPORT_H

#include <vector>
#include "CAveragine.h"
#include "CMercury8.h"
#include "HardklorTypes.h"

/* ******************************************************* */
/*                   CPeakPrediction                       */
/* ******************************************************* */
class CPeakPrediction{
 private:
  //Data Members
  double mz;                //Peak m/z in spectrum
  float intensity;          //Peak intensity in spectrum
  std::vector<int> *charges;     //All charge state predictions for this peak
 protected:
 public:
  //Constructors and Destructors
  CPeakPrediction();
  CPeakPrediction(const CPeakPrediction& p);
  ~CPeakPrediction();

  //Overloaded operators
  CPeakPrediction& operator=(const CPeakPrediction& p);

  //Functions
  void AddCharge(int charge);         //Add charge state to the vector
  void Clear();                       //Clear all information
  void EraseCharge(int index);        //Remove charge state
  int GetCharge(int index);           //Get charge state from index
  float GetIntensity();               //Return intensity
  double GetMZ();                     //Return m/z
  void SetIntensity(float newIntens); //Set intensity
  void SetMZ(double newMZ);           //Set m/z
  int Size();                         //Get number of charge state predictions
};


/* ******************************************************* */
/*                  CPeptideVariant                        */
/* ******************************************************* */
class CPeptideVariant {
 private:
  Peak_T *match;            //Array of mz and intensities that match the peak list at index n
  int matchSize;
  Peak_T *mismatch;         //Array of mz and intensities that do not match the peak list at index n
  int mismatchSize;

  std::vector<Peak_T> *extra;    //Temporary storage of mismatch data until mismatch array size is determined

  int charge;               //Charge state of variant
  char formula[64];         //Molecular formula (averagine-based)
  double monoMass;          //monoisotopic mass
  double distArea;          //distribution area (sum of centroided peaks)
  CHardklorVariant variant; //Variant settings
 protected:
 public:
  //Constructors & Destructors
  CPeptideVariant();
  CPeptideVariant(const CPeptideVariant& p);
  ~CPeptideVariant();

  //Operator overloads
  CPeptideVariant& operator=(const CPeptideVariant& p);

  //Functions
  void AddExtra(double& mz, float& intensity);
  void AddMatch(int& index, double& mz, float& intensity);
  void AddMismatch(int& index, double& mz, float& intensity);
  void Clear();
  void DeleteExtra();
  double& GetArea();
  int& GetCharge();
  Peak_T& GetExtra(int& index);
  char* GetFormula();
  CHardklorVariant& GetHKVariant();
  Peak_T& GetMatch(int& index);
  Peak_T& GetMismatch(int& index);
  double& GetMonoMass();
  void SetArea(double& d);
  void SetCharge(int ch);
  void SetFormula(char form[64]);
  void SetMatchSize(int& num);
  void SetMismatchSize(int& num);
  void SetMonoMass(double mass);
  void SetVariant(CHardklorVariant& v);
  int SizeExtra();

};


/* ******************************************************* */
/*                CPeptidePrediction                       */
/* ******************************************************* */
class CPeptidePrediction {
 private:
   std::vector<CPeptideVariant> *variantList;   //vector of all variants for this peptide
  double mz;                              //mz peak from which variants are derived
  float intensity;                        //intensity of mz peak
  int bestVar;                            //Index to variant with best match/mismatch score
  int maxPeakIndex;
 protected:
 public:
  //Make variants here?
  /* Pass vector of variants (CHardklorVariant) and make Mercury distributions */
  /* Or calculate the variant externally and add it. */

  //Constructors & Destructors
  CPeptidePrediction();
  CPeptidePrediction(const CPeptidePrediction& p);
  ~CPeptidePrediction();

  //Overloaded operators
  CPeptidePrediction& operator=(const CPeptidePrediction& p);

  //Functions
  void AddVariant(CPeptideVariant& var);
  void Clear();
  float GetIntensity();
  int GetMaxPeakIndex();
  double GetMZ();
  CPeptideVariant& GetVariant(int& index);
  void SetIntensity(float intens);
  void SetMaxPeakIndex(int index);
  void SetMZ(double mass);
  int VariantListSize();

  float& GetIntensityAt(int& variant, int& mz);
};

#endif

