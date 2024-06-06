/*
Program       : QuanticProteinRatio
Author        : David Shteynberg
Date          : 11.03.2020


Copyright (C) 2020 David Shteynberg

Based on ASAPRatioProteinRatio
Copyright (C) 2003 Andrew Keller

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#ifndef QUANTIC_PROT_RAT_H
#define QUANTIC_PROT_RAT_H


#include "QuanticMSMSSeqData.h"
#include "Common/Array.h"
#include "Common/constants.h"
//#include "Quantitation/ASAPRatio/ASAPRatioPeptideParser/ASAPRatioPeptideParser.h"
// #include "Quantitation/ASAPRatio/ASAPRatio_structs.h"

// structure of protein 
typedef struct {
  int indx_; // 0 new, 1 analyzed, 2 verified, -1 invalid STATUS
  vector<double> quantmeans_; //
  vector<double> quantstdvs_; // 

  //vector<int> dataIndx_; // count data?
  bool dataInc_; // count data?
  string pep_; // peptide sequence
  string specName_;

  long scan_;
  int chg_;
  double weight_;
  int msms_run_idx_;
  unsigned bofIndx_;
  
} psmQuantStrct;

// structure of protein 
typedef struct {
  int indx_; // 0 new, 1 analyzed, 2 verified, -1 invalid STATUS
  vector<double> quantmeans_; //
  vector<double> quantstdvs_; // 

  vector<psmQuantStrct*> psms_;
  bool dataInc_; // count this data? one for each psm
  

  string pep_; // peptide sequence


  int chg_;
  double weight_;

} pepQuantStrct;


// structure for protein quant
typedef struct {
  int indx_; // 0 new, 1 analyzed, 2 verified, -1 invalid STATUS
  vector<double> quantmeans_; //
  vector<double> quantstdvs_; //
  

  vector<pepQuantStrct*> peps_; // unique sequences data
  bool dataInc_; // count this data? one for each pep
  
} proQuantStrct;


class QuanticProteinRatio {
  friend class QuanticGroupPeptideParser;
 public:

  QuanticProteinRatio(double minprob, double minwt);
  ~QuanticProteinRatio();
  void setRunInfo();
  void enter(const char* mod_pep,const ModificationInfo* modinfo, QuanticPSMData  *data, int index, int xml_index, double wt, double prob, int msms_run_idx, int chg);
  void computeProteinQuants();


 protected:

  char** getStrSects(int *sectNum, char *string, char sep);
  void freeMtrx(void **mtrx, int size);
  void getRidOfSpace(char *string);
  void sortMSMSSeqData();
  
  
  //void initProDataStrct();
  
  void initProQuantStrct(); 
  
  //void initSeqDataStrct(seqDataStrct *seqData, int start_index, int end_index, double mnProb);

  void initPepQuantStrct(pepQuantStrct *seqData, int start_index, int end_index, double mnProb);


  //  void getDataStrctRatio(dataStrct *data);
  //void getPSMQuantStrct(psmQuantStrct *data);
  // void getSeqDataStrctRatio(seqDataStrct *data);
   void getSeqQuantStrct(pepQuantStrct *data);
  //void getProDataStrct(proDataStrct *data);
  void getProQuantStrct(proQuantStrct *data);

  proQuantStrct* getProQuantStrct();
  void getDataQuants(vector<double>* quants, vector<double> *stdvs,
		     vector<vector<double>*>* dataQuant,
		     vector<vector<double>*>* dataStdev,
		     vector<double> *dataWghs);

  //void getDataRatio(double *ratio, double *error, double *inv_ratio, double *inv_error, double confL, 
  //		    double *data, double *dataErrs, double *inv_data, double *inv_dataErrs, 
  //		    double *dataWghs, int *dataIndx, int dataSize, int testType);

  void findMeanAndStdDevWeight(double *mean, double *error,
			       double *data, double *inv_mean, double *inv_error,
			       double *inv_data, double *weight, int size);


  void DixonTest(double *data, int *outliers, int size);
  double PadeApprx(double x, double *xa, double *ya, int size);
  void searchPepDataStrct(pepDataStrct *pepData, int pepIndx, int xmlIndx);

  vector<QuanticMSMSSeqData*>* runseqdata_;

  
  int prtnAANum_;
  Boolean set_;
  //char* basename_;
  
  //proDataStrct* pro_;

  proQuantStrct* pro_;

  int start_;
 
  double min_probability_;
  double min_weight_;

};

#endif
