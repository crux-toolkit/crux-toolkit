/*
Program       : ASAPRatioProteinRatio
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPProteinRatio.h 8163 2020-06-02 09:10:54Z real_procopio $

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

#ifndef ASAP_PROT_RAT_H
#define ASAP_PROT_RAT_H


#include "ASAPRatioMSMSSeqData.h"
#include "Common/Array.h"
#include "Common/constants.h"
#include "Quantitation/ASAPRatio/ASAPRatioPeptideParser/ASAPRatioPeptideParser.h"
#include "Quantitation/ASAPRatio/ASAP_structs.h"

//int pairStrctCmp(void const *a, void const *b);
int MSMSSeqDataCmp(void const *a, void const *b);

class ASAPProteinRatio {

 public:

  ASAPProteinRatio(double minprob, double minwt);
  ~ASAPProteinRatio();
  void setRunInfo(/* char* basename, */ char* lightlabels, char* heavylabels);
#ifdef USE_STD_MODS
  void enter(const char *peptide, double lightmass, const ModificationInfo* modinfo, const char *quant_labels, 
	     const pepDataStrct &data, int index, int xml_index, double wt, double prob, int msms_run_idx);
#else
  void enter(const char* peptide, const pepDataStrct &data, int index, int xml_index, double wt, double prob, int msms_run_idx);
#endif
  void computeProteinRatio();
  proDataStrct* getProDataStruct(); 
  pairStrct* collectPrtnAAStrct(int *prtnAANum, char *lightString, char *heavyString);
  char* getLightSequence(const char* sequence, pairStrct *prtnAAs, int prtnAANum);


 protected:

  char** getStrSects(int *sectNum, char *string, char sep);
  void freeMtrx(void **mtrx, int size);
  void getRidOfSpace(char *string);
  ASAPRatioMSMSSeqData** getSortedMSMSSeqData();
  void sortMSMSSeqData();
  void initProDataStrct();
  void initSeqDataStrct(seqDataStrct *seqData, int start_index, int end_index, double mnProb);
  void getDataStrctRatio(dataStrct *data);
  void getSeqDataStrctRatio(seqDataStrct *data);
  void getProDataStrct(proDataStrct *data);
  void getDataRatio(double *ratio, double *error, double *inv_ratio, double *inv_error, double confL, 
		    double *data, double *dataErrs, double *inv_data, double *inv_dataErrs, 
		    double *dataWghs, int *dataIndx, int dataSize, int testType);

  void findMeanAndStdDevWeight(double *mean, double *error,
			       double *data, double *inv_mean, double *inv_error,
			       double *inv_data, double *weight, int size);


  void DixonTest(double *data, int *outliers, int size);
  double PadeApprx(double x, double *xa, double *ya, int size);
  void searchPepDataStrct(pepDataStrct *pepData, int pepIndx, int xmlIndx);

  Array<ASAPRatioMSMSSeqData*>* runseqdata_;
  pairStrct* label_partners_;
  int prtnAANum_;
  Boolean set_;
  //char* basename_;
  proDataStrct* pro_;
  int start_;
  char* lightlabels_;
  char* heavylabels_;
  double min_probability_;
  double min_weight_;

};

#endif
