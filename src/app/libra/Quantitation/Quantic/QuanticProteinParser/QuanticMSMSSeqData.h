/*
Program       : QuanticMSMSSeqData
Author        : David Shteynberg
Date          : 11.03.20
SVN info      :

Copyright (C) 2020 David Shteynberg

Based on ASAPRatioMSMSSeqData
Copyright (C) 2003 Andrew Keller, 

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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#ifndef QUANTIC_MSMS_SEQ_H
#define QUANTIC_MSMS_SEQ_H


#include "Common/sysdepend.h"
#include <string.h>

#include "Common/Array.h"
#include "Common/ModificationInfo/ModificationInfo.h"

class QuanticPSMData {
public:
  QuanticPSMData() {
    rawquants_ = new vector<double>;
    quantmeans_ = new vector<double>;
    quantstdvs_ = new vector<double>;
  }
  QuanticPSMData(int indx, long scan, int chrg, int cidIndx_, const ModificationInfo* modinfo) {
    rawquants_ = new vector<double>;
    quantmeans_ = new vector<double>;
    quantstdvs_ = new vector<double>;
  }

  QuanticPSMData(const QuanticPSMData* data) {
    rawquants_ = new vector<double>;
    quantmeans_ = new vector<double>;
    quantstdvs_ = new vector<double>;
    for (int i=0; i< data->rawquants_->size(); i++) {
      rawquants_->push_back((*data->rawquants_)[i]);
      quantmeans_->push_back((*data->quantmeans_)[i]);
      quantstdvs_->push_back((*data->quantstdvs_)[i]);
    }
   
    specName_ = data->specName_;
    chg_ = data->chg_;
    scan_ = data->scan_;
    cidIndx_ = data->cidIndx_;
    indx_ = data->indx_;
  }
  
  ~QuanticPSMData() {
    delete rawquants_;
    delete quantmeans_;
    delete quantstdvs_;
  }

  
  void enterQuant( double quant, double mean, double stdev) {
    rawquants_->push_back(quant);
    quantmeans_->push_back(mean);
    quantstdvs_->push_back(stdev);
  }
  

  double getTotalQuant() {
    double rtn = 0;
    for (int i=0; i< rawquants_->size(); i++) {
      rtn += (*rawquants_)[i];
    }
    return rtn;
    
  }
  
  vector<double>* rawquants_;

  vector<double>* quantmeans_;
  vector<double>* quantstdvs_;

  int indx_;
  long scan_;
  int chg_;
  string specName_;
  int cidIndx_;

};

class QuanticMSMSSeqData {

 public:

  QuanticMSMSSeqData(const char* modseq, const ModificationInfo* modinfo, QuanticPSMData *data, int index, int xml_index, double wt, double prob, int msms_run_idx, int chg);

  ~QuanticMSMSSeqData();

  vector<QuanticPSMData*>* data_;
  //  Array<pepDataStrct>* data_;
  vector<int>* indices_;
  //char* basename_;
  int xml_index_;

  char* modseq_;

  ModificationInfo* mod_info_;
  char* quant_labels_;

  double weight_;
  vector<double>* probs_; // probability_;

  int msms_run_idx_;
  int chg_;

 protected:


};


#endif
