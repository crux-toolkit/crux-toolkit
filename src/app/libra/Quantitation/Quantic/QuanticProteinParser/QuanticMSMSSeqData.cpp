/*
Program       : QuanticMSMSSeqData
Author        : David Shteynberg <dshteynb@systemsbiology.org>

Date          : 12.02.20
SVN info      : $Id: QuanticMSMSSeqData.cpp 


Copyright (C) 2020 David Shteynberg

Based on ASAPRatioMSMSSeqData
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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include "QuanticMSMSSeqData.h"


QuanticMSMSSeqData::QuanticMSMSSeqData(const char* modseq, const ModificationInfo* modinfo, QuanticPSMData * data, int index, int xml_index, double wt, double prob, int msms_run_idx, int chg) {

  data_ = new vector<QuanticPSMData*>;
  indices_ = new vector<int>;
  data_->push_back(new QuanticPSMData(data));
  indices_->push_back(index);
  modseq_ = new char[strlen(modseq)+1];
  strcpy(modseq_, modseq);
  xml_index_ = xml_index;
  msms_run_idx_ = msms_run_idx;
  weight_ = wt;
  probs_ = new vector<double>;
  probs_->push_back(prob);

  mod_info_ = new ModificationInfo(modinfo);
  chg_ = chg;
}

QuanticMSMSSeqData::~QuanticMSMSSeqData() {

  if(data_ != NULL) {
    data_->clear();
    delete data_;

  }
  if(indices_ != NULL)
    delete indices_;
  if(probs_ != NULL)
    delete probs_;
  if(modseq_ != NULL)
    delete[] modseq_;
  if(mod_info_ != NULL)
    delete mod_info_;
  
}
