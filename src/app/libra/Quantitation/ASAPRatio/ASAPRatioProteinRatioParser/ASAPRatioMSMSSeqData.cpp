/*
Program       : ASAPRatioMSMSSeqData
Author        : Andrew Keller <akeller@systemsbiology.org>
                David Shteynberg <dshteynb@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioMSMSSeqData.cpp 7996 2019-12-25 00:16:42Z real_procopio $

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

#include "ASAPRatioMSMSSeqData.h"

#ifdef USE_STD_MODS
ASAPRatioMSMSSeqData::ASAPRatioMSMSSeqData(const char* lightseq, double lightmass, const ModificationInfo* modinfo, const char* quant_labels, 
					   const pepDataStrct &data, int index, int xml_index, double wt, double prob, int msms_run_idx) {
  light_mass_ = lightmass;
  mod_info_ = new ModificationInfo(modinfo);
  quant_labels_ = new char[strlen(quant_labels)+1];
  strcpy(quant_labels_, quant_labels);
#else
ASAPRatioMSMSSeqData::ASAPRatioMSMSSeqData(const char* lightseq, const pepDataStrct &data, int index, /* char* basename,*/ int xml_index, 
					   double wt, double prob, int msms_run_idx) {
#endif
  data_ = new Array<pepDataStrct>;
  indices_ = new Array<int>;
  data_->insertAtEnd(data);
  indices_->insertAtEnd(index);
  lightsequence_ = new char[strlen(lightseq)+1];
  strcpy(lightsequence_, lightseq);
  xml_index_ = xml_index;
  msms_run_idx_ = msms_run_idx;
  weight_ = wt;
  probs_ = new Array<double>;
  probs_->insertAtEnd(prob);
}

ASAPRatioMSMSSeqData::~ASAPRatioMSMSSeqData() {
  if(data_ != NULL)
    delete data_;
  if(indices_ != NULL)
    delete indices_;
  if(probs_ != NULL)
    delete probs_;
  if(lightsequence_ != NULL)
    delete[] lightsequence_;

#ifdef USE_STD_MODS
  if(mod_info_ != NULL)
    delete mod_info_;
  if(quant_labels_ != NULL)
    delete[] quant_labels_;
#endif

}
