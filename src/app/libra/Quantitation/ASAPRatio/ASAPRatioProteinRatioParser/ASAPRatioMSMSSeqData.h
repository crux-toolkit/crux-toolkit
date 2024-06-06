/*
Program       : ASAPRatioPeptideParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioMSMSSeqData.h 7996 2019-12-25 00:16:42Z real_procopio $

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

#ifndef ASAP_MSMS_SEQ_H
#define ASAP_MSMS_SEQ_H


#include "Common/sysdepend.h"
#include <string.h>

#include "Common/Array.h"
#include "Common/ModificationInfo/ModificationInfo.h"


class ASAPRatioMSMSSeqData {

 public:

#ifdef USE_STD_MODS
  ASAPRatioMSMSSeqData(const char* lightseq, double lightmass, const ModificationInfo* modinfo, const char* quant_labels, 
		       const pepDataStrct &data, int index, int xml_index, double wt, double prob, int msms_run_idx);
#else
  ASAPRatioMSMSSeqData(const char* lightseq, const pepDataStrct &data, int index, /* char* basename,*/ int xml_index, double wt, double prob, int msms_run_idx);
#endif
  ~ASAPRatioMSMSSeqData();

  Array<pepDataStrct>* data_;
  Array<int>* indices_;
  //char* basename_;
  int xml_index_;
  char* lightsequence_;
#ifdef USE_STD_MODS
  double light_mass_;
  ModificationInfo* mod_info_;
  char* quant_labels_;
#endif
  double weight_;
  Array<double>* probs_; // probability_;

  int msms_run_idx_;


 protected:


};


#endif
