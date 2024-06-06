/*
Program       : ASAPCGIDisplayParser
Authors       : Andrew Keller <akeller@systemsbiology.org>
                Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPCGIDisplayParser.h 7996 2019-12-25 00:16:42Z real_procopio $

Displays ASAPRatio protein information from ProteinProphet XML

Copyright (C) 2003 Andrew Keller, Xiao-jun Li

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

#ifndef ASAP_CGI_DISP_PARSER_H
#define ASAP_CGI_DISP_PARSER_H


#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/ASAPRatio/Ratio.h"
#include "Quantitation/ASAPRatio/Normalization/Normalization.h"
#include "Quantitation/ASAPRatio/ASAP_structs.h"

class ASAPCGIDisplayParser : public Parser {

 public:

  // given protein, goes into xml and creates a proData structure for displaying html
  ASAPCGIDisplayParser(const char* xmlfile, const char* protein);
  ASAPCGIDisplayParser(const char* xmlfile, const char* protein, unsigned int index);
  ASAPCGIDisplayParser(const char* xmlfile, unsigned int index);
  proDataStrct* getProDataStrct();
  Normalization* getNormalized();
  char* getPvalLink();
  char** getInputFiles();
  Array<char*>* getInputFilesArray();
  const char* getName(void) {return this->protein_;};

 protected:

  void parse(const char* xmlfile);
  void setFilter(Tag* tag);
  void freeProDataStrct(proDataStrct* data);
  proDataStrct* getProDataStruct(Array<Tag*>* tags);
  int setProSeqStruct(proDataStrct* data, int seq_ind, Array<Tag*>* tags, int tag_ind);
  int setProPeakStruct(seqDataStrct* seq, int peak_ind, Array<Tag*>* tags, int tag_ind);
  int setProDtaStruct(dataStrct* peak, int dta_ind, Array<Tag*>* tags, int tag_ind);

  Boolean heavy2light_;
  proDataStrct* pro_ratio_;
  char* protein_;
  unsigned int index_;
  char** inputfiles_;
  Array<char*>* inputfiles_array_;
  Normalization* norm_;
  char* pval_link_;
};

#endif
