/*
Program       : ASAPRatioGroupPeptideParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioGroupPeptideParser.h 7996 2019-12-25 00:16:42Z real_procopio $

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

#ifndef ASAP_GROUP_PEP_PARSER_H
#define ASAP_GROUP_PEP_PARSER_H


#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/Option.h"
//#include "Validation/MixtureModel/MixtureModel.h"
#include "Parsers/mzParser/mzParser.h"
#include "ASAPProteinRatio.h"
#include "Common/constants.h"
#include "Quantitation/ASAPRatio/UniquePeptide/UniquePeptide.h"
#include "Common/ModificationInfo/ModificationInfo.h"

class ASAPRatioGroupPeptideParser : public Parser {

 public:

  ASAPRatioGroupPeptideParser(const char* pepxmlfile, Array<UniquePeptide*>* peptides, double minprob, double minwt, Boolean heavy2light);
  ASAPRatioGroupPeptideParser(Array<const char*> *pepxmlfiles, Array<UniquePeptide*>* peptides, double minprob, double minwt, Boolean heavy2light);
  ASAPRatioGroupPeptideParser(Array<const char*>* pepxmlfiles, double minprob, double minwt, Boolean heavy2light);
  proDataStrct* getProDataStruct(Array<UniquePeptide*>* peptides);

  ~ASAPRatioGroupPeptideParser();

  void setFilter(Tag* tag);
  double getRatioSum();
  double getRatioSquareSum();
  int getRatioNum();
  RatioStruct getRatio();
  void clearRatio();
  void setPepDataStruct(Array<Tag*>* tags, int elution, long scan, int precursor_charge);
  proDataStrct* getProDataStruct();

 protected:

  void parse(const char* xmlfile);
  Boolean peptideListMember(const char* pep, double* wt);
  Boolean possiblePeptideListMember(const char* data);

  struct psm {
    std::string sequence;
    double lightmass;
    ModificationInfo* modinfo;
    std::string quant_labels;
    pepDataStrct data;
    int index;
    int xml_index;
    double wt;
    double prob;
    int msms_run_idx;
  };

  Array<UniquePeptide*>* peptides_;
  map <string, vector<psm> > all_peptides_map_; // used when parsing the entire pepXML file(s) only once

  double min_probability_;
  double min_weight_;

  double ratio_sum_;
  double ratio_square_sum_;
  double inv_ratio_sum_;
  double inv_ratio_square_sum_;
  int ratio_num_;
  Boolean heavy2light_;

  Array<const char*> *pepxmlfiles_;
  Boolean single_input_;
  ASAPProteinRatio* ratio_;
  pepDataStrct data_;

  Boolean parse_all_;
};


#endif
