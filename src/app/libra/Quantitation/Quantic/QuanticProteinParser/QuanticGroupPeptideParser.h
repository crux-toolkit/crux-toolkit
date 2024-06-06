/*
Program       : QuanticGroupPeptideParser
Author        : David Shteynberg
Date          : 11.07.20
SVN info      : $Id: QuanticGroupPeptideParser.h 

Copyright (C) 2020 David Shteynberg

Based on ASAPRatioGroupPeptideParser
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

#ifndef QUANTIC_GROUP_PEP_PARSER_H
#define QUANTIC_GROUP_PEP_PARSER_H


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Quantitation/Option.h"
#include "Parsers/mzParser/mzParser.h"
#include "QuanticProteinRatio.h"
#include "Common/constants.h"
#include "Common/ModificationInfo/ModificationInfo.h"

using namespace std;

class UniqPeptide {
public:
  UniqPeptide() {
    peptide_ = "";
    weight_ = 0.;
    probability_ = 0.;
  }

  UniqPeptide(string p, double w, double pr, int charge) {
    peptide_ = p;
    weight_ = w;
    probability_ = pr;
    charge_ = charge;
    modpeps_ = new map<string, bool>();
  }

  void insertModPep(string p) {
    modpeps_->insert(make_pair(p, true));
  }

  
  bool findModPep(string p) {
    map<string, bool>::iterator it = modpeps_->find(p);

    if (it == modpeps_->end()){
      return false;
    }

    return it->second;    
  }

  
  
  ~UniqPeptide() {
    modpeps_->clear();
    delete modpeps_;
  }

  
  map<string,bool> * modpeps_;
  int charge_;
  string peptide_;
  double weight_;
  double probability_;
};

class QuanticGroupPeptideParser : public Parser {

 public:

  QuanticGroupPeptideParser(const char* pepxmlfile, Array<UniqPeptide*>* peptides, double minprob, double minwt);
  QuanticGroupPeptideParser(Array<const char*> *pepxmlfiles, Array<UniqPeptide*>* peptides, double minprob, double minwt);
  QuanticGroupPeptideParser(Array<const char*>* pepxmlfiles, double minprob, double minwt, bool ch = true);
  proQuantStrct* getProQuantStruct(Array<UniqPeptide*>* peptides);

  ~QuanticGroupPeptideParser();

  void setUseCharge(bool use);
  
  void setFilter(Tag* tag);
  double getRatioSum();
  double getRatioSquareSum();
  int getRatioNum();
  void parseTuple(vector<double>* result, string input);
  RatioStruct getRatio();
  void clearRatio();
  void setPepDataStruct(Array<Tag*>* tags, string spec_name,
			long scan, int precursor_charge);
  proQuantStrct* getProQuantStruct();

 protected:

  void parse(const char* xmlfile);
  Boolean peptideListMember(const char* pep, double* wt);
  Boolean possiblePeptideListMember(const char* data);

  struct psm {
    std::string sequence;
    double lightmass;
    ModificationInfo* modinfo;
    std::string quant_labels;
    QuanticPSMData* data;
    //psmQuantStrct data;
    int charge;
    int index;
    int xml_index;
    double wt;
    double prob;
    int msms_run_idx;
  };

  Array<UniqPeptide*>* peptides_;
  map <string, vector<psm> > all_peptides_map_; // used when parsing the entire pepXML file(s) only once

  double min_probability_;
  double min_weight_;


  Array<const char*> *pepxmlfiles_;
  Boolean single_input_;

  QuanticProteinRatio* ratio_;

  QuanticPSMData* data_;
  
  //  psmDataStrct data_;

  bool use_charge_;

  Boolean parse_all_;
};


#endif
