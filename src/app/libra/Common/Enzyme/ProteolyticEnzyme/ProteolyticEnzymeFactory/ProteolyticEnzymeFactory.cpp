#include "ProteolyticEnzymeFactory.h"

/*

Program       : EnzymeSpecificity for PeptideProphet
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ProteolyticEnzymeFactory.cpp 8309 2020-12-03 03:17:24Z dshteyn $

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

#define CATALOG_NAMES_INDEX 0
#define CATALOG_INDEPENDENT_INDEX 1
#define CATALOG_DESCRIPTION_INDEX 2
#define CATALOG_SPECIFICITYRULES_INDEX 3

const char *ProteolyticEnzymeFactory::CATALOG_ENZYMES[][4] = 
  {
    /*
      For a more complete "standard" list see:

         http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MS&termId=MS:1001045&termName=cleavage%20agent%20name

      {<name(s)>, <indepedent>, <description>, <specificity rule(s)>}
      <name(s)> : multiple names are separated by "\n"
      <specificity rule(s)> : multiple rules are separated by ";"
                              each rules is made up of 
                              "<cut>:<no cut>:<n term sense>:[<min spacing>]"
                              where <min spacing> = -1 indicate used default
                                                    blank used factory settings
                                                    numeric value to override
    */
      {"trypsin", "true", NULL, "KR|P|false|"}                            // default & original
    , {"stricttrypsin|trypsin/p|trypsinp", "true", NULL, "KR||false|"}             // added 2007-10-12, DDS
    , {"argc|arg-c|arg_c", "true", NULL, "R|P|false|"}                    // added 2007-02-07, WCH
    , {"aspn|asp-n|asp_n", "true", NULL, "D||true|"}                      // original
    , {"chymotrypsin", "true", NULL, "FLWY|P|false|"}                     // original 
    , {"clostripain", "true", NULL, "R|-|false|"}                         // original 
    , {"cnbr", "true", NULL, "M||false|"}                                 // original 
    , {"elastase", "true", NULL, "AGILV|P|false|"}                        // original 
    , {"formicacid|formic_acid", "true", NULL, "D|P|false|"}              // added 2007-02-07, WCH
    , {"gluc|glu_c|v8-de", "true", NULL, "DE|P|false|"}                   // original
    , {"gluc_bicarb|v8-e", "true", NULL, "E|P|false|"}                    // original
    , {"iodosobenzoate", "true", NULL, "W|-|false|"}                      // original
    , {"lysc|lys-c|lys_c", "true", NULL, "K|P|false|"}                    // added 2007-02-07, WCH
    , {"lysc-p|lys-c/p", "true", NULL, "K||false|"}                       // added 2007-02-07, WCH
    , {"lysn|lys_n", "true", NULL, "K||true|"}                            // original
    , {"lysn_promisc", "true", NULL, "KR||true|"}                         // original
    , {"lysarginase", "true", NULL, "KR|P|true|"}                       // added 2019-12-18, DDS
    , {"ralphtrypsin", "true", NULL, "KRST|P|false|"}                     // original
    , {"nonspecific|no enzyme|no cleavage", "true", NULL, ""}             // updated
    , {"pepsina", "true", NULL, "FL|-|false|"}                            // added 2007-02-07, WCH
    , {"proline_endopeptidase", "true", NULL, "P|-|false|"}               // original
    , {"trypsin/chymotrypsin|trypsinchymotrypsin|trypchymo|trypsin_chymotrypsin", "true", NULL, "FYWLKR|P|false|"}
    , {"staph_protease", "true", NULL, "E|-|false|"}                      // original
    , {"tca", "true", NULL, "KR|P|false|,YWFM|P|false|,D||true|-1"}       // original
    , {"trypsin/cnbr|tryp-cnbr", "true", NULL, "KR|P|false|,M|P|false|"}  // original
    , {"trypsin_gluc", "true", NULL, "DEKR|P|false|"}                     // original
    , {"trypsin_k", "true", NULL, "K|P|false|"}                           // original
    , {"trypsin_r", "true", NULL, "R|P|false|"}                           // original
    , {"thermolysin", "true", NULL, "ALIVFM|DE|true|"}                    // added 2007-08-04, WCH, [Bill Vensel]
  };


ProteolyticEnzymeFactory::ProteolyticEnzymeFactory() { 
  min_spacing_ = 1;
  user_enzyme_ = false;
}

ProteolyticEnzymeFactory::ProteolyticEnzymeFactory(int min_spacing) { 
  user_enzyme_ = false;
  if (min_spacing < 0) 
    min_spacing_ = 0;
  else  
    min_spacing_ = min_spacing;
}

void ProteolyticEnzymeFactory::showEnzymesCatalog() {
  for (int i=0; i<(sizeof(ProteolyticEnzymeFactory::CATALOG_ENZYMES)/sizeof(ProteolyticEnzymeFactory::CATALOG_ENZYMES[0])); i++) {
    ProteolyticEnzyme* enzyme = getEnzymeFromCatalog(i, NULL);
    enzyme->write(cout);
    delete enzyme;
  }
}

// factory for producing enzyme digestions
// register all new enzyme digestions here
// please include/remove the name(s) into char *enzymes[] of showEnzymesCatalog()
ProteolyticEnzyme* ProteolyticEnzymeFactory::getProteolyticEnzyme(const char* input_name) {
  ProteolyticEnzyme* enz = NULL;
  // go to lower case
  Boolean isSemiSpecific = False;
  char *name = NULL;
  ProteolyticEnzyme::parseEnzymeName(input_name, &name, isSemiSpecific);
  const char *fidelity = (isSemiSpecific?ProteolyticEnzyme::CONSTANT_FIDELITY_SEMISPECIFIC_:NULL);

  if(name == NULL) {// default
    enz = getEnzymeFromCatalog("trypsin", NULL);
  } else {
    enz = getEnzymeFromCatalog(name, fidelity);
  }

  delete[] name;
  return enz;
}


ProteolyticEnzyme* ProteolyticEnzymeFactory::getProteolyticEnzyme(const string* input_name) {
  if (input_name->find(":") == std::string::npos)
    return getProteolyticEnzyme(input_name->c_str());

  user_enzyme_ = true;
  
  std::istringstream lin(*input_name); 

  int col = 0;
  std::string name="";
  std::string desc=""; //FUTURE USE???
  std::string fidel="";
  std::string indep="";
  std::string specs="";

  while (lin) {
    std::string dat;
    if (!std::getline (lin, dat, ':')) break;
    
    switch (col++) {
    case 0:
      name = dat;
      break;
    case 1:
      fidel = dat;
      break;
    case 2:
      indep = dat;
      break;
    case 3:
      specs = dat;	
      break;
    default:
      break;
      
    }
  }

  return getProteolyticEnzyme(name.c_str(), fidel.c_str(),indep.c_str(), desc.c_str(), specs.c_str());
  
}

ProteolyticEnzyme* ProteolyticEnzymeFactory::getEnzymeFromCatalog(const char* name, const char* fidelity) const {
  ProteolyticEnzyme* enz = NULL;

  // let's look up the entry
  for (int i=0; i<(sizeof(ProteolyticEnzymeFactory::CATALOG_ENZYMES)/sizeof(ProteolyticEnzymeFactory::CATALOG_ENZYMES[0])); i++) {
    if (isEnzymeMatched(name, CATALOG_ENZYMES[i][CATALOG_NAMES_INDEX])) {
      enz = getEnzymeFromCatalog (i, fidelity);
      break;
    }
  }

  return enz;
}

ProteolyticEnzyme* ProteolyticEnzymeFactory::getEnzymeFromCatalog(int index, const char* fidelity) const {
  if (index >= (sizeof(ProteolyticEnzymeFactory::CATALOG_ENZYMES)/sizeof(ProteolyticEnzymeFactory::CATALOG_ENZYMES[0]))) {
    return NULL;
  }

  return getProteolyticEnzyme (CATALOG_ENZYMES[index][CATALOG_NAMES_INDEX], 
			       fidelity, 
			       CATALOG_ENZYMES[index][CATALOG_INDEPENDENT_INDEX],
			       CATALOG_ENZYMES[index][CATALOG_DESCRIPTION_INDEX],
			       CATALOG_ENZYMES[index][CATALOG_SPECIFICITYRULES_INDEX]);
}

ProteolyticEnzyme* ProteolyticEnzymeFactory::getProteolyticEnzyme(const char *name, 
  const char* fidelity, const char* independent, const char *description, 
  const char *specificities) const {

  // {"trypsin", "True", NULL, "KR:P:False:"}
  string standardName = getEnzymeStandardName(name);
  ProteolyticEnzyme* enz = new ProteolyticEnzyme(standardName.c_str(), fidelity, 
    strcmp("true", independent) ? False : True, description);

  // load the specificity rules
  bool endOfRules=false;
  string rules(specificities);
  do {
    //"KR|P|false|,YWFM|P|false|,D||true|-1"

    string rule;
    size_t separatorPos = rules.find(',');
    endOfRules = (separatorPos == rules.npos);

    if (endOfRules) {
      rule = rules;
    } else {
      rule = rules.substr(0, separatorPos);
      rules = rules.substr(separatorPos+1);
    }

    string cut;
    separatorPos = rule.find('|');
    if (separatorPos==rule.npos) {
      cut = rule;
    } else {
      cut = rule.substr(0, separatorPos);
      rule = rule.substr(separatorPos+1);
    }

    string nocut;
    separatorPos = rule.find('|');
    if (separatorPos==rule.npos) {
      nocut = rule;
    } else {
      nocut = rule.substr(0, separatorPos);
      rule = rule.substr(separatorPos+1);
    }

    string value;
    separatorPos = rule.find('|');
    if (separatorPos==rule.npos) {
      value = rule;
    } else {
      value = rule.substr(0, separatorPos);
      rule = rule.substr(separatorPos+1);
    }
    Boolean n_sense = strcmp("true", value.c_str()) ? False : True;

    separatorPos = rule.find('|');
    if (separatorPos==rule.npos) {
      value = rule;
    } else {
      value = rule.substr(0, separatorPos);
      rule = rule.substr(separatorPos+1);
    }

    if (value == "-1") {
      enz->enterSpecificity(cut.c_str(), nocut.c_str(), n_sense);
    } else {
      if (value == "") {
        enz->enterSpecificity(cut.c_str(), nocut.c_str(), n_sense, min_spacing_);
      } else {
        enz->enterSpecificity(cut.c_str(), nocut.c_str(), n_sense, atoi(value.c_str()));
      }
    }
  } while (!endOfRules);
  enz->fixSpecificity();

  return enz;
}

string ProteolyticEnzymeFactory::getEnzymeStandardName(const char *names) const {
  string standardName(names);
  size_t separatorPos = standardName.find('|');
  if (separatorPos != standardName.npos) {
    standardName = standardName.substr(0, separatorPos);
  }

  return standardName;
}

bool ProteolyticEnzymeFactory::isEnzymeMatched(const char *findName, const char *names) const {
  bool nameFound=false;
  bool endOfList=false;
  string nameList(names);
  do {
    string singleName;
    size_t separatorPos = nameList.find('|');
    endOfList = (separatorPos == nameList.npos);

    if (endOfList) {
      singleName = nameList;
    } else {
      singleName = nameList.substr(0, separatorPos);
      nameList = nameList.substr(separatorPos+1);
    }
    nameFound = (singleName==findName);
  } while (!nameFound && !endOfList);

  return nameFound;
}
