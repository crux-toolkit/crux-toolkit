#ifndef _PEPXMLSTRUCTS_H
#define _PEPXMLSTRUCTS_H

#include <string>

typedef struct PepXMLError {
  double error;
  double prob;
} PepXMLError;

typedef struct PepXMLPepMod {
  char index;
  char pos; //0-based
} PepXMLPepMod;

typedef struct PepXMLPepScore {
  char index;
  double value;
} PepXMLPepScore;

typedef struct PepXMLProtein {
  char DBindex;
  std::string name;
  std::string description;
} PepXMLProtein;

typedef struct PepXMLSearchMod {
  char aa;
  bool fixed;
  bool proteinTerm;
  double mass;
} PepXMLSearchMod;

typedef struct PepXMLXL {
  std::string ID;
  double mass;
} PepXMLXL;

typedef struct PepXMLMod {
  char aa;              //the amino acid, stored as letter in pepXML, returned as position when requested
  double massSearch;    //mass that was in the search parameters
  double massDiff;      //mass difference that was in the search parameters
  double massStd;       //standard representation of the searched mass
  double massDiffStd;   //standard representation of the diffential mass
  std::string label;         //Give it a name
} PepXMLMod;

typedef struct PepXMLProtTable {
  char DBIndex;
  size_t protIndex;
  std::string name;
} PepXMLProtTable;

typedef struct PepXMLParam {
  std::string name;
  std::string value;
} PepXMLParam;

typedef struct PepXMLPTMMod {
  char position;
  double pVal;
  double probability;
  double oscore;
  double mscore;
  double ctermscore;
  double ntermscore;
} PepXMLPTMMod;

#endif