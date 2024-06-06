#ifndef _PROTXMLPARSER_H
#define _PROTXMLPARSER_H

#include "expat.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <cstdio>

#define XMLCLASS		
#ifndef XML_STATIC
#define XML_STATIC	// to statically link the expat libraries
#endif

using namespace std;

typedef struct ProtXMLPeptide{
	bool nonDegenerate;
	bool evidence;
  double initProb;
  double nspProb;
  double calcNeutralPepMass;
  int instances;
  int charge;
  int length;
	string sequence;
  string modSequence;
  int cMod;
  int nMod;
  int* mods;

  ProtXMLPeptide(){
    nonDegenerate=true;
    evidence=false;
    initProb=0;
    nspProb=0;
    sequence="";
    instances=0;
    charge=0;
    length=0;
    calcNeutralPepMass=0;
    cMod=0;
    nMod=0;
    mods = NULL;
  }
  ProtXMLPeptide(const ProtXMLPeptide& p){
    nonDegenerate=p.nonDegenerate;
    evidence=p.evidence;
    sequence=p.sequence;
    initProb=p.initProb;
    nspProb=p.nspProb;
    instances=p.instances;
    modSequence=p.modSequence;
    charge=p.charge;
    length=p.length;
    calcNeutralPepMass=p.calcNeutralPepMass;
    cMod=p.cMod;
    nMod=p.nMod;
    if(length==0) mods=NULL;
    else {
      mods = new int[length];
      for(int i=0;i<length;i++) mods[i]=p.mods[i];
    }
  }
  ~ProtXMLPeptide(){
    if(mods!=NULL) delete [] mods;
  }
  ProtXMLPeptide& operator=(const ProtXMLPeptide& p){
    if(this!=&p){
      nonDegenerate=p.nonDegenerate;
      evidence=p.evidence;
      sequence=p.sequence;
      modSequence=p.modSequence;
      initProb=p.initProb;
      nspProb=p.nspProb;
      instances=p.instances;
      charge=p.charge;
      calcNeutralPepMass=p.calcNeutralPepMass;
      cMod=p.cMod;
      nMod=p.nMod;
      if(mods!=NULL) delete [] mods;
      if(length==0) mods=NULL;
      else {
        mods = new int[length];
        for(int i=0;i<length;i++) mods[i]=p.mods[i];
      }
    }
    return *this;
  }
} ProtXMLPeptide;

typedef struct ProtXMLStPeter{
  double SI;
  double SIn;
  double dSIn;
  int counts;
  double NSAF;
  double ng;
  double ngCounts;
  ProtXMLStPeter(){
    SI = 0;
    SIn = 0;
    dSIn = 0;		   
    counts = 0;
    NSAF = 0;
    ng = 0;
    ngCounts = 0;
  }
  void clear(){
    SI=0;
    SIn=0;
    dSIn=0;		   
    counts=0;
    NSAF=0;
    ng=0;
    ngCounts=0;
  }
} ProtXMLStPeter;

typedef struct ProtXMLEntry{
	double                  probability;
  ProtXMLStPeter          stPeter;
  string                  proteinDescription;
	string                  proteinName;
	int                     groupID;
  int                     length;
	int                     subID;
  float                   coverage;
	vector<string>*         altProteinNames;
	vector<ProtXMLPeptide>* peptides;

	ProtXMLEntry(){
		peptides = new vector<ProtXMLPeptide>;
		altProteinNames = new vector<string>;
    coverage=0;
		probability=0.0;
    stPeter.clear();
    proteinDescription="";
		proteinName="";
		groupID=0;
    length=0;
		subID=0;
	}
	ProtXMLEntry(const ProtXMLEntry& p){
		unsigned int i;
    coverage=p.coverage;
		probability=p.probability;
    stPeter=p.stPeter;
    proteinDescription=p.proteinDescription;
		proteinName=p.proteinName;
		peptides = new vector<ProtXMLPeptide>;
		for(i=0;i<p.peptides->size();i++) peptides->push_back(p.peptides->at(i));
		altProteinNames = new vector<string>;
		for(i=0;i<p.altProteinNames->size();i++) altProteinNames->push_back(p.altProteinNames->at(i));
		groupID=p.groupID;
    length=p.length;
		subID=p.subID;

	}
	~ProtXMLEntry(){
		delete peptides;
		delete altProteinNames;
	}

	ProtXMLEntry& operator=(const ProtXMLEntry& p){
		if (this != &p) {
			unsigned int i;
      coverage=p.coverage;
			probability=p.probability;
      stPeter=p.stPeter;
      proteinDescription=p.proteinDescription;
			proteinName=p.proteinName;
			delete peptides;
			peptides = new vector<ProtXMLPeptide>;
			for(i=0;i<p.peptides->size();i++) peptides->push_back(p.peptides->at(i));
			delete altProteinNames;
			altProteinNames = new vector<string>;
			for(i=0;i<p.altProteinNames->size();i++) altProteinNames->push_back(p.altProteinNames->at(i));
			groupID=p.groupID;
      length=p.length;
			subID=p.subID;
		}
		return *this;
	}

} ProtXMLEntry;

typedef struct ProtXMLError{
  double error;
  double prob;
} ProtXMLError;

class ProtXMLParser {

public:
	ProtXMLParser();
	~ProtXMLParser();

	ProtXMLEntry& operator[ ](const size_t& i);

	void characters(const XML_Char *s, int len);
	void endElement(const XML_Char *el);
	void startElement(const XML_Char *el, const XML_Char **attr);

  string  getDatabase     ();
  string  getPepXML       ();
  double  getProbability  (double err);
	bool    parse           (const char* fileName);
	int     size            ();

	void sortProbabilityRev();

	// -----------------------------------------------------------------------
	//  SAXHandler helper functions
	// -----------------------------------------------------------------------
	inline const char* getAttrValue(const char* name, const XML_Char **attr) {
		for (int i = 0; attr[i]; i += 2) {
			if (isAttr(name, attr[i])) return attr[i + 1];
		}
		return "";
	}
	inline bool isAttr(const char *n1, const XML_Char *n2) {	return (strcmp(n1, n2) == 0); }
	inline bool isElement(const char *n1, const XML_Char *n2)	{	return (strcmp(n1, n2) == 0); }

protected:

	bool									bIndistinguishablePeptide;
  bool                  bIndistinguishableProtein;

  string                database;
  //double                modTable[128];
  XML_Parser						parser;
	ProtXMLPeptide				peptide;
  string                pepXML;
  double								probability;
	string								proteinName;
	ProtXMLEntry					protein;
	int										rank;
  vector<ProtXMLError>  vError;
	vector<ProtXMLEntry>	vProteins;

  //void        createModTable();
  //int         modLookup(char aa, double mass);
  void        parseModPeptide(string& s, ProtXMLPeptide& p);
	static int  compareProbRev(const void *p1,const void *p2);

};

#endif 
