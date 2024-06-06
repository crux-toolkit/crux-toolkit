#ifndef _PEPXMLPARSER_H
#define _PEPXMLPARSER_H

#include "expat.h"
#include <cstring>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

#define XMLCLASS		
#ifndef XML_STATIC
#define XML_STATIC	// to statically link the expat libraries
#endif

using namespace std;

typedef struct PepXMLEntry{
  char prevAA;
  char nextAA;
	int charge;
  int fileID;
	int scanNum;
	float RTime;
	double expect;
	double iProphetProbability;
	double monoMass;
	double precursorMonoMass;
	double probability;
  double xcorr;
  string label;
  string xLabel;
  string modifiedPeptide;  //is this still used? maybe peptide has the mods included already?
  string modifiedPeptidePlus;
	string peptide;
	string protein;
} PepXMLEntry;

typedef struct PepXMLMod{
  double mass;
  int pos;
} PepXMLMod;

typedef struct PepXMLError{
  double error;
  double prob;
} PepXMLError;

class PepXMLParser {

public:
	PepXMLParser();
	~PepXMLParser();

	PepXMLEntry& operator[ ](const unsigned int& i);

	void characters(const XML_Char *s, int len);
	void endElement(const XML_Char *el);
	void startElement(const XML_Char *el, const XML_Char **attr);

  bool    addMod          (char aa, double massDiff, int pos);
  double  calcMonoMass    (char *seq, bool water=true);
  void    getFile         (int index, char* str);
  int     getFileCount    ();
  void    getFileFromList (int index, char* str);
  bool    getIprophet     ();
  double  getProbability  (double err);
	int     parse           (const char* fileName);
  void    setMassCalc     (bool b);
	int     size            ();
  void    sortModPep      ();
  void    sortScanNum     ();

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

  double                aaMass[128];
  bool                  bIprophet;
  bool                  bMassCalc;  //indicates true mass should be recalculated (do not trust pepXML).
  bool                  bHit;       //indicates a rank=1 hit was already added for the spectrum.
	PepXMLEntry						peptide;
	XML_Parser						parser;
	int										rank;
  int                   currentFileID;
  vector<PepXMLError>   vError;
	vector<PepXMLEntry>		vPeptides;
	//vector<int>						vStaticMods;
	//vector<int>						vVarMods;
  unsigned int          nTerminalMod;
  unsigned int          cTerminalMod;
  vector<string>        vFiles;
  vector<PepXMLMod>     aaModList;

  static int compareModPep  (const void *p1, const void *p2); 
  static int compareScanNum (const void *p1, const void *p2); 

};

#endif 
