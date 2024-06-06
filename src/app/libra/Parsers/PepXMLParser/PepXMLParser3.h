#ifndef _PEPXMLPARSER3_H
#define _PEPXMLPARSER3_H

#include "expat.h"
#include "PepXMLStructs.h"
#include "CPepXMLAnalysis.h"
#include "CPepXMLSearch.h"
#include "CPepXMLSpectrum.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <string.h>

#define XMLCLASS		
#ifndef XML_STATIC
#define XML_STATIC	// to statically link the expat libraries
#endif

using namespace std;

enum eAnalysisState {
  asNone=0,
  asPeptideProphet,
  asInterProphet,
  asPTMProphet,
  asPTMProphetResult
};

class PepXMLParser3 {

public:
	PepXMLParser3();
	~PepXMLParser3();

	CPepXMLSpectrum& operator[ ](const size_t& i);

	void characters(const XML_Char *s, int len);
	void endElement(const XML_Char *el);
	void startElement(const XML_Char *el, const XML_Char **attr);

  //Assuming functions are usually index-centric (i.e. will be given a peptide array index)

  CPepXMLAnalysis getAnalysis           (size_t index);
  size_t          getAnalysisCount      ();
  string          getFile               (int pepIndex);
  int             getFileCount          ();
  void            getFileFromList       (int index, char* str);
  PepXMLXL        getLinker             (size_t index, char rank=1);
  PepXMLXL        getLinkerFromPos      (size_t index, size_t nth=0);
  bool            getLinkSites          (size_t index, char& a, char& b, char rank=1);
  bool            getLinkSitesFromPos   (size_t index, char& a, char& b, size_t nth=0);
  char            getLinkType           (size_t index, char rank=1);
  string          getPeptide            (size_t index, bool mod=false, char rank=1, bool link=false);
  string          getPeptideFromPos     (size_t index, bool mod = false, size_t nth=0, bool link = false);
  PepXMLMod       getPeptideMod         (size_t pepIndex, size_t modIndex, char rank=1, bool link=false);
  PepXMLMod       getPeptideModFromPos  (size_t pepIndex, size_t modIndex, size_t nth=0, bool link = false);
  size_t          getPeptideModCount    (size_t index, char rank=1, bool link=false);
  size_t          getPeptideModCountFromPos (size_t index, size_t nth=0, bool link = false);
  PepXMLMod       getPeptideModFromList (int modIndex);
  double          getProbability        (double err);
  string          getProtein            (size_t index, size_t protIndex, char rank=1, bool link=false);
  string          getProteinFromPos     (size_t index, size_t protIndex, size_t nth=0, bool link = false);
  string          getProteinDB          (size_t index, size_t protIndex, char rank = 1, bool link = false);
  string          getProteinDesc        (size_t index, size_t protIndex, char rank = 1, bool link = false);
  string          getProteinDescFromPos (size_t index, size_t protIndex, size_t nth=0, bool link = false);
  string          getProteinFromList    (int protIndex);
  CPepXMLPSM      getPSM                (size_t index, char rank = 1, bool link = false);
  CPepXMLPSM      getPSMFromPos         (size_t index, size_t nth=0, bool link=false);
  char            getScoreIndex         (string label);
  string          getScoreLabel         (char scoreIndex);
  CPepXMLSearch   getSearchParams       (size_t pepIndex);

  bool    hasIProphet     ();
  bool    hasPepProphet   ();

	bool    readFile        (const char* fileName, double probFilter=-1,double iProbFilter=-1);

	size_t  size            ();
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

  bool                  bIProphet;    //Indicates data has been iPropheted
  bool                  bPepProphet;  //Indicates data has been PeptidePropheted

  CPepXMLAnalysis       analysis;
  CPepXMLSearch         search;
	CPepXMLSpectrum 			spec;
	XML_Parser						parser;
  CPepXMLPeptide        peptide;
  CPepXMLPSM            psm;

  CPepXMLProphetResult*     prophetResult;
  CPepXMLPTMProphetResult*  ptmProphetResult;

  char                  currentDBID;
  int                   currentFileID;
  char                  currentSearchSummary;

  double                probabilityFilter;
  double                iProbabilityFilter;

  eAnalysisState        analysisState;
  
  //Internal, global lists
  vector<CPepXMLAnalysis> vAnalysis;
  vector<string>          vDBs;
  vector<PepXMLError>     vError;
  vector<string>          vFiles;
  vector<PepXMLMod>       vMods;
  vector<PepXMLProtein>   vProteins;
  vector<string>          vScores;
  vector<CPepXMLSearch>   vSearch;
  vector<CPepXMLSpectrum> vSpectra;
  vector<PepXMLXL>        vXL;

  vector<PepXMLProtTable> vProtTable;

  bool    addMod        (char aa, double mass, double massDiff=0);
  double  calcMonoMass  (const char *seq, bool water=true);
  char    findMod       (char aa, double mass);
  size_t  findProtein   (string& s, char index, string desc="");
  char    findScore     (string& s);
  char    findXL        (string& s, double mass);

  static bool compareProt(const PepXMLProtTable& a, const PepXMLProtTable& b);

};

#endif 