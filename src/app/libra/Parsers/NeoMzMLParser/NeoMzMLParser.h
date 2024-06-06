#ifndef _NEOMZMLPARSER_H
#define _NEOMZMLPARSER_H

#include "expat.h"
#include "NeoMzMLStructs.h"
#include "CnmzIndexedmzML.h"
#include "CnmzMzML.h"

#include <iostream>
#include <vector>
#include <stdio.h>

#define XMLCLASS		
#ifndef XML_STATIC
#define XML_STATIC	// to statically link the expat libraries
#endif

class NeoMzMLParser {
public:
  NeoMzMLParser();
  ~NeoMzMLParser();

  CnmzIndexedmzML indexedmzML;
  CnmzMzML mzML;

  CnmzSpectrum* nextSpectrum(bool buffered=false);
  bool read(const char* fn, bool iterative=false);
  bool write(const char* fn, bool tabs = true, bool iterative = false);
  //bool writeChromat()
  mzMLWriteState writeNextState();
  bool writeSpectrum(CnmzSpectrum* spec);
  mzMLWriteState writeState();

  //Functions for XML Parsing
  void characters(const XML_Char *s, int len);
  void endElement(const XML_Char *el);
  void startElement(const XML_Char *el, const XML_Char **attr);

protected:
  bool                killRead;
  XML_Parser				  parser;
  std::vector<mzMLElement> activeEl;
  int version;  //1=1.1, 2=1.2, etc.


  //Functions for XML Parsing
  inline const char* getAttrValue(const char* name, const XML_Char **attr) {
    for (int i = 0; attr[i]; i += 2) {
      if (isAttr(name, attr[i])) return attr[i + 1];
    }
    return "";
  }
  inline bool isAttr(const char *n1, const XML_Char *n2) { return (strcmp(n1, n2) == 0); }
  inline bool isElement(const char *n1, const XML_Char *n2)	{ return (strcmp(n1, n2) == 0); }

private:

  size_t scanIndex;
  size_t spectrumIndexID;
  bool bDoParseIndex;
  bool bIndexed;
  bool bIterative;
  bool bParseAbort;
  bool bParseIndexList;
  mzMLElement elSkip;
  std::string elements[MZML_NUM_ELEMENTS];
  std::string fileName;
  FILE* xml_fptr;


  //For writing
  bool bWIterative;
  bool bWTabs;
  int WSpecCount;
  mzMLWriteState WState;
  FILE* Wptr;
  CnmzIndexList WindexList;

  void init();
  void parse(f_off_nmz offset=0);
  void parseIndex();
  void processCvParam(CnmzCvParam& c);
  void processUserParam(CnmzUserParam& c);
  void xmlError(XML_Error error);

};

#endif
