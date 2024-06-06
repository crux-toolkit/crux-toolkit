#ifndef _NEOPROTXMLPARSER_H
#define _NEOPROTXMLPARSER_H

#include "expat.h"
#include "NeoProtXMLStructs.h"

#include "CnprProteinGroup.h"
#include "CnprProteinSummary.h"

#include <iostream>
#include <vector>

//#define XMLCLASS		
//#ifndef XML_STATIC
//#define XML_STATIC	// to statically link the expat libraries
//#endif

#define NPR_VERSION "1.0.0"
#define NPR_DATE "2022 MAY 19"


class NeoProtXMLParser {
public:
  NeoProtXMLParser();
  ~NeoProtXMLParser();

  CnprProteinSummary protein_summary;

  bool read(const char* fn);
  std::string versionNeo(); //returns version information
  bool write(const char* fn, bool tabs=true);

  //Functions for XML Parsing
  void characters(const XML_Char *s, int len);
  void endElement(const XML_Char *el);
  void startElement(const XML_Char *el, const XML_Char **attr);

protected:
  bool                killRead;
  XML_Parser				  parser;
  std::vector<protXMLElement> activeEl;
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

  std::string elements[PROTXML_NUM_ELEMENTS];

  void init();

};

#endif