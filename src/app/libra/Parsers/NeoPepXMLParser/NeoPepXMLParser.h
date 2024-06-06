#ifndef _NEOPEPXMLPARSER_H
#define _NEOPEPXMLPARSER_H

#include "CnpxMSMSPipelineAnalysis.h"
#include "expat.h"
#include "NeoPepXMLStructs.h"

#include "CnpxUIPipeline.h"
#include "CnpxUIPSM.h"
#include "CnpxUIRunSummary.h"
#include "CnpxUISpectra.h"

#include <iostream>
#include <vector>
#include <stdio.h>

#define XMLCLASS		
#ifndef XML_STATIC
#define XML_STATIC	// to statically link the expat libraries
#endif

#define NPX_VERSION "1.0.5"
#define NPX_DATE "2023 MAY 15"


class NeoPepXMLParser {
public:
  NeoPepXMLParser();
  ~NeoPepXMLParser();

  //User interface
  CnpxUIPipeline uiPipelines;
  CnpxUIRunSummary uiRunSummaries;
  CnpxUISpectra uiSpectra; //shortcuts to set of spectrum_query of the slected msms_peipline_analysis and msms_run_summary

  std::vector<CnpxMSMSPipelineAnalysis> msms_pipeline_analysis;

  //Functions for user interface
  CnpxUIPSM& operator[](const size_t& index);

  void addMSMSPipelineAnalysis(std::string date, std::string summary_xml);
  void setFilterProbability(double probability);
  void setFilterRunSummary(std::string str);
  void setFilterSearchHit(std::string str);
  bool setRunSummaries(const size_t pipeIndex);
  bool setSpectra(const size_t pipeIndex, const size_t runIndex);
  size_t size();
  bool read(const char* fn);
  std::string versionNeo(); //returns version information
  bool write(const char* fn, bool tabs=false);

  //Functions for XML Parsing
  void characters(const XML_Char *s, int len);
  void endElement(const XML_Char *el);
  void startElement(const XML_Char *el, const XML_Char **attr);

protected:
  bool                killRead;
  XML_Parser				  parser;
  std::vector<pepXMLElement> activeEl;
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

  CnpxUIPSM psm;
  size_t sz;
  double probFilter;
  std::string rsFilter;
  std::string shFilter;

  double pProb;
  double iProb;

  std::string elements[PEPXML_NUM_ELEMENTS];
  
  void calcSize();
  void init();
  
};

#endif
