/*
Copyright 2017, Michael R. Hoopmann, Institute for Systems Biology
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef _PROTXMLPARSER2_H
#define _PROTXMLPARSER2_H

#include "cpxProteinGroup.h"
#include "expat.h"
#include <iostream>
#include <ctime>
#include <string>
#include <cstring>
#include <vector>

#define XMLCLASS		
#ifndef XML_STATIC
#define XML_STATIC	// to statically link the expat libraries
#endif

enum protXMLElement{
  indistinguishable_peptide,
  peptide
};

typedef struct sSTPSummary{
  std::string version;
  bool bDegen;
  double FDR;
  double probability;
  double sampleLoad;
  double tolerance;
  sSTPSummary(){
    version="0";
    bDegen=false;
    FDR=0;
    probability=0;
    sampleLoad=0;
    tolerance=0;
  }
}sSTPSummary;

class ProtXMLParser2 {
public:

  //Constructor & Destructor
  ProtXMLParser2();
  ~ProtXMLParser2();

  //operators
  cpxProteinGroup& operator[](const size_t index);

  //Data members
  std::vector<cpxProteinGroup> proteinGroup;

  //Functions
  std::string getAlgorithm();
  std::string getSourceFile(size_t index);
  size_t getSourceFileCount();
  sSTPSummary getStPeterSummary();
  std::string getVersion();
  bool   hasStPeter();
  bool   readFile(const char* fn);
  size_t size();
  bool   writeFile(const char* fn);

  //Functions for XML Parsing
  void characters(const XML_Char *s, int len);
  void endElement(const XML_Char *el);
  void startElement(const XML_Char *el, const XML_Char **attr);

protected:
  std::string                 algorithmName;
  bool                        bFixMods;
  bool                        killRead;
  XML_Parser				          parser;
  std::vector<std::string>*   sourceFiles;
  std::string                 version;
  std::vector<protXMLElement> elements;
  sSTPSummary                 STPSummary;


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
};

#endif