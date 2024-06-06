// TandemParamsParser.h
//     Handles parsing of X!Tandem input parameters from a
//     results file, into a set of name-value maps.

#include <iostream>
#include <string>
#include <map>

using namespace std;

typedef map<string,string,less<string> > xMap;

#include "Parsers/mzParser/mzParser.h"
#include "Parsers/Algorithm2XML/saxtandemhandler.h"


class TandemParamsParser : public TANSAXHandler {
 public:
  TandemParamsParser() {
    initEnzymes();

    pParamsCur = NULL;
  }

  virtual void startElement(const XML_Char *el, const XML_Char **attr);
  virtual void endElement(const XML_Char *el);
  virtual void characters(const XML_Char *s, int len);

  bool writePepXML(ostream& out);

  string getUsedParam(const char* name)
  {return getParam(name, used); }
  string getPerfParam(const char* name)
  {return getParam(name, perf); }

  string getSearchEnzyme();
  string getSampleEnzyme();

  void usageEnzymes();

 protected:
  void writePepXML(ostream& out, xMap& params);

  static string getParam(const char* name, xMap& params) {
    string s;
    xMap::iterator it = params.find(name);
    if (it != params.end())
      s = it->second;
    return s;
  }

 public:
  xMap used;
  xMap unused;
  xMap perf;

 protected:
  xMap* pParamsCur;
  string paramCurName;
  string paramCurValue;

 protected:
  static xMap enzymes;
  static void initEnzymes();
  static void addEnzyme(const char* name, const char* site);
  static string normalizeCleavageSite(string site);
  static string normalizeCleavagePart(string part);
};
