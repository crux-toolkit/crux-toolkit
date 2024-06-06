// TandemParamsParser.cxx
//     Handles parsing of X!Tandem input parameters from a
//     results file, into a set of name-value maps.

#include <vector>
#include <algorithm>

#include "Tandem2xml.h"
#include "TandemParamsParser.h"

static const char EL_NOTE[] = "note";
static const char EL_GROUP[] = "group";

void TandemParamsParser::startElement(const XML_Char *el, const XML_Char **attr) {
  if (pParamsCur != NULL) {
    if (isElement(EL_NOTE, el)) {
      paramCurName.assign(getAttrValue("label", attr));
    }
  }
  else if (isElement(EL_GROUP, el)) {
    if (strcmp("parameters", getAttrValue("type", attr)) != 0)
      return;

    const char* label = getAttrValue("label", attr);
    if (strcmp("input parameters", label) == 0)
      pParamsCur = &used;
    else if (strcmp("unused input parameters", label) == 0)
      pParamsCur = &unused;
    else if (strcmp("performance parameters", label) == 0)
      pParamsCur = &perf;
  }
}

void TandemParamsParser::endElement(const XML_Char *el) {
  if (isElement(EL_GROUP, el))
    pParamsCur = NULL;
  else if (isElement(EL_NOTE, el)) {
    if (pParamsCur != NULL && !paramCurName.empty()) {
      pParamsCur->insert(make_pair(paramCurName, paramCurValue));
      paramCurName.clear();
      paramCurValue.clear();
    }
  }
}

void TandemParamsParser::characters(const XML_Char *s, int len) {
  if (pParamsCur != NULL && !paramCurName.empty()) {
    paramCurValue.append(s, len);
  }
}

bool TandemParamsParser::writePepXML(ostream& out) {
  // Write parameters
  out << nl() << "<!-- Input parameters -->";
  writePepXML(out, used);
  out << nl() << "<!-- Unused input parameters -->";
  writePepXML(out, unused);
  out << nl() << "<!-- Performance parameters -->";
  writePepXML(out, perf);
  return true;
}

void TandemParamsParser::writePepXML(ostream& out, xMap& params) {
  pair<string,string> pairValue;
  xMap::iterator it = params.begin();
  xMap::iterator end = params.end();
  while(it != end) {
    out << nl() << "<parameter name=\"" << XMLEscape(it->first) << "\""
	<< " value=\"" << XMLEscape(it->second) << "\"/>";
    it++;
  }
}

string TandemParamsParser::getSearchEnzyme() {
  const char* name = "protein, cleavage site";
  string site = getParam(name, used);
  if (site.length() == 0)
    site = getParam(name, unused);
  site = normalizeCleavageSite(site);

  string result;

  xMap::iterator it = enzymes.find(site);
  if (it != enzymes.end())
    result = it->second;

  return result;
}

string TandemParamsParser::getSampleEnzyme() {
  const char* name = "pipeline prophet, sample cleavage site";
  string site = getParam(name, used);
  if (site.length() == 0)
    site = getParam(name, unused);

  string result;

  if (site.length() == 0) {
    // If no explicit sample enzyme param, then use the search enzyme,
    // or default to trypsin, if using unconstrained cleavage.
    result = getSearchEnzyme();
    if (result.compare("unspecified") == 0)
      result = "trypsin";
    return result;
  }

  site = normalizeCleavageSite(site);
  xMap::iterator it = enzymes.find(site);
  if (it != enzymes.end())
    result = it->second;

  return result;
}

xMap TandemParamsParser::enzymes;

void TandemParamsParser::initEnzymes() {
  if (enzymes.size() == 0) {
    addEnzyme("argc", "[R]|{P}");
    addEnzyme("aspn", "[X]|[D]");
    addEnzyme("chymotrypsin", "[FLWY]|{P}");
    addEnzyme("trypsin/chymotrypsin", "[FLMWYKR]|{P}");
    addEnzyme("clostripain", "[R]|[X]");
    addEnzyme("cnbr", "[M]|[X]");
    addEnzyme("elastase", "[AGILV]|{P}");
    addEnzyme("formicacid", "[D]|{P}");
    addEnzyme("gluc", "[DE]|{P}");
    addEnzyme("gluc_bicarb", "[E]|{P}");
    addEnzyme("iodosobenzoate", "[W]|[X]");
    addEnzyme("lysc", "[K]|{P}");
    addEnzyme("lysc-p", "[K]|[X]");
    addEnzyme("lysn", "[X]|[K]");
    addEnzyme("lysn_promisc", "[X]|[KR]");
    addEnzyme("nonspecific", "[X]|[X]");
    addEnzyme("pepsina", "[FL]|[X]");
    addEnzyme("proline_endopeptidase", "[P]|[X]");
    addEnzyme("staph_protease", "[E]|[X]");
    addEnzyme("tca", "[KR]|{P},[FMWY]|{P},[X]|[D]");
    addEnzyme("trypsin", "[KR]|{P}");
    addEnzyme("trypsin/cnbr", "[KR]|{P},[M]|{P}");
    addEnzyme("trypsin_gluc", "[DEKR]|{P}");
    addEnzyme("trypsin_k", "[K]|{P}");
    addEnzyme("trypsin_r", "[R]|{P}");
  }
}

void TandemParamsParser::usageEnzymes() {
  cerr << "Use one of the following cleavage site specifications:\n";
  vector<string> enzymeNames;
  for (xMap::iterator it = enzymes.begin(); it != enzymes.end(); it++)
    enzymeNames.push_back(it->second + " - " + it->first);
  sort(enzymeNames.begin(), enzymeNames.end());
  for (vector<string>::iterator itName = enzymeNames.begin(); itName != enzymeNames.end(); itName++) {
    cerr << "\t" << *itName << "\n";
  }
}

void TandemParamsParser::addEnzyme(const char* name, const char* site) {
  // Look-up enzyme name by cleavage site spec.
  string nameString(name);
  string siteString(normalizeCleavageSite(site));
  enzymes.insert(make_pair(siteString, nameString));
}

string TandemParamsParser::normalizeCleavageSite(string site) {
  string result;

  size_t tStart = 0;
  size_t tComma = site.find(',');
  if (tComma != site.npos) {
    vector<string> sites;
    while (tStart != site.npos) {
      sites.push_back(normalizeCleavageSite(site.substr(tStart, tComma - tStart)));
      if ((tStart = tComma) != site.npos) {
	tStart++;
	tComma = site.find(',', tStart);
      }
    }
    sort(sites.begin(), sites.end());

    for (int i = 0; i < (int) sites.size(); i++) {
      if (i > 0)
	result.append(",");
      result.append(sites[i]);
    }
    return result;
  }

  size_t tPipe = site.find('|');
  if (tPipe == site.npos)
    return site;

  result.append(normalizeCleavagePart(site.substr(0, tPipe)));
  result.append("|");
  result.append(normalizeCleavagePart(site.substr(tPipe + 1)));
  return result;
}

string TandemParamsParser::normalizeCleavagePart(string part) {
  size_t last = part.length() - 1;
  if (part.length() < 4 ||
      (part.at(0) != '{' && part.at(0) != '[') ||
      (part.at(last) != '}' && part.at(last) != ']'))
    return part;

  string result = part;

  // n^2 sort for very short string of amino acid chars
  for (int i = 1; i < (int) last; i++) {
    for (int j = i + 1; j < (int) last; j++) {
      if (result[i] > result[j]) {
	char ch = result[i];
	result[i] = result[j];
	result[j] = ch;
      }
    }
  }

  return result;
}
