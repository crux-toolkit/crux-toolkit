// TandemResultsParser.cxx
//     Handles parsing of X!Tandem output file, and writing
//     to pepXML.

#include "Tandem2xml.h"
#include "TandemResultsParser.h"
#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#include "Parsers/Algorithm2XML/masscalc.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/util.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "gzstream.h"

#ifdef _MSC_VER
	#include <direct.h>
	#define getcwd(buff, n) _getcwd(buff, n)
#else
	#include <unistd.h>
#endif
#include <algorithm>

#ifdef COMET_EXACT
const double dProtonMass = 1.007825;
#else
const double dProtonMass = 1.007276;
#endif

// Allow slightly larger than 0.0001 tolerance to account for both rounding
// and the inability of double to represent values exactly.
// e.g. In double math
//		114.04300000000001
//    - 114.04290000000000
//    = 0.00010000000000331966
const double dModTolerance = 0.00010001;

using namespace mzParser;

bool compare_mods(const TandemResultsParser::ModData& i,  const TandemResultsParser::ModData& j) {
  if (i.aa < j.aa) {
    return true;
  }
  else if (i.aa > j.aa) {
    return false;
  }
  else {
    if (i.loc < j.loc) {
      return true;
    }
    else if (i.loc > j.loc) {
      return false;
    }
    else {
      if (i.modMass < j.modMass) {
	return true;
      }
      else if (i.modMass > j.modMass) {
	return false;
      }
    }
  }
}

#define makeupperc(c) ('a' <= c && c <= 'z' ? c += 'A' - 'a' : c)

/*
 * initializes the mass of amino acid residues (less water)
 *     chemical formulas found at http://haven.isb-sib.ch/tools/isotopident/htdocs/aa-list.html
 */
static void initAminoAcidMasses(double* pdMasses,
				masscalc::massType type = masscalc::monoisotopic, bool n15 = false)
{
  masscalc calc(n15);
  if (type == masscalc::monoisotopic)
    {
      pdMasses['a'] = pdMasses['A'] = calc.calcMass("C3H5ON");
      pdMasses['b'] = pdMasses['B'] = calc.calcMass("C4H6O2N2");// Same as N
      pdMasses['c'] = pdMasses['C'] = calc.calcMass("C3H5ONS");
      pdMasses['d'] = pdMasses['D'] = calc.calcMass("C4H5O3N");
      pdMasses['e'] = pdMasses['E'] = calc.calcMass("C5H7O3N");
      pdMasses['f'] = pdMasses['F'] = calc.calcMass("C9H9ON");
      pdMasses['g'] = pdMasses['G'] = calc.calcMass("C2H3ON");
      pdMasses['h'] = pdMasses['H'] = calc.calcMass("C6H7ON3");
      pdMasses['i'] = pdMasses['I'] = calc.calcMass("C6H11ON");
      pdMasses['j'] = pdMasses['J'] = 0.0;
      pdMasses['k'] = pdMasses['K'] = calc.calcMass("C6H12ON2");
      pdMasses['l'] = pdMasses['L'] = calc.calcMass("C6H11ON");
      pdMasses['m'] = pdMasses['M'] = calc.calcMass("C5H9ONS");
      pdMasses['n'] = pdMasses['N'] = calc.calcMass("C4H6O2N2");
      pdMasses['o'] = pdMasses['O'] = calc.calcMass("C4H6O2N2");// Same as N
      pdMasses['p'] = pdMasses['P'] = calc.calcMass("C5H7ON");
      pdMasses['q'] = pdMasses['Q'] = calc.calcMass("C5H8O2N2");
      pdMasses['r'] = pdMasses['R'] = calc.calcMass("C6H12ON4");
      pdMasses['s'] = pdMasses['S'] = calc.calcMass("C3H5O2N");
      pdMasses['t'] = pdMasses['T'] = calc.calcMass("C4H7O2N");
      pdMasses['u'] = pdMasses['U'] = 150.953640;// Why?
      pdMasses['v'] = pdMasses['V'] = calc.calcMass("C5H9ON");
      pdMasses['w'] = pdMasses['W'] = calc.calcMass("C11H10ON2");
      pdMasses['x'] = pdMasses['X'] = 111.060000;// Why?
      pdMasses['y'] = pdMasses['Y'] = calc.calcMass("C9H9O2N");
      pdMasses['z'] = pdMasses['Z'] = calc.calcMass("C5H8O2N2");// Same as Q
    }
  else    // masscalc doesn't work for average amino acid masses
    {
      /*
       * unfortunately, the average masses for amino acids do not
       * seem to be straight sums of the average masses for the atoms
       * they contain.
       *
       * instead of using the mass calculator, these numbers are taken
       * as constants from the web page referenced above.
       */
      pdMasses['a'] = pdMasses['A'] = 71.0788;
      pdMasses['b'] = pdMasses['B'] = 114.1038;	// Same as N
      pdMasses['c'] = pdMasses['C'] = 103.1388;
      pdMasses['d'] = pdMasses['D'] = 115.0886;
      pdMasses['e'] = pdMasses['E'] = 129.1155;
      pdMasses['f'] = pdMasses['F'] = 147.1766;
      pdMasses['g'] = pdMasses['G'] = 57.0519;
      pdMasses['h'] = pdMasses['H'] = 137.1411;
      pdMasses['i'] = pdMasses['I'] = 113.1594;
      pdMasses['j'] = pdMasses['J'] = 0.0;
      pdMasses['k'] = pdMasses['K'] = 128.1741;
      pdMasses['l'] = pdMasses['L'] = 113.1594;
      pdMasses['m'] = pdMasses['M'] = 131.1926;
      pdMasses['n'] = pdMasses['N'] = 114.1038;
      pdMasses['o'] = pdMasses['O'] = 114.1038;	// Same as N
      pdMasses['p'] = pdMasses['P'] = 97.1167;
      pdMasses['q'] = pdMasses['Q'] = 128.1307;
      pdMasses['r'] = pdMasses['R'] = 156.1875;
      pdMasses['s'] = pdMasses['S'] = 87.0782;
      pdMasses['t'] = pdMasses['T'] = 101.1051;
      pdMasses['u'] = pdMasses['U'] = 0.0;	// Why?
      pdMasses['v'] = pdMasses['V'] = 99.1326;
      pdMasses['w'] = pdMasses['W'] = 186.2132;
      pdMasses['x'] = pdMasses['X'] = 113.1594;	// Why?
      pdMasses['y'] = pdMasses['Y'] = 163.1760;
      pdMasses['z'] = pdMasses['Z'] = 128.1307;	// Same as Q
    }

#ifdef COMET_EXACT
  // COMET had different values for the non-amino-acid characters.
  pdMasses['b'] = pdMasses['B'] = 114.5349350;
  pdMasses['o'] = pdMasses['O'] = 114.0793126;
  pdMasses['u'] = pdMasses['U'] = 0.0;
  pdMasses['x'] = pdMasses['X'] = 113.0840636;
  pdMasses['z'] = pdMasses['Z'] = 128.5505850;

  // Deal with rounding error.
  for (int i = 0; i < 128; i++)
    {
      char buffer[20];
      double d = pdMasses[i];
      sprintf(buffer, "%0.7lf", d);
      pdMasses[i] = atof(buffer);
    }
#endif
}

static const char EL_GROUP[] = "group";
static const char EL_PROTEIN[] = "protein";
static const char EL_PEPTIDE[] = "peptide";
static const char EL_DOMAIN[] = "domain";
static const char EL_NOTE[] = "note";
static const char EL_AA[] = "aa";

static const char STATE_MODEL[] = "model";
static const char STATE_SUPPORT[] = "support";
static const char STATE_PROT_DESC[] = "protein_description";
static const char STATE_SCAN_DESC[] = "scan_description";


TandemResultsParser::TandemResultsParser()
{
  instrument = NULL;
  searchEnzyme = NULL;
  sampleEnzyme = NULL;

  masscalc calc;
  dNTerm = dNTermDefault = calc.calcMass("H");
  dCTerm = dCTermDefault = calc.calcMass("OH");
  dNTermProt = 0.0;
  dCTermProt = 0.0;

  memset(dMassAA, 0, sizeof(dMassAA));

  m_bCometScoring = false;
  m_pout = &cout;

  scanCount = 0;
  groupDepth = 0;
  state = NULL;

  pf = NULL;
  mzXMLindices = NULL;
  lastIndex = 0;
  useN15 = false;
}

TandemResultsParser::~TandemResultsParser()
{
  if (instrument != NULL)
    free(instrument);
  if (searchEnzyme != sampleEnzyme)
    delete searchEnzyme;
  delete sampleEnzyme;

  if (pf != NULL) 
    rampCloseFile(pf);
  if (mzXMLindices != NULL) 
    free(mzXMLindices); 

}

void TandemResultsParser::startElement(const XML_Char *el, const XML_Char **attr) {
  if (isElement(EL_GROUP, el))
    startGroup(attr);
  else if (isElement(EL_PROTEIN, el)) {
    if (state == STATE_MODEL)
      startProtein(attr);
  }
  else if (isElement(EL_PEPTIDE, el)) {
    if (state == EL_PROTEIN)
      state = EL_PEPTIDE;
  }
  else if (isElement(EL_DOMAIN, el)) {
    if (state == EL_PEPTIDE)
      startDomain(attr);
  }
  else if (isElement(EL_AA, el)) {
    if (state == EL_DOMAIN)
      startAA(attr);
  }
  else if (isElement(EL_NOTE, el)) {
    if (state == STATE_SUPPORT || state == EL_PROTEIN)
      startNote(attr);
  }
}

void TandemResultsParser::endElement(const XML_Char *el) {
  if (isElement(EL_GROUP, el))
    endGroup();
  else if (isElement(EL_PROTEIN, el)) {
    if (state == EL_PROTEIN)
      endProtein();
  }
  else if (isElement(EL_PEPTIDE, el)) {
    if (state == EL_PEPTIDE)
      state = EL_PROTEIN;
  }
  else if (isElement(EL_DOMAIN, el)) {
    if (state == EL_DOMAIN)
      endDomain();
  }
  else if (isElement(EL_NOTE, el))
    endNote();
}

void TandemResultsParser::startGroup(const XML_Char **attr)
{
  const char* groupType = getAttrValue("type", attr);

  if (groupDepth == 0 && strcmp(groupType, "model") == 0)
    {
      state = STATE_MODEL;

      // Reset all values
      scanCur.clear();

      // Look for values on this tag.  Other values on child tags.
      for (int i = 0; attr[i]; i += 2)
	{
	  const char* name = attr[i];
	  const char* value = attr[i + 1];

	  if (isAttr("id", name))
	    {
	      scanCur.id = ++scanCount;

	      int index = atoi(value);
	      // X!Tandem has no built in facility for noting
	      // start and end scan, so default to ID.
	      scanCur.scanStart = index+indexOffset;
	      scanCur.scanEnd = index+indexOffset;
	    }
	  else if (isAttr("mh", name))
	    {
	      scanCur.mh = atof(value);
	    }
	  else if (isAttr("rt", name))
	    {
	      scanCur.retentionTime = value;
	    }
	  else if (isAttr("z", name))
	    {
	      scanCur.charge = (short) atoi(value);
	    }
	  else if (isAttr("label", name))
	    {
	      scanCur.prot.assign(value);
	    }
	}
    }
  else if (groupDepth == 1 && state == STATE_MODEL &&
	   strcmp(groupType, "support") == 0)
    {
      state = STATE_SUPPORT;
    }

  groupDepth++;
}

void TandemResultsParser::endGroup()
{
  groupDepth--;

  if (groupDepth == 0 && state == STATE_MODEL)
    {
#if 0
      if (scanCur.vScores.size() > 1)
	{
	  /* Sort peptide scores matching this peptide to the top. */
	  vector<ScoreData>::iterator begin = scanCur.vScores.begin();
	  vector<ScoreData>::iterator it = scanCur.vScores.begin() + 1;
	  while (it != scanCur.vScores.end())
	    {
	      if (it->seq.compare(begin->seq) == 0 &&
		  lessThanScoreMods(*it, *begin))
		{
		  
		  ScoreData score = *it;
		  *it = *begin;
		  *begin = score;
		}
	      it++;
	    }
	}
#endif
      writePepXML(*m_pout, scanCur);

      state = NULL;
    }
  else if (groupDepth == 1 && state == STATE_SUPPORT)
    {
      state = STATE_MODEL;
    }
}

void TandemResultsParser::startProtein(const XML_Char **attr)
{
  state = EL_PROTEIN;

  scoreCur.clear();

  scoreCur.prot.assign(getAttrValue("label", attr));
}

void TandemResultsParser::endProtein()
{
  state = STATE_MODEL;
}

void TandemResultsParser::startDomain(const XML_Char **attr)
{
  state = EL_DOMAIN;

  string proteinName(scoreCur.prot);
  scoreCur.clear();
  scoreCur.prot.assign(proteinName);

  for (int i = 0; attr[i]; i += 2)
    {
      const char* name = attr[i];
      const char* value = attr[i + 1];

      if (isAttr("pre", name))
	{
	  scoreCur.seqPre.assign(value);
	  if (scoreCur.seqPre.compare("[") == 0)
	    scoreCur.seqPre.assign("-");
	}
      else if (isAttr("post", name))
	{
	  scoreCur.seqPost.assign(value);
	  if (scoreCur.seqPost.compare("]") == 0)
	    scoreCur.seqPost.assign("-");
	}
      else if (isAttr("seq", name))
	{
	  scoreCur.seq.assign(value);
	}
      else if (isAttr("start", name))
	{
	  scoreCur.seqStart = atoi(value);
	}
      else if (isAttr("end", name))
	{
	  scoreCur.seqEnd = atoi(value);
	}
      else if (isAttr("expect", name))
	{
	  scoreCur.expect = atof(value);
	}
      else if (isAttr("hyperscore", name))
	{
	  scoreCur.hyperScore.assign(value);
	}
      else if (isAttr("nextscore", name))
	{
	  scoreCur.nextScore.assign(value);
	}
      else if (isAttr("mh", name))
	{
	  scoreCur.seqMh = atof(value);
	}
      else if (isAttr("y_score", name))
	{
	  scoreCur.yScore = atof(value);
	}
      else if (isAttr("y_ions", name))
	{
	  scoreCur.yIons = atoi(value);
	}
      else if (isAttr("b_score", name))
	{
	  scoreCur.bScore = atof(value);
	}
      else if (isAttr("b_ions", name))
	{
	  scoreCur.bIons = atoi(value);
	}
      else if (isAttr("z_score", name))
	{
	  scoreCur.zScore = atof(value);
	}
      else if (isAttr("z_ions", name))
	{
	  scoreCur.zIons = atoi(value);
	}
      else if (isAttr("c_score", name))
	{
	  scoreCur.cScore = atof(value);
	}
      else if (isAttr("c_ions", name))
	{
	  scoreCur.cIons = atoi(value);
	}
      else if (isAttr("x_score", name))
	{
	  scoreCur.xScore = atof(value);
	}
      else if (isAttr("x_ions", name))
	{
	  scoreCur.xIons = atoi(value);
	}
      else if (isAttr("a_score", name))
	{
	  scoreCur.aScore = atof(value);
	}
      else if (isAttr("a_ions", name))
	{
	  scoreCur.aIons = atoi(value);
	}
      else if (isAttr("missed_cleavages", name))
	{
	  scoreCur.missedCleavages = atoi(value);
	}
    }
}

void TandemResultsParser::endDomain()
{
  sort(scoreCur.vMods.begin(), scoreCur.vMods.end(),
       lessThanModLoc);
  scanCur.vScores.push_back(scoreCur);

  state = EL_PEPTIDE;
}

void TandemResultsParser::startAA(const XML_Char **attr)
{
  ModData modCur;

  for (int i = 0; attr[i]; i += 2)
    {
      const char* name = attr[i];
      const char* value = attr[i + 1];

      if (isAttr("type", name))
	{
	  modCur.aa = *value;
	}
      else if (isAttr("at", name))
	{
	  modCur.loc = atoi(value);
	}
      else if (isAttr("modified", name))
	{
	  modCur.modMass = atof(value);
	}
      else if (isAttr("pm", name))
	{
	  modCur.pm = *value;
	}
    }

  scoreCur.vMods.push_back(modCur);
}

void TandemResultsParser::startNote(const XML_Char **attr)
{
  const char *noteLabel = getAttrValue("label", attr);
  if (state == EL_PROTEIN && strcmp("description", noteLabel) == 0)
    {
      state = STATE_PROT_DESC;

      // At this point prot should contain the group label
      // which may be abreviated.  The note should contain
      // the full text.
      scoreCur.prot.clear();
    }
  else if (state == STATE_SUPPORT && strcmp("Description", noteLabel) == 0)
    {
      state = STATE_SCAN_DESC;

      // The spectrum buffer is probably alread clear, but
      // just to be sure.
      scanCur.description.clear();
    }
}

void TandemResultsParser::endNote()
{
  if (state == STATE_PROT_DESC)
    {
      // Remove trailing white space.  Sometimes X! Tandem appends a
      // newline to the end, which can be problematic later on in the TPP.
      string::iterator it = scoreCur.prot.end();
      string::iterator begin = scoreCur.prot.begin();
      while (it != begin && isspace(*(it - 1)))
	it--;
      if (it != scoreCur.prot.end())
	scoreCur.prot.erase(it, scoreCur.prot.end());

      // Set state to parent tag.
      state = EL_PROTEIN;
    }

  if (state != STATE_SCAN_DESC)
    return;

  // Erase leading whitespace.
  string::iterator it = scanCur.description.begin();
  string::iterator end = scanCur.description.end();
  while (it != end && isspace(*it))
    it++;
  scanCur.description.erase(scanCur.description.begin(), it);

  size_t dot = scanCur.description.find('.');
  size_t dta = scanCur.description.find(".dta");
  // handle all RAMP supported types
  //const char **rampListSupportedFileTypes();
  size_t xml = string::npos;
  size_t xmlAfter = string::npos;
  int longest = 0;
  const char **RAMPFileTypes = rampListSupportedFileTypes();
  for (int i=0;RAMPFileTypes[i];i++) {
    size_t ixml = scanCur.description.find(RAMPFileTypes[i]);
    if (ixml != string::npos) {
      int len = (int)strlen(RAMPFileTypes[i]);
      if (len > longest) {
	xmlAfter = ixml + len + 1;
	longest = len;
	xml = ixml;
      }
    }
  }

  // See if it is special .dta format.
  if (dot != string::npos && dta != string::npos && dot < dta)
    {
      int scanStart, scanEnd, charge;

      dot = dta;
      for (int i = 0; i < 3 && dot != string::npos; i++)
	dot = scanCur.description.rfind('.', dot - 1);

      if (dot != string::npos)
	{
	  int count = sscanf(scanCur.description.data() + dot + 1, "%d.%d.%d",
			     &scanStart, &scanEnd, &charge);
	  if (count == 3)
	    {
	      scanCur.scanStart = scanStart;
	      scanCur.scanEnd = scanEnd;
	    }
	}

      scanCur.description.erase(dta);
    }
  // Try special mzXML / mzData format.
  else if (xml != string::npos)
    {
      int scan, charge;
      int count = sscanf(scanCur.description.data() + xmlAfter,
			 "scan %d (charge %d)", &scan, &charge);
      scanCur.description.erase(xml + 1);
      char szBuff[128];
      sprintf(szBuff, "%05d.%05d.%d", scan, scan, charge);
      scanCur.description.append(szBuff);
    }

  state = STATE_SUPPORT;
}

void TandemResultsParser::characters(const XML_Char *s, int len)
{
  if (len <= 0)
    return;

  if (state == STATE_PROT_DESC)
    {
      scoreCur.prot.append(s, len);
    }
  else if (state == STATE_SCAN_DESC)
    {
      scanCur.description.append(s, len);
    }
}

bool TandemResultsParser::writePepXML()
{
  paramsHandler.setFileName(m_strFileName.data());
  if (!paramsHandler.parse())
    return false;

  // Initialize masses
  
 

  if (useN15) {
    h2o = 18.0105647;
    nh3 = 18.02358311;
  }
  
  masscalc::massType type = masscalc::monoisotopic;
  const char* fragmentType = "monoisotopic";
  string strType = paramsHandler.getUsedParam("spectrum, fragment mass type");
  if (strType == "average") {
    type = masscalc::average;
    fragmentType = "average";
  }
  initAminoAcidMasses(dMassAA, type, useN15);

  ProteolyticEnzymeFactory factory;

  // Initialize enzymes
  searchEnzymeName = paramsHandler.getSearchEnzyme();
  if (searchEnzymeName.length() != 0)
    searchEnzyme = factory.getProteolyticEnzyme(searchEnzymeName.data());
  /* TODO: Add else with warning, and construct custom PE where possible. */

  if (searchEnzyme == NULL)
    {
      cerr << "ERROR: Unsupported search enzyme.\n";	// TODO: Do better here.
      paramsHandler.usageEnzymes();
      return false;
    }

  sampleEnzymeName = paramsHandler.getSampleEnzyme();
  if (searchEnzymeName.compare(sampleEnzymeName) == 0)
    sampleEnzyme = searchEnzyme;
  else
    sampleEnzyme = factory.getProteolyticEnzyme(sampleEnzymeName.data());

  if (sampleEnzyme == NULL)
    {
      cerr << "ERROR: Unsupported sample enzyme.\n";
      return false;
    }

  // Try to get machine information.
  //RAMPFILE* pf = NULL;
  string strBaseName;
  const char *rawDataXMLext = NULL;

  string strMzXMLFile = paramsHandler.getUsedParam("spectrum, path");
  if (strMzXMLFile.length() > 0)
    {
      rawDataXMLext = rampValidFileType(strMzXMLFile.data());
      pf = rampOpenFile(strMzXMLFile.data());
    }

  // Try making a good guess in the local directory.
  //Fred Hutch specific code.  If the mzXML file isn't found in the  same directory as the xtan.xml file
  //try looking in .., ../.., and ../../xml

  if (pf == NULL)
    {
      strMzXMLFile = m_strFileName;
      unsigned int pepXMLExtLen = strlen(get_pepxml_dot_ext());
      if (strMzXMLFile.length() >= pepXMLExtLen && strMzXMLFile.compare(strMzXMLFile.length() - pepXMLExtLen, pepXMLExtLen, get_pepxml_dot_ext()) == 0)
	strMzXMLFile.erase(strMzXMLFile.length() - pepXMLExtLen);
      if (strMzXMLFile.length() >= 4 && strMzXMLFile.compare(strMzXMLFile.length() - 4, 4, ".xml") == 0)
	strMzXMLFile.erase(strMzXMLFile.length() - 4);
      if (strMzXMLFile.length() >= 5 && strMzXMLFile.compare(strMzXMLFile.length() - 5, 5, ".xtan") == 0)
	strMzXMLFile.erase(strMzXMLFile.length() - 5);
      int len;
      char *fname=new char[len=strMzXMLFile.length()+21];
      // DCT fix: buffer needs to have additional 21 chars
      //          longest prefix tried is /../../xml/         11 chars
      //          longest suffix tried is .mzXML.gz           9 chars 

      rampConstructInputFileName(fname,len,strMzXMLFile.data());
      pf = rampOpenFile(fname);

      if(pf == NULL) {
	string oneUp = "../"+strMzXMLFile;
	rampConstructInputFileName(fname,len,oneUp.data());
	pf = rampOpenFile(fname);                
      }
      if(pf == NULL) {
	string twoUp = "../../"+strMzXMLFile;
	rampConstructInputFileName(fname,len,twoUp.data());
	pf = rampOpenFile(fname);                
      }
      if(pf == NULL) {
	string xmlUp = "../../xml/"+strMzXMLFile;
	rampConstructInputFileName(fname,len,xmlUp.data());
	pf = rampOpenFile(fname);                
      }

      if(pf != NULL) {
	strMzXMLFile = fname;
	rawDataXMLext = rampValidFileType(strMzXMLFile.data());
      }

      delete [] fname;
    }

  if (pf != NULL)
    {
      instrument = getInstrumentStruct(pf);
      mzXMLindices = readIndex(pf,getIndexOffset(pf),&lastIndex);

      if (mzXMLindices == NULL)
	{
	  cerr << "ERROR: Failed reading scan indices from " << strMzXMLFile << " file.\n";
	  return false;
	}
      //rampCloseFile(pf);
    }
  else
    {
      cerr << "WARNING: Failed to open mzML file.\n"
	   << "         Output will not contain retention times.\n";
    }
  if (rawDataXMLext == NULL)
    {
      strBaseName = m_strFileName;
      rawDataXMLext = ".?";
    }
  else
    {
      strBaseName = strMzXMLFile.substr(0, strMzXMLFile.length() - strlen(rawDataXMLext));
    }
 
  // Unless an output file was specified
  ogzstream ofOut; // write as gzip if filename so indicates
  if (!m_strOutputFile.empty())
    {
      ofOut.open(m_strOutputFile.data()); // write as gzip if filename so indicates
      if (!ofOut.is_open())
	{
	  cerr << "Failed to open " << m_strOutputFile << " for write.\n";
	  return false;
	}

      m_pout = &ofOut;
    }

  bool success = writePepXML(*m_pout, strBaseName.data(), rawDataXMLext, fragmentType);

  if (ofOut.is_open())
    ofOut.close();

  return success;
}

bool TandemResultsParser::writePepXML(ostream& out, const char *baseName, const char *rawDataXMLext, const char *fragmentType)
{

  
  string scoringEngine = "X! Tandem";
  string scoringAlgorithm = paramsHandler.getUsedParam("scoring, algorithm");
  if (scoringAlgorithm.length() != 0)
    {
      scoringEngine.append(" (");
      scoringEngine.append(scoringAlgorithm);
      scoringEngine.append(")");
    }

  m_bCometScoring = (scoringAlgorithm.compare("comet") == 0);

  string dateStart = paramsHandler.getPerfParam("process, start time");
  if (string::npos == dateStart.find('T')) { // not xsd:dateTime
    // convert, for example, "2011:07:14:16:51:35" to "2011-07-14T16:51:35"
    istringstream date(dateStart);
    int y,m,d,hh,mm,ss;
    char c; // soaks up the :
    date >> y >> c >> m >> c >> d >> c >> hh >> c >> mm >> c >> ss;
    ostringstream datetime;
    datetime << setfill('0') << setw(4) << y << "-" << setw(2) << m <<  "-" << setw(2) << d << "T" 
	     << setw(2) << hh << ":" << setw(2) << mm << ":" << setw(2) << ss;
    dateStart = datetime.str();
  }

  string pathSummary;
  if (!m_strOutputFile.empty())
    {
      // Deal with absolute paths
      if (!isAbsolutePath(m_strOutputFile))
	{
	  char buffer[4096];
	  char *gret=getcwd(buffer, _countof(buffer) - 1);
	  pathSummary.append(buffer);
	  char c = pathSummary.at(pathSummary.length() - 1);
	  if (c != '\\' && c != '/')
	    {
	      if (pathSummary.find('\\') != string::npos)
		pathSummary.append("\\");
	      else
		pathSummary.append("/");
	    }
	}
      pathSummary.append(m_strOutputFile);
    }

  char buffer[256];

  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << nl();
  //out << "<!DOCTYPE msms_pipeline_analysis SYSTEM \"msms_analysis3.dtd\">" << nl();
  out << "<?xml-stylesheet type=\"text/xsl\" href=\"pepXML_std.xsl\"?>" << nl();
  out << "<msms_pipeline_analysis date=\"" << dateStart << "\""
      << " summary_xml=\"" << pathSummary <<"\""
      << " xmlns=\"http://regis-web.systemsbiology.net/pepXML\""
      << " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
      << " xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/pepXML/" << PEPXML_SCHEMA << "\""
      << ">" << nlIn();

  string strBaseName = m_strOutputFile;
  const char *ext = hasValidPepXMLFilenameExt(strBaseName.c_str());
  if (ext)
    strBaseName.erase(strBaseName.length() - strlen(ext));
  // Just in case ".xml" becomes not valid pepXML ext.
  else if (strBaseName.length() >= 4 && strBaseName.compare(strBaseName.length() - 4, 4, ".xml") == 0)
    strBaseName.erase(strBaseName.length() - 4);

  out << "<msms_run_summary base_name=\"" << baseName << "\"";

  if (instrument != NULL)
    {
      out			<< " msManufacturer=\"" << instrument->manufacturer << "\""
				<< " msModel=\"" << instrument->model << "\"" 
				<< " msIonization=\"" << instrument->ionisation << "\""
				<< " msMassAnalyzer=\"" << instrument->analyzer << "\""
				<< " msDetector=\"" << instrument->detector << "\"";
    }

  out 				<< " raw_data_type=\"raw\""
				<< " raw_data=\"" << rawDataXMLext << "\">" << nlIn();

  sampleEnzyme->writePepXMLTags(out);

  out	<< "<search_summary base_name=\"" << baseName << "\""
	<< " search_engine=\"" << scoringEngine << "\""
	<< " precursor_mass_type=\"monoisotopic\""
	<< " fragment_mass_type=\"" << fragmentType << "\""
    //					<< " out_data_type=\"xml\""
    //					<< " out_data=\".xml\""
	<< " search_id=\"1\">" << nlIn();

  // Write database files
  string paramCurName = "list path, sequence source #";
  int len = (int) paramCurName.length();
  int i;
  for (i = 1; ; i++)
    {
      paramCurName.erase(len);
      sprintf(buffer, "%d", i);
      paramCurName.append(buffer);

      string source = paramsHandler.getPerfParam(paramCurName.data());
      if (source.empty())
	break;

      out << "<search_database local_path=\"" << source << "\""
	  << " type=\"AA\"/>" << nl();
    }

  // Search enzyme
  string missedCleavageMaxParam =
    paramsHandler.getUsedParam("scoring, maximum missed cleavage sites");
  int missedCleavageMax = atoi(missedCleavageMaxParam.data());

  int minCleavageTermini = 2;

  string semiCleavage =
    paramsHandler.getUsedParam("protein, cleavage semi");

  if (semiCleavage.compare("yes") == 0)
    minCleavageTermini = 1;
  if (searchEnzymeName.compare("unspecified") == 0)
    minCleavageTermini = 0;

  semiCleavage =
    paramsHandler.getUsedParam("refine");

  if (semiCleavage.compare("yes") == 0) {
    int term = 2;
    semiCleavage =
      paramsHandler.getUsedParam("refine, cleavage semi");

    if (semiCleavage.compare("yes") == 0)
      term = 1;

    semiCleavage =
      paramsHandler.getUsedParam("refine, unanticipated cleavage");
	  
    if (searchEnzymeName.compare("yes") == 0)
      term = 0;

    if (term < minCleavageTermini) 
      minCleavageTermini = term;

  }
		
  out		<< "<enzymatic_search_constraint"
		<<		" enzyme=\"" << searchEnzymeName << "\""
		<<		" max_num_internal_cleavages=\"" << missedCleavageMax << "\""
		<<		" min_number_termini=\"" << minCleavageTermini << "\""
		<<		" />" << nl();

  // Modifications
  // Peptide static terminal
  string modTerm = paramsHandler.getUsedParam("protein, cleavage N-terminal mass change");
  if (!modTerm.empty())
    dNTerm = atof(modTerm.data());

  modTerm = paramsHandler.getUsedParam("protein, cleavage C-terminal mass change");
  if (!modTerm.empty())
    dCTerm = atof(modTerm.data());

  // Protein static terminal
  modTerm = paramsHandler.getUsedParam("protein, N-terminal residue modification mass");
  if (!modTerm.empty())
    dNTermProt = atof(modTerm.data());

  modTerm = paramsHandler.getUsedParam("protein, C-terminal residue modification mass");
  if (!modTerm.empty())
    dCTermProt = atof(modTerm.data());

  // Static modifications
  readMods(paramsHandler.getUsedParam("residue, modification mass"), false, 0);

  // Special X! Tandem n-terminal AA variable modifications
  ModSpecData modSpec;
  modSpec.symbol = '^';
  modSpec.comment = "X! Tandem n-terminal AA variable modification";

  modSpec.aa = 'E';
  modSpec.mass = 0-h2o;
  addModSpec(modSpec, true);
  modSpec.aa = 'Q';
  modSpec.mass = 0-nh3;
  addModSpec(modSpec, true);
  if ((long) modsStatic['C'].mass == 57)
    {
      modSpec.aa = 'C';
      modSpec.mass = 0-nh3;
      addModSpec(modSpec, true);
    }
  // deal with quick acetyl
  modSpec.aa = '[';
  modSpec.mass = 42.0106;
  addModSpec(modSpec, true);

  // Variable modifications
  readMods(paramsHandler.getUsedParam("residue, potential modification mass"), true, 0);

  char residueParam[256];
  int pass = 1;
  for (;;)
    {
      sprintf(residueParam, "residue, modification mass %d", pass);
      if (!readMods(paramsHandler.getUsedParam(residueParam), true, pass++))
	break;
    }

  // Refine variable modifications
  readMods(paramsHandler.getUsedParam("refine, potential modification mass"), true, 0);

  char refineParam[256];
  pass = 1;
  for (;;)
    {
      sprintf(refineParam, "refine, potential modification mass %d", pass++);
      if (!readMods(paramsHandler.getUsedParam(refineParam), true, 0))
	break;
    }
  readMods(paramsHandler.getUsedParam("refine, potential N-terminus modifications"), true, 0);
  readMods(paramsHandler.getUsedParam("refine, potential C-terminus modifications"), true, 0);

  // Motif variable modifications
  readMotifs(paramsHandler.getUsedParam("residue, potential modification motif"));

  // Refine motif variable modifications
  readMotifs(paramsHandler.getUsedParam("refine, potential modification motif"));
	
  pass = 1;
  for (;;)
    {
      sprintf(refineParam, "refine, potential modification motif %d", pass++);
      if (!readMotifs(paramsHandler.getUsedParam(refineParam)))
	break;
    }

  writeModSpecs(out);

  // X!Tandem parameters
  if (!paramsHandler.writePepXML(out))
    return false;

  out << nlOut();
  out << "</search_summary>" << nl();

  if (!parse())
    return false;

  out << nlOut();
  out << "</msms_run_summary>";
  out << nlOut();
  out << "</msms_pipeline_analysis>"
      << nl();

  return true;
}

bool TandemResultsParser::readMods(string mod, bool variable, int pass)
{
  if (mod.length() < 3)
    return false;

  ModSpecData modSpec;

  size_t tStart = 0;
  while (tStart != mod.npos)
    {
      string next = mod.substr(tStart, mod.find(',', tStart) - tStart);
      modSpec.mass = atof(next.data());

      if (modSpec.mass == 0.0)
	break;	// X!Tandem stops in this case.



      size_t tAt = next.find('@');
      if(tAt == next.npos)
	break;
      tAt++;

      modSpec.aa = makeupperc(next.at(tAt));


      if (variable && pass > 0) {
	//This is a second or higher pass of static mods, subtract eh first pass static mod 
	//and create a variable new mod 
	if (modsStatic[modSpec.aa].mass != 0.0) {
	  modSpec.mass -= modsStatic[modSpec.aa].mass;

	  if (modSpec.mass >= 0.0000001) {
	    addModSpec(modSpec, variable);
	  } 
	}
      }
      else {
	addModSpec(modSpec, variable);
      }

      tStart = mod.find(',', tStart);
      if (tStart != mod.npos)
	tStart++;	// advance past comma
    }

  return true;
}

bool TandemResultsParser::readMotifs(string motif)
{
  if (motif.length() < 3)
    return false;

  ModSpecData modSpec;

  size_t tStart = 0;
  while (tStart != motif.npos)
    {
      string next = motif.substr(tStart, motif.find(',', tStart) - tStart);
      modSpec.mass = atof(next.data());
      if (modSpec.mass == 0.0)
	break;	// X!Tandem stops in this case.

      // TODO: Deal with colon

      size_t tAt = next.find('@');
      if(tAt == next.npos)
	break;
      tAt++;

      size_t tAA = tAt;
      size_t tBang = next.find('!', tAt);
      if (tBang != next.npos)
        {
	  // Walk successive alpha characters looking for an open bracket
	  for (tAA = tBang - 1; tAA >= tAt; tAA--)
            {
	      if (next.at(tAA) == '[')
		break;
	      if (!isalpha(next.at(tAA)))
		tAA = tAt;
            }

	  if (tAA < tAt)
	    tAA = tBang - 1;
        }

      if (next.at(tAA) == '[')
        {
	  tAA++;
	  while (tAA < next.npos && next.at(tAA) != ']')
            {
	      char aa = next.at(tAA++);
	      if (isalpha(aa))
                {
		  modSpec.aa = makeupperc(aa);
		  modSpec.comment.assign("motif ").append(next);
		  addModSpec(modSpec, true);
                }
	      else if (aa != '!')
                {
		  cerr << "WARNING: Failed reading motif '" << next << "'.\n";
		  break;
                }
            }
        }
      else if (isalpha(next.at(tAA)))
        {
	  modSpec.aa = makeupperc(next.at(tAA));
	  modSpec.comment.assign("motif ").append(next);
	  addModSpec(modSpec, true);
        }
      else
        {
	  cerr << "WARNING: Failed to interpret motif '" << next << "'.\n";
        }

      tStart = motif.find(',', tStart);
      if (tStart != motif.npos)
	tStart++;	// advance past comma
    }
  return true;
}

void TandemResultsParser::writeModSpecs(ostream& out)
{
  int c;
  // Write static amino acid modification specifications
  for (c = 'A'; c < 'Z'; c++)
    writeAAModSpec(out, modsStatic[c], false);

  // Write variable amino acid modification specifications
  for (c = 'A'; c < 'Z'; c++)
    {
      vector<ModSpecData>& mods = vModsVariable[c];
      for (int i = 0; i < (int) mods.size(); i++)
	writeAAModSpec(out, mods[i], true);
    }

  // Write static terminal modifications
  ModSpecData mod;
  double dModProt = dNTermProt;
  if (modsStatic['['].mass != 0.0)
    {
      dModProt += modsStatic['['].mass;
      writeTermModSpec(out, modsStatic['['], false, true);
    }
  else if (fabs(dNTerm - dNTermDefault) > dModTolerance)
    {
      mod.aa = '[';
      mod.mass = 0.0;
      writeTermModSpec(out, mod, false, true);
    }
  if (dNTermProt != 0.0)
    {
      mod.aa = '[';
      mod.mass = dModProt;
      writeTermModSpec(out, mod, false, false);
    }

  dModProt = dCTermProt;
  if (modsStatic[']'].mass != 0.0)
    {
      dModProt = modsStatic[']'].mass;
      writeTermModSpec(out, modsStatic[']'], false, true);
    }
  else if (fabs(dCTerm - dCTermDefault) > dModTolerance)
    {
      mod.aa = ']';
      mod.mass = 0.0;
      writeTermModSpec(out, mod, false, true);
    }
  if (dCTermProt != 0.0)
    {
      mod.aa = ']';
      mod.mass = dModProt;
      writeTermModSpec(out, mod, false, false);
    }

  // Write peptide variable terminal modifications
  int i;
  vector<ModSpecData>& modsN = vModsVariable['['];
  for (i = 0; i < (int) modsN.size(); i++)
    writeTermModSpec(out, modsN[i], true, true);
  vector<ModSpecData>& modsC = vModsVariable[']'];
  for (i = 0; i < (int) modsC.size(); i++)
    writeTermModSpec(out, modsC[i], true, true);
}

void TandemResultsParser::writeAAModSpec(ostream& out, const ModSpecData& mod, bool variable)
{
  if (mod.mass == 0.0)
    return;

  double massTotal = dMassAA[mod.aa] + mod.mass;
  if (variable)
    massTotal += modsStatic[mod.aa].mass;

  char buffer[128];
  out << "<aminoacid_modification"
      << " aminoacid=\"" << mod.aa << "\"";
  sprintf(buffer, "%0.4f", mod.mass);
  out		<< " massdiff=\"" << buffer << "\"";
  sprintf(buffer, "%0.4f", massTotal);
  out		<< " mass=\"" << buffer << "\"";
  out     << " variable=\"" << (variable ? "Y" : "N") << "\"";
  if (mod.symbol != '\0')
    out	<< " symbol=\"" << mod.symbol << "\"";
  out		<< " />";
  if (mod.comment.length() > 0)
    out << "<!--" << mod.comment << "-->";
  out << nl();
}

void TandemResultsParser::writeTermModSpec(ostream& out, const ModSpecData& mod, bool variable, bool peptide)
{
  double massTotal = mod.mass + (mod.aa == '[' ? dNTerm : dCTerm);
  double massMod = massTotal - (mod.aa == '[' ? dNTermDefault : dCTermDefault);

  if (variable)
    {
      // If it is variable, make sure all the static-mod mass is added in.
      // And, make the mod mass the value the user specified.
      massTotal += modsStatic[mod.aa].mass;
      massMod = mod.mass;
    }

  char buffer[128];
  out << "<terminal_modification"
      << " terminus=\"" << (mod.aa == '[' ? 'n' : 'c') << "\"";
  sprintf(buffer, "%0.4f", massMod);
  out		<< " massdiff=\"" << buffer << "\"";
  sprintf(buffer, "%0.4f", massTotal);
  out		<< " mass=\"" << buffer << "\"";
  out     << " protein_terminus=\"" << (peptide ? "N" : "Y") << "\"";
  out     << " variable=\"" << (variable ? "Y" : "N") << "\"";
  if (mod.symbol != '\0')
    out	    << " symbol=\"" << mod.symbol << "\"";
  out		<< " />";
  if (mod.comment.length() > 0)
    out << "<!--" << mod.comment << "-->";
  out << nl();
}

void TandemResultsParser::addModSpec(ModSpecData mod, bool variable)
{
  // YUCK!  X!Tandem writes the modification mass directly to the output
  //        stream, which always yields 6-digit output.  This makes the
  //        modification tolerance too small for numbers with greater than
  //        2 digits left of the decimal point.  The ugly solution here
  //        is to use a string stream to generate the string as X!Tandem
  //        would write it to the output stream, and then convert that
  //        to a double value that should match what we see in the <aa>
  //        output.
  stringstream outputPredictorStream;
  outputPredictorStream << mod.mass;
  mod.massPredictedOutput = atof(outputPredictorStream.str().data());

  if (!variable)
    {
      // As in X! Tandem, repeated static mods overwrite the original.
      modsStatic[mod.aa] = mod;
    }
  // Add only variable modifications not already seen.  Refinement
  // will often repeat modifications.
  else if (matchMod(mod.aa, mod.mass, true) == 0.0)
    {
      vModsVariable[mod.aa].push_back(mod);
    }
}

double TandemResultsParser::matchMod(char c, double massMod, bool variable)
{
  // Below, we check both the originally declared modification mass,
  // as well as the predicted output, if X! Tandem continues to write
  // these masses directly to an out stream.  Not great, but the best we
  // can do.

  if (variable)
    {
      vector<ModSpecData>& mods = vModsVariable[c];
      for (size_t i = 0; i < mods.size(); i++)
        {
	  if (fabs(massMod - (mods[i].mass+modsStatic[c].mass)) <= dModTolerance) 
	    return mods[i].mass+modsStatic[c].mass;
	  
	  if (fabs(massMod - mods[i].mass) <= dModTolerance ||
	      fabs(massMod - mods[i].massPredictedOutput) <= dModTolerance)
	    return mods[i].mass;
        }
    }
  else
    {
      if (fabs(massMod - modsStatic[c].mass) <= dModTolerance ||
	  fabs(massMod - modsStatic[c].massPredictedOutput) <= dModTolerance)
	return modsStatic[c].mass;
    }

  return 0.0;
}

void TandemResultsParser::writePepXML(ostream& out, ScanData& scan)
{
  // If no search results, then skip this scan, per Jimmy Eng.
  if (scan.vScores.size() == 0)
    return;

  char buffer[256];
  ScoreData score;
  int num_tot_proteins = 0;

  vector<ScoreData>::iterator it = scan.vScores.begin();
  vector<ScoreData>::iterator end = scan.vScores.end();
  while (it != end)
    {
      if (score.hyperScore.empty())
	score = *it;

      //CONSIDER: Use a map to ensure not specified multiple times?
      if (it->isEqualPeptide(score))
	{
	  num_tot_proteins++;
	}
      else
	{
	  // Another sequence with top score should yield delta of 0.0.
	  score.nextScore = score.hyperScore;
	}
      it++;
    }

  string prot = score.prot;
  string prot_desc;

  size_t space = min(prot.find(' '), prot.find('\t'));
  if (space != prot.npos)
    {
#ifndef COMET_EXACT
      prot_desc = prot.substr(space + 1);
#endif
      prot.erase(space);
    }

  out << nl();
  out << "<spectrum_query spectrum=\"";
  // jmt: scan description may be blank; if so, try to assemble
  // one from other data
  bool gotscan = false;
  if (useDescription && scan.description != "") {
    size_t scanpos;
    if ((scanpos = scan.description.find("scan=")) != string::npos) {
      string scanstr = scan.description.substr(scanpos + 5);
      space = min(scanstr.find(' '), scanstr.find('\t'));
      if (space != scanstr.npos) {
	scanstr.erase(space);
      }
      scan.scanEnd = scan.scanStart = atol(scanstr.c_str());
      gotscan = true;
    }
  } 

  if (!gotscan && useDescription && scan.description != "") {
    out << scan.description.substr(0, scan.description.find_last_of("0123456789")+1);
  }
  else {
    // check to make sure we have the required values, leave blank otherwise
    string strMzXMLFile = paramsHandler.getUsedParam("spectrum, path");
    // remove all but the filename
    int pos = findRightmostPathSeperator(strMzXMLFile);
    if (pos != string::npos) {
      strMzXMLFile = strMzXMLFile.substr(pos+1);
    }
    if ((strMzXMLFile != "") && (scan.charge > 0) && (scan.id > 0)) {
      char *ext = rampValidFileType(strMzXMLFile.c_str()); // case insensitive
      if (ext != NULL) { // will point to .mzXML or .mzml or .mzdata or whatever
	pos = strlen(strMzXMLFile.c_str())-strlen(ext);
	strMzXMLFile = strMzXMLFile.substr(0,pos);
      }
      sprintf(buffer, "%s.%05d.%05d.%d", strMzXMLFile.c_str(), scan.scanStart, scan.scanEnd, scan.charge);
      out << buffer;
    }
    // else, it's blank
  }


  out << "\""
      << " start_scan=\"" << scan.scanStart << "\""
      << " end_scan=\"" << scan.scanEnd << "\"";
  sprintf(buffer, "%0.4f", scan.mh - dProtonMass);
  out				<< " precursor_neutral_mass=\"" << buffer << "\""
				<< " assumed_charge=\"" << scan.charge << "\""
				<< " index=\"" << scan.id << "\"";
  if (scan.retentionTime.length()) { // passed through by X!Tandem "Tornado" release
    string rt = scan.retentionTime;
    // don't output XML/XSD timestamp, just output float "time in seconds"
    if (rt.substr(0,2) == "PT" && rt[scan.retentionTime.length() - 1] == 'S') {
      rt=rt.substr(2, scan.retentionTime.length() - 3);
    }
    out       << " retention_time_sec=\"" << rt << "\"";
  } else // dig it out of mzXML file
    if(pf != NULL && scan.scanStart <= lastIndex) { // last index inclusive
      struct ScanHeaderStruct shs;
      readHeader(pf,mzXMLindices[scan.scanStart],&shs); 
      out                           << " retention_time_sec=\"" << shs.retentionTime << "\"";
    }
 
  out 				<< ">";

  if (scan.vScores.size() == 0)
    {
      out	<< nl();
    }
  else
    {
      char prevAA = *(score.seqPre.end() - 1);
      char nextAA = *(score.seqPost.begin());

      out << nl();
      out << "<search_result>" << nlIn();
      out      << "<search_hit hit_rank=\"1\" peptide=\"" << score.seq << "\""
	       << " peptide_prev_aa=\"" << *(score.seqPre.end() - 1) << "\""
	       << " peptide_next_aa=\"" << *(score.seqPost.begin()) << "\""
	       << " protein=\"" << XMLEscape(prot) << "\"";
      if (prot_desc.length() > 0)
	{
	  out  << " protein_descr=\"" << XMLEscape(prot_desc) << "\"";
	}

      string used_a_ions = paramsHandler.getUsedParam("scoring, a ions");
      string used_b_ions = paramsHandler.getUsedParam("scoring, b ions");
      string used_c_ions = paramsHandler.getUsedParam("scoring, c ions");
      string used_x_ions = paramsHandler.getUsedParam("scoring, x ions");
      string used_y_ions = paramsHandler.getUsedParam("scoring, y ions");
      string used_z_ions = paramsHandler.getUsedParam("scoring, z ions");

      int num_ion_series = (used_a_ions.compare("yes") == 0) + (used_b_ions.compare("yes") == 0) + (used_c_ions.compare("yes") == 0) + (used_x_ions.compare("yes") == 0) + (used_y_ions.compare("yes") == 0) + (used_z_ions.compare("yes") == 0);
      if (num_ion_series < 1) num_ion_series = 2; // failsafe

      out      << " num_tot_proteins=\"" << num_tot_proteins << "\""
	       << " num_matched_ions=\"" << (score.bIons + score.yIons + score.zIons + score.cIons + score.xIons + score.aIons) << "\""
	       << " tot_num_ions=\"" << (int)(score.seq.length() - 1) * num_ion_series * max(1, scan.charge - 1) << "\"";

      sprintf(buffer, "%0.4f", score.seqMh - dProtonMass);
      out      << " calc_neutral_pep_mass=\"" << buffer << "\"";
      sprintf(buffer, "%0.3f", scan.mh - score.seqMh);
      out      << " massdiff=\"" << buffer << "\""
	//     << " mass_diff=\"" << buffer << "\""
	       << " num_tol_term=\"" << sampleEnzyme->getNumTolTerm(prevAA, score.seq.data(), nextAA) << "\"";


/*
 * X!Tandem specified missed cleavages may not be compatible with
 * min_spacing parameter in the search_enzyme specificity tag.
 * Always use the searchEnzyme ProteolyticEnzyme instance for
 * num_missed_cleavages, since it produced the search_enzyme tag.
 *
 * ProteolyticEnzyme should produce a result equal to X!Tandem, if
 * min_spacing is set to 0.

		if (score.missedCleavages >= 0)
		{
			out				<< " num_missed_cleavages=\"" << score.missedCleavages << "\"";
		}
		else
		{
*/
      out      << " num_missed_cleavages=\"" << sampleEnzyme->getNumMissedCleavages(score.seq.data()) << "\"";
//		}
      out      << " is_rejected=\"" << 0 << "\"" // TODO: What is this?
	       << ">" << nlIn();

#ifndef COMET_EXACT
      // Add alternative proteins for this peptide
      it = scan.vScores.begin();
      while (it != end)
	{
	  if (it->prot != score.prot && it->isEqualPeptide(score))
	    {
	      string altProt = it->prot;
	      string altProt_desc;

	      space = min(altProt.find(' '), altProt.find('\t'));
	      if (space != altProt.npos)
		{
		  altProt_desc = altProt.substr(space + 1);
		  altProt.erase(space);
		}

	      out << "<alternative_protein protein=\"" << XMLEscape(altProt) << "\"";
	      if (altProt_desc.length() > 0)
		{
		  out <<	" protein_descr=\"" << XMLEscape(altProt_desc) << "\"";
		}
	      // Hypothetically possible for the same peptide in different proteins
	      // to differ in the number of properly cleaved termini.
	      prevAA = *(it->seqPre.end() - 1);
	      nextAA = *(it->seqPost.begin());
	      out <<      " num_tol_term=\"" << sampleEnzyme->getNumTolTerm(prevAA, it->seq.data(), nextAA) << "\"";
	      out << 		"/>" << nl();
	    }
	  it++;
	}
#endif

      // Make a copy of the modifications to keep track of what has
      // been written.
      vector<ModData> vMods;

      double savMass = 0.0;

      std::sort(score.vMods.begin(), score.vMods.end(), compare_mods);

      vector<ModData>::iterator it = score.vMods.begin();
      vector<ModData>::iterator sav_it = it;
      vector<ModData>::iterator prev_it = it;

      while (it != score.vMods.end())
        {
	  sav_it++;
	  if (it->aa == '[' || it->aa == ']') {
	    it->terminal = true;
	    //continue;
	  }
	  else if (it->aa == 'E' && fabs(it->modMass+h2o) < 0.0001) {
	    it->terminal = false;
	  }
	  else if (it->aa == 'Q' && fabs(it->modMass+nh3) < 0.0001) {
	    it->terminal = false;
	  }
	  else if ((long) modsStatic['C'].mass == 57 && 
		   it->aa == 'C' && fabs(it->modMass+nh3) < 0.0001) {
	    it->terminal = true;
	  }
	  else if ((score.seqStart == it->loc && (matchMod('[', it->modMass, false)!=0.0 || matchMod('[', it->modMass, true)!=0.0) ) || 
		   (score.seqEnd == it->loc &&  (matchMod(']', it->modMass, false)!=0.0 || matchMod(']', it->modMass, true)!=0.0) ) ) {
	    
	    if (prev_it == it || prev_it->loc != it->loc || !prev_it->terminal) {
	      it->terminal = true;
	    }
	    else if (prev_it->loc == it->loc && prev_it->modMass == it->modMass && prev_it->terminal) {
	      it->terminal = false;
	    }

	  }
	  else {
	    it->terminal = false;
	  }

	  vMods.push_back(*it);
	  prev_it = it;
	  it++;
        }

      // Terminal modifications - first add static delta
      double termNDiff = dNTerm - dNTermDefault;
      if (score.seqStart == 0)
	termNDiff += dNTermProt;
      double termCDiff = dCTerm - dCTermDefault;
      if (score.seqEnd == score.seq.length() - 1)
	termCDiff += dCTermProt;

      // Then look for terminal modifications
      for (int v = 0; v < 2; v++)
        {
	  // Static modifications (v == 0), and variable modifications (v == 1).
	  bool variable = (v == 1);

	  int locLast = -1;
	  it = vMods.begin();
	  while (it != vMods.end())
            {
	      int loc = it->loc;
	      double modMass = it->modMass;

	      // DDS: COMMENTED OUT - sum together variable term mods when found :::::: Allow only 1 static and 1 variable terminal modification.
	      //	      if (loc != locLast && loc == score.seqStart && ( matchMod('[', modMass, variable) != 0 ) )
	      if (loc == score.seqStart && ( matchMod('[', modMass, variable) != 0 ) )
                {
		  if (it->terminal) {
		    locLast = loc;
		    termNDiff += matchMod('[', modMass, variable);
		    it = vMods.erase(it);
		  }
		  else {
		    it++;
		  }
                }
	      //	      else if (loc != locLast && loc == score.seqEnd && ( matchMod(']', modMass, variable) != 0 ) )
	      else if (loc == score.seqEnd && ( matchMod(']', modMass, variable) != 0 ) )
                {
		  if (it->terminal) {
		    locLast = loc;
		    termCDiff += matchMod(']', modMass, variable);
		    it = vMods.erase(it);
		  }
		  else {
		    it++;
		  }
		}
	      else
                {
		  it++;
                }
            }
        }

      bool infoTagStarted = false;
      if (fabs(termNDiff) > dModTolerance)
        {
	  if (!infoTagStarted)
            {
	      infoTagStarted = true;
	      out		<< "<modification_info";
            }
	  sprintf(buffer, "%0.4f", termNDiff + dNTermDefault);
	  out << " mod_nterm_mass=\"" << buffer << "\"";
        }
      if (fabs(termCDiff) > dModTolerance)
        {
	  if (!infoTagStarted)
            {
	      infoTagStarted = true;
	      out		<< "<modification_info";
            }
	  sprintf(buffer, "%0.4f", termCDiff + dCTermDefault);
	  out << " mod_cterm_mass=\"" << buffer << "\"";
        }
      int nAAMods = 0;
      it = vMods.begin();
      while (it != vMods.end())
	{
	  int loc = it->loc;
	  char aa = makeupperc(it->aa);

	  // Sum all modifictions to this location, ensuring validity.
	  double modMass = 0.0;
	  savMass = 0.0;
	  double modMassStatic = 0.0;
	  int nVMods = 0;
	  while (it < vMods.end() && loc == it->loc)
	    {
	      if (aa == makeupperc(it->aa))
		{
		  // If no static mod found yet, check for a static mod.
		  if (modMassStatic == 0.0 && matchMod(aa, it->modMass, false) != 0.0)
		    {
		      modMassStatic = matchMod(aa, it->modMass, false);
		      modMass += modMassStatic;
		      it = vMods.erase(it);
		      continue;
		    }
		  // Then check for variable.
		  else if (matchMod(aa, it->modMass, true) != 0.0)
		    {
		      if (nVMods++ == 1)
			{
			  cerr << "WARNING: Multiple variable modifications on "
			       << it->aa << " for scan " << scan.scanStart << ".\n";
			}
		      modMass += matchMod(aa, it->modMass, true);
		      it = vMods.erase(it);
		      continue;
		    }
		}
              sav_it = it;  
	      savMass += it->modMass;
	      it++;   // Reported later as unmatched
	    }

	  if (modMass != 0.0)
	    {
	      if (nAAMods++ > 0)
		out << nl();
	      else		{
		  if (!infoTagStarted)
		    {
		      infoTagStarted = true;
		      out		<< "<modification_info";
		    }
		  out     << ">" << nlIn();
		}
	      out		<< "<mod_aminoacid_mass"
				<< " position=\"" << loc - score.seqStart + 1 << "\"";	// 1-based index
	      sprintf(buffer, "%0.4f", modMass + dMassAA[aa]);
	      out			<< " mass=\"" << buffer << "\""
					<< " />";
	    }
	  else {
	    modMass = savMass;
	    if (modMass != 0.0 ) {
	      aa = sav_it->aa;
	      if (nAAMods++ > 0)
		out << nl();
	      else
		{
		  if (!infoTagStarted)
		    {
		      infoTagStarted = true;
		      out		<< "<modification_info";
		    }
		  out     << ">" << nlIn();
		}
	      out		<< "<mod_aminoacid_mass"
				<< " position=\"" << loc - score.seqStart + 1 << "\"";	// 1-based index
	      sprintf(buffer, "%0.4f", modMass + dMassAA[aa]);
	      out			<< " mass=\"" << buffer << "\"";
		
	      if (sav_it->pm != '\0') {
		out			<< " alt_aa=\"" << sav_it->pm << "\"";
	      }
	      
	      out			<< " />";
	      
	    }
	  }
	}

      // Anything that hasn't yet matched causes a warning.
      it = vMods.begin();
      while (it != vMods.end())
	{
	  cerr << "WARNING: Unknown modification '" << it->aa
	       << "' (" << it->modMass << ") for scan " << scan.scanStart << ".\n";
	  it++;
	}

      if (infoTagStarted)
	{
	  if (nAAMods == 0)
	    out << " />" << nl();
	  else
	    {
	      out << nlOut();
	      out << "</modification_info>" << nl();
	    }
	}

      if (m_bCometScoring)
	{
	  double delta = 1.0;
	  int score1 = atoi(score.hyperScore.data());
	  int score2 = atoi(score.nextScore.data());
	  delta -= ((double) score2) / ((double) score1);

	  char buffer[20];
	  sprintf(buffer, "%.03f", delta);
	  out		<< "<search_score name=\"dotproduct\" value=\"" << score.hyperScore << "\" />" << nl()
			<< "<search_score name=\"delta\" value=\"" << buffer << "\" />" << nl()
			<< "<search_score name=\"deltastar\" value=\"0\" />" << nl()
			<< "<search_score name=\"zscore\" value=\"" << 0 << "\" />" << nl()
			<< "<search_score name=\"expect\" value=\"" << score.expect << "\" />";
	}
      else
	{
	  out		<< "<search_score name=\"hyperscore\" value=\"" << score.hyperScore << "\"/>" << nl()
			<< "<search_score name=\"nextscore\" value=\"" << score.nextScore << "\"/>" << nl()
			<< "<search_score name=\"bscore\" value=\"" << score.bScore << "\"/>" << nl()
			<< "<search_score name=\"yscore\" value=\"" << score.yScore << "\"/>" << nl()
			<< "<search_score name=\"cscore\" value=\"" << score.cScore << "\"/>" << nl()
			<< "<search_score name=\"zscore\" value=\"" << score.zScore << "\"/>" << nl()
			<< "<search_score name=\"ascore\" value=\"" << score.aScore << "\"/>" << nl()
			<< "<search_score name=\"xscore\" value=\"" << score.xScore << "\"/>" << nl()
			<< "<search_score name=\"expect\" value=\"" << score.expect << "\"/>";
	}

      out		<< nlOut();
      out		<< "</search_hit>";
      out		<< nlOut();
      out		<< "</search_result>" << nl();
    }

  out	<< "</spectrum_query>";
}
