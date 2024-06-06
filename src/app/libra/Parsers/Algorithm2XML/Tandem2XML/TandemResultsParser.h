// TandemResultsParser.h
//     Handles parsing of X!Tandem output file, and writing
//     to pepXML.

#include <iostream>
#include <string>
#include <vector>
#include <set>

using namespace std;

//#include "TandemResultsParser.h"
#include "Parsers/Algorithm2XML/saxtandemhandler.h"
#include "TandemParamsParser.h"
#include "Parsers/mzParser/mzParser.h"

struct InstrumentStruct;
class ProteolyticEnzyme;

class TandemResultsParser : public TANSAXHandler {
 public:
  TandemResultsParser();
  virtual ~TandemResultsParser();

  virtual void startElement(const XML_Char *el, const XML_Char **attr);
  virtual void endElement(const XML_Char *el);
  virtual void characters(const XML_Char *s, int len);	

  bool writePepXML();

  inline void setIndexOffset(int off) {
    indexOffset = off;
  }

  inline void setUseDescription(bool useDesc) {
    useDescription = useDesc ;
  }
  
  inline void setUseN15(bool n15) {
    useN15 = n15 ;
  }

  inline void setOutputFile(const char* fileName) {
    m_strOutputFile = fileName;
  }

  inline void setSampleEnzyme(const char* enzyme) {
    sampleEnzymeName = enzyme;
  }

  //protected:
  class ModData {
  public:
    ModData() {
      clear();
    }

    ModData(const ModData& rhs) {
      *this = rhs;
    }

    void clear() {
      aa = '\0';
      pm = '\0';
      loc = -1;
      modMass = 0.0;
      terminal = false;
    }

    ModData& operator=(const ModData& rhs) {
      pm = rhs.pm;
      aa = rhs.aa;
      loc = rhs.loc;
      modMass = rhs.modMass;
      terminal = rhs.terminal;

      return *this;
    }
    char pm;
    char aa;
    int loc;
    double modMass;
    bool terminal;
  };
	
  class ModSpecData {
  public:
    ModSpecData() {
      clear();
    }

    ModSpecData(const ModSpecData& rhs) {
      *this = rhs;
    }

    void clear() {
      aa = '\0';
      symbol = '\0';
      mass = 0.0;
      massPredictedOutput = 0.0;
      comment.clear();
    }

    ModSpecData& operator=(const ModSpecData& rhs) {
      aa = rhs.aa;
      symbol = rhs.symbol;
      mass = rhs.mass;
      massPredictedOutput = rhs.massPredictedOutput;
      comment = rhs.comment;

      return *this;
    }

    char aa;
    char symbol;
    double mass;
    double massPredictedOutput;
    string comment;
  };

  static bool lessThanModLoc(const ModData& l, const ModData& r) {
    return (l.loc < r.loc);
  }

  class ScoreData {
  public:
    ScoreData() {
      clear();
    }

    ScoreData(const ScoreData& rhs) {
      *this = rhs;
    }

    void clear() {
      prot.clear();
      expect = 0;
      hyperScore.clear();
      nextScore.clear();
      seq.clear();
      seqPre.clear();
      seqPost.clear();
      seqMh = 0.0;
      yScore = 0.0;
      yIons = 0;
      bScore = 0.0;
      bIons = 0;
      zScore = 0.0;
      zIons = 0;
      cScore = 0.0;
      cIons = 0;
      xScore = 0.0;
      xIons = 0;
      aScore = 0.0;
      aIons = 0;
      missedCleavages = -1;
      vMods.clear();
    }

    ScoreData& operator=(const ScoreData& rhs) {
      prot.assign(rhs.prot);
      expect = rhs.expect;
      hyperScore.assign(rhs.hyperScore);
      nextScore.assign(rhs.nextScore);
      seq.assign(rhs.seq);
      seqPre.assign(rhs.seqPre);
      seqPost.assign(rhs.seqPost);
      seqStart = rhs.seqStart;
      seqEnd = rhs.seqEnd;
      seqMh = rhs.seqMh;
      yScore = rhs.yScore;
      yIons = rhs.yIons;
      bScore = rhs.bScore;
      bIons = rhs.bIons;
      zScore = rhs.zScore;
      zIons = rhs.zIons;
      cScore = rhs.cScore;
      cIons = rhs.cIons;
      xScore = rhs.xScore;
      xIons = rhs.xIons;
      aScore = rhs.aScore;
      aIons = rhs.aIons;
      missedCleavages = rhs.missedCleavages;

      vMods.clear();
      vector<ModData>::const_iterator it = rhs.vMods.begin();
      vector<ModData>::const_iterator end = rhs.vMods.end();
      while (it != end) {
	vMods.push_back(*it);
	it++;
      }

      return *this;
    }

    bool isEqualPeptide(ScoreData& score2) {
      if (seq.compare(score2.seq) != 0)
	return false;
      if (vMods.size() != score2.vMods.size())
	return false;

      vector<ModData>::iterator it = vMods.begin();
      vector<ModData>::iterator it2 = score2.vMods.begin();
      while (it != vMods.end() && it2 != score2.vMods.end()) {
	if (it->aa != it2->aa)
	  return false;
	if (it->loc - seqStart != it2->loc - score2.seqStart)
	  return false;
	it++;
	it2++;
      }

      return true;
    }

    string prot;
    double expect;
    string hyperScore;	// Keep original formatting
    string nextScore; // Keep original formatting
    string seq;
    string seqPre;
    string seqPost;
    int seqStart;
    int seqEnd;
    double seqMh;
    double yScore;
    int yIons;
    double bScore;
    int bIons;
    double zScore;
    int zIons;
    double cScore;
    int cIons;
    double xScore;
    int xIons;
    double aScore;
    int aIons;
    int missedCleavages;

    vector<ModData> vMods;
  };

  static bool lessThanScoreMods(const ScoreData& l, const ScoreData& r) {
    vector<ModData>::const_iterator lit = l.vMods.begin();
    vector<ModData>::const_iterator rit = r.vMods.begin();

    for (;;) {
      if (lit == l.vMods.end())
	return false;
      else if (rit == r.vMods.end())
	return true;
      else {
	int loc = lit->loc - l.seqStart;
	int roc = rit->loc - r.seqStart;
	if (loc != roc)
	  return (loc < roc);
      }
      lit++;
      rit++;
    }
  }

  class ScanData {
  public:
    ScanData() {
      clear();
    }

    void clear() {
      mh = 0.0;
      charge = 0;
      id = 0;
      prot.clear();
      vScores.clear();
    }

    double mh;
    int charge;
    int id;
    int scanStart;
    int scanEnd;
    string description;
    string retentionTime;

    string prot;

    vector<ScoreData> vScores;
  };

 protected:
  bool writePepXML(ostream& out, const char *baseName, const char *rawDataXMLext, const char *fragmentType);
  void writePepXML(ostream& out, ScanData& scan);

  bool readMods(string mod, bool variable, int pass);
  bool readMotifs(string motif);
  void writeModSpecs(ostream& out);
  void writeAAModSpec(ostream& out, const ModSpecData& mod, bool variable);
  void writeTermModSpec(ostream& out, const ModSpecData& mod, bool variable, bool peptide);

  void addModSpec(ModSpecData mod, bool variable);
  double matchMod(char c, double massMod, bool variable);

 protected:
  void startGroup(const XML_Char **attr);
  void endGroup();
  void startProtein(const XML_Char **attr);
  void endProtein();
  void startDomain(const XML_Char **attr);
  void endDomain();
  void startAA(const XML_Char **attr);
  void startNote(const XML_Char **attr);
  void endNote();

 protected:
  TandemParamsParser paramsHandler;

  mzParser::RAMPFILE *pf;
  mzParser::ramp_fileoffset_t *mzXMLindices;
  int lastIndex;

  mzParser::InstrumentStruct* instrument;
  string sampleEnzymeName;
  ProteolyticEnzyme* sampleEnzyme;
  string searchEnzymeName;
  ProteolyticEnzyme* searchEnzyme;

  bool m_bCometScoring;

  string m_strOutputFile;
  ostream* m_pout;

  double dNTermDefault;
  double dCTermDefault;
  double dNTerm;
  double dCTerm;
  double dNTermProt;
  double dCTermProt;

  bool useDescription;
  bool useN15;

  double h2o = 18.0105647;
  double nh3 = 17.02655;
  
  double dMassAA[256];
  ModSpecData modsStatic[256];
  vector<ModSpecData> vModsVariable[256];

  size_t scanCount;
  ScanData scanCur;
  ScoreData scoreCur;

  int indexOffset;
  int groupDepth;
  const char* state;
};

