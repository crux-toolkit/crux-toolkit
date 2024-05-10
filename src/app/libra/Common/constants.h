#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "sysdepend.h"

#define SIZE_BUF             8192
#define SIZE_FILE            1024
#define SIZE_PEP             128
#define INIT_INTERACT_LINES  1000
#define PERCENTAGE           0.75
#define EXTRAITRS            5

#define USE_LOCAL_TIME 1

#define MAX_CHARGE 7

//
// CONFIGURATION DEFINITIONS
//
// Provide default values for the filesystem/url configuration which may be
// overridden during compile.
//

// Path to the TPP directory, or "home".  Paths such as cgi-bin, html, conf
// are all relative to this path.  
#ifndef TPP_HOME
#if defined(_MSC_VER) || defined(__MINGW32__)
#define TPP_HOME "C:/TPP"
#else
#define TPP_HOME "/usr/local/tpp"
#endif
#endif

// TPP path to the user data.
#ifndef TPP_DATADIR
#define TPP_DATADIR TPP_HOME"/data"
#endif

#ifndef TPP_PORT
#if defined(_MSC_VER) || defined(__MINGW32__)
#define TPP_PORT "10401"
#else
#define TPP_PORT "80"
#endif
#endif


// Prefix used for all TPP URL's.  Not prefixed with a '/' in case its being
// overridden at compile time from MinGW/MSYS which mangles paths
#ifndef TPP_BASEURL
#define TPP_BASEURL "tpp/"
#endif

// Prefix used for TPP URL to data. Not prefixed with a '/' in case its being
// overridden at compile time from MinGW/MSYS which mangles paths
#ifndef TPP_DATAURL
#define TPP_DATAURL "tpp/data"
#endif

//
// Older pepxml path definitions 
//

// Linux
#if !(defined(_WIN32))

// schema
#define DEFAULT_PEPXML_STD_XSL TPP_HOME"/schema/"  // a webserver reference, to the directory, filename is hard coded

// schema namespace constants
#define PEPXML_NAMESPACE "http://regis-web.systemsbiology.net/pepXML"

// Win32 Installation
#else

// these are for pepXML specific functionality
#define PEPXML_NAMESPACE "http://regis-web.systemsbiology.net/pepXML"
#define DEFAULT_PEPXML_STD_XSL "http://localhost/"
#define DEFAULT_PEPXML_STD_XSL_WEB_PATH "http://localhost/"

#endif


//
// pepXML stuff that's platform independent
//
#define DEFAULT_PEPXML_FILENAME_DOTEXT ".pep.xml" // was .xml before Jan 2008, can override with env var PEPXML_EXT
#define PEPXML_NAMESPACE_PX "pepx"
#define PEPXML_SCHEMA "pepXML_v123.xsd"

//
// protXML stuff that's platform independent
//
#define DEFAULT_PROTXML_FILENAME_DOTEXT ".prot.xml" // was -prot.xml before Jan 2008,, can override with env var PROTXML_EXT
#define PROTXML_SCHEMA "protXML_v10.xsd"

#define _ASAPRATIO_MXQ_ 5
#define CON_SIZE_FILE        256
#define CON_SIZE_PATH       4096
#define LABEL_MASS_SIZE 1000
#define AMINO_ACIDS "ARNDCEQGHILKMFPSTWYV"
#define LINE_WIDTH 1000000
#define USE_STD_MODS 1  // comment this out if don't use standard way of representing modifications
#define MOD_ERROR 0.5
# define _ISOMASSDIS_ 1.0033548378

//Comet constant (moved from Comet.h)
#define COMETLINKSFILE  "cometlinks.def"

struct RatioStruct
{
   int iNumPeptides;
   double dRatio;
   double dStdDev;
   double dh2lRatio;
   double dh2lStdDev;
};

#include <string.h>

class InputStruct
{
public:
   InputStruct() {
      memset(this,0,sizeof(InputStruct));
   }
   int bPpmMassTol;               // command line option
   int bForceMS1Endpoints;        // command line option
   int iAnalysisFirstScan;
   int iAnalysisLastScan;
   int iScanCount;                // scan count; typically but not necessarily iAnalysisLastScan+1
   int bZeroAllBackGrnd;          // command line option
   int bQuantHighBackGrnd;        // command line option
   int bUseSameScanRange;         // command line option
   int bUseFixedScanRange;        // command line option (not really a bool, but a 3 way flag)
   int bIntensityRatio;           // command line option
   int bLabelFreeMode;
   int iFixedScanRange;           // command line option
   int bQuantCIDChrgOnly;         // command line option
   int bXpressLight1;             // command line option
   int bUseWaveletSmoothing;      // command line option, currently experimental for ASAPRatioPeptideParser
   int iMetabolicLabeling;
   int iMinNumChromatogramPoints; // command line option
   int iMinNumIsotopePeaks;       // command line option
   int iStartCharge;              // -C label free command line option for xpress
   int iEndCharge;                // -C label free command line option for xpress
   char szXpressResidues1[256];   // command line option
   char szXpressResidues2[256];   // command line option
   char szXpressResidues3[256];   // command line option
   char szXpressResidues[256];    // total of above
   double dXpressMassDiff;        // TO BE DEPRECATED AFTER ASAP BROUGHT ONBOARD
   double dXpressMassDiff1;       // command line option
   double dXpressMassDiff2;       // command line option
   double dXpressMassDiff3;       // command line option
   double dMassTol;               // command line option

   double dMinPprob;
   double dMinIprob;

   char szMzXMLDir[CON_SIZE_PATH]; // command line option

   int iFirstScan;             // from xml
   int iLastScan;              // from xml
   int iChargeState;           // from xml
   int iModCount;              // from xml, for metabolic labeling count # of mods in peptide to determine if light or heavy
   char szPeptide[128];        // from xml
   double dPeptideMass;        // from xml Measured Mass + H
   double dCalcPeptideMass;        // from xml Search Computed
   char szOutFile[CON_SIZE_FILE];  // from xml
   char szSpectrumName[CON_SIZE_FILE];  // from xml
   bool staticQuant; // whether or not quant is based on static variables
   char labelMasses[LABEL_MASS_SIZE];
};

struct lcPeakStrct {
  int indx;
  int peak;
  int valley[2];
  double bckgrnd;
  double area[2];
  double time[2];
};


struct pepDataStrct {
  int indx;
  long scan;
  int chrg;
  int cidIndx;
  double msLight;
  double msHeavy;
  int eltn;
  int areaFlag;
  struct lcPeakStrct peaks[_ASAPRATIO_MXQ_][2];
  double pkRatio[_ASAPRATIO_MXQ_];
  double pkError[_ASAPRATIO_MXQ_];
  int pkCount[_ASAPRATIO_MXQ_];
  double pepRatio[2];
  double pepH2LRatio[2];
  double pepTime[2][2];
  double pepArea;
};


struct pairedLabel {
  char labelA[3];
  char labelB[3];
};

struct Pep3D_dataStrct {
  double prob;
  double score;
  int startScan;
  int endScan;
  int charge;
};

class spectStrct; // defined in spectStrct.h

#endif
