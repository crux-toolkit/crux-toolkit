#ifndef NEOPEPXMLSTRUCTS_H
#define NEOPEPXMLSTRUCTS_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>

//For Windows
#ifdef _MSC_VER
#define __inline__ __inline
typedef _int64  __int64_t;
typedef unsigned _int32 uint32_t;
typedef unsigned _int64 uint64_t;
typedef __int64 f_off;
#define npxfseek(h,p,o) _fseeki64(h,p,o)
#define npxftell(h) _ftelli64(h)
#define npxatoi64(h) _atoi64(h)
#endif

//For MinGW toolset, which lacks the ftello, fseeko, etc functions
#ifdef MINGW
//typedef __int64 f_off;
//#define __int64_t int64_t
#define npxfseek(h,p,o) fseeko64(h,p,o)
#define npxftell(h) ftello64(h)
#define npxatoi64(h) _atoi64(h)
//#include <stdexcept>
#endif

#if defined(GCC) || defined(__LINUX__) || defined(GNUC) || defined(__MINGW32__)
#include <stdint.h>
#include <stdexcept>
#ifndef _LARGEFILE_SOURCE
#error "need to define _LARGEFILE_SOURCE!!"
#endif    /* end _LARGEFILE_SOURCE */
#if _FILE_OFFSET_BITS<64
#error "need to define _FILE_OFFSET_BITS=64"
#endif

typedef off_t f_off;
#define npxfseek(h,p,o) fseeko(h,p,o)
#define npxftell(h) ftello(h)
#define npxatoi64(h) atoll(h)
#endif

#ifdef OSX
#define __inline__ inline
#ifndef OSX_TIGER
#define __int64_t int64_t
#endif
#endif

// this define for the INTEL-based OSX platform is untested and may not work
#ifdef OSX_INTEL
#define __inline__ inline
#endif

static std::string npx_xmlns = "http://regis-web.systemsbiology.net/pepXML";
static std::string npx_xmlns_xsi = "http://www.w3.org/2001/XMLSchema-instance";
static std::string npx_xsi_schemaLocation = "http://regis-web.systemsbiology.net/pepXML /tools/bin/TPP/tpp/schema/pepXML_v123.xsd";

enum pepXMLElement:int{
  pxAffectedChannel,
  pxAlternativeProtein,
  pxAminoAcidModification,
  pxAminoAcidSubstitution,
  pxAnalysisResult,
  pxAnalysisSummary,
  pxAnalysisTimestamp,
  pxBin,
  pxContributingChannel,
  pxDatabaseRefreshTimestamp,
  pxDatasetDerivation,
  pxDecoyAnalysis,
  pxDecoyAnalysisSummary,
  pxDistributionPoint,
  pxEnzymaticSearchConstraint,
  pxErrorPoint,
  pxFragmentMasses,
  pxInputfile,
  pxIntensity,
  pxInteractSummary,
  pxInterprophetResult,
  pxInterprophetSummary,
  pxIsotopicContributions,
  pxLability,
  pxLibraResult,
  pxLibraSummary,
  pxLinkedPeptide,
  pxMixture_Model,
  pxMixturemodel,
  pxMixturemodelDistribution,
  pxModAminoAcidMass,
  pxModAminoAcidProbability,
  pxModificationInfo,
  pxModTerminalProbability,
  pxMSMSPipelineAnalysis,
  pxMSMSRunSummary,
  pxNegmodelDistribution,
  pxParameter,
  pxPeptideProphetResult,
  pxPeptideprophetSummary,
  pxPepXMLQuantResult,
  pxPoint,
  pxPosmodelDistribution,
  pxPTMProphetResult,
  pxPTMProphetSummary,
  pxQuanticResult,
  pxQuanticSummary,
  pxROCDataPoint,
  pxROCErrorData,
  pxSampleEnzyme,
  pxSearchDatabase,
  pxSearchHit,
  pxSearchResult,
  pxSearchScore,
  pxSearchScoreSummary,
  pxSearchSummary,
  pxSpecificity,
  pxSpectrumQuery,
  pxTerminalModification,
  pxXLink,
  pxXLinkScore,
  pxXpressLabelFreeResult,
  pxXpressLabelFreeSummary,
  PEPXML_NUM_ELEMENTS
};

typedef struct npxDate {
  int year;
  int month;
  int day;
  npxDate() {
    year = 0;
    month = 0;
    day = 0;
  }
} npxDate;

typedef struct npxTime {
  int hour;
  int minute;
  int second;
  npxTime() {
    hour = 0;
    minute = 0;
    second = 0;
  }
} npxTime;

typedef struct npxDateTime{
  npxDate date;
  npxTime time;
  void clear(){
    date.day=0;
    date.month=0;
    date.year=0;
    time.hour=0;
    time.minute=0;
    time.second=0;
  }
  void parseDateTime(const char* dt){
    if(strlen(dt)<2){
      clear();
      return;
    }
    int x = sscanf(dt, "%d-%d-%dT%d:%d:%d", &date.year, &date.month, &date.day, &time.hour, &time.minute, &time.second);
  }
  void parseDateTime(std::string s){
    parseDateTime(s.c_str());
  }
  std::string write(){
    std::string s;
    char str[64];
    sprintf(str, "%4d-%02d-%02dT%02d:%02d:%02d", date.year, date.month, date.day, time.hour, time.minute, time.second);
    s=str;
    return s;
  }
} npxDateTime;

static void NPXerrMsg(std::string el, std::string attr){
  std::cerr << el << "::" << attr << " required." << std::endl;
  exit(69);
}

static void NPXprintTabs(FILE* f, int tabs){
  for (int i = 0; i<tabs; i++) fprintf(f, " ");
}

#endif
