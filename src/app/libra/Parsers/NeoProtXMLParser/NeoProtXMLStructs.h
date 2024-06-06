#ifndef NEOPROTXMLSTRUCTS_H
#define NEOPROTXMLSTRUCTS_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>

static std::string xmlns="http://regis-web.systemsbiology.net/protXML";
static std::string xmlns_xsi="http://www.w3.org/2001/XMLSchema-instance";
static std::string xsi_schemaLocation="http://sashimi.sourceforge.net/schema_revision/protXML/protXML_v8.xsd";

enum protXMLElement :int{
  prAffectedChannel,
  prAnalysisResult,
  prAnalysisSummary,
  prAnnotation,
  prContributingChannel,
  prDatasetDerivation,
  prDecoyAnalysis,
  prDecoyAnalysisSummary,
  prErrorPoint,
  prFragmentMasses,
  prIndistinguishablePeptide,
  prIndistinguishableProtein,
  prIntensity,
  prIsotopicContributions,
  prLibraResult,
  prLibraSummary,
  prModAminoacidMass,
  prModificationInfo,
  prNSPInformation,
  prNSPDistribution,
  prParameter,
  prPeptide,
  prPeptideParentProtein,
  prPoint,
  prProgramDetails,
  prProtein,
  prProteinGroup,
  prProteinProphetDetails,
  prProteinSummary,
  prProteinSummaryDataFilter,
  prProteinSummaryHeader,
  prStPeterAnalysisSummary,
  prStPeterQuant,
  prStPeterQuantPeptide,
  PROTXML_NUM_ELEMENTS
};


static void NPRprintTabs(FILE* f, int tabs){
  for(int i=0;i<tabs;i++) fprintf(f," ");
}

static void NPRerrMsg(std::string el, std::string attr){
  std::cerr << el << "::" << attr << " required." << std::endl;
  exit(69);
}

typedef struct nprDate {
  int year;
  int month;
  int day;
  nprDate() {
    year = 0;
    month = 0;
    day = 0;
  }
} nprDate;

typedef struct nprTime {
  int hour;
  int minute;
  int second;
  nprTime() {
    hour = 0;
    minute = 0;
    second = 0;
  }
} nprTime;

typedef struct nprDateTime{
  nprDate date;
  nprTime time;
  void clear(){
    date.day = 0;
    date.month = 0;
    date.year = 0;
    time.hour = 0;
    time.minute = 0;
    time.second = 0;
  }
  void parseDateTime(const char* dt){
    if (strlen(dt)<2){
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
    s = str;
    return s;
  }
} nprDateTime;

#endif
