#ifndef _KLOG_H
#define _KLOG_H

#include <cstdint>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#define LOGSZ 500

//Warning Codes:
//0 = MS/MS Spectrum is centroid, but parameter is profile

typedef struct kWarning {
  int count;
  std::string msg;
} kWarning;

class KLog {
public:
  KLog();
  
  void addError(std::string msg);
  void addDBWarning(std::string msg);
  void addMessage(std::string msg, bool silent=false);
  void addParameter(std::string msg);
  void addParameterWarning(std::string msg);
  void addWarning(size_t id,std::string msg);
  
  void clear();
  void exportLog();

  void setDBinfo(std::string fn, int prot, int pep, int linkPep);
  void setLog(char* fn);
  void setLog(std::string fn);

private:
  size_t idIndex[LOGSZ]; //room for defined number of unique warning messages
  
  std::string dbInfo;
  std::string logFile;
  std::vector<std::string> vDBWarnings; //special log for database;
  std::vector<std::string> vMsg;
  std::vector<std::string> vParams; //special log for params;
  std::vector<std::string> vParamWarnings; //special log for params;
  std::vector<kWarning> vWarnings;
  std::string strError;


};

#endif
