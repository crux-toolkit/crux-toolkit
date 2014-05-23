#ifndef _CHARDKLORPARSER_H
#define _CHARDKLORPARSER_H

#include <fstream>
#include <iostream>
#include <vector>

#include "CHardklorSetting.h"
#include "CPeriodicTable.h"
#include "MSToolkitTypes.h"

using namespace std;

class CHardklorParser {

 public:
  //Constructors & Destructors
  CHardklorParser();
  ~CHardklorParser();

  //Methods
  void parse(char*);
  void parseConfig(char*);
	MSFileFormat getFileFormat(char* c);
  CHardklorSetting& queue(int);
  int size();

 protected:

 private:
  //Data Members
  CHardklorSetting global;
  vector<CHardklorSetting> *vQueue;
};


#endif
