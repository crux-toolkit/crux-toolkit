#ifndef COMBINEOUT_H_
#define COMBINEOUT_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "Common/tpp_hashmap.h"
#ifndef _MSC_VER
#include <unistd.h>
#include <dirent.h>
#endif

#include <string.h>

#include "Parsers/Algorithm2XML/SearchResult/SequestResult.h"
#include "Common/sysdepend.h"
#include "Common/Array.h"
#include "SequestParams.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Common/ModifiedResidueMass/ModifiedResidueMass.h"
#include "Parsers/Algorithm2XML/pICalculator.h"
#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#include "SequestOut.h"
#include "SequestHit.h"
#include "Parsers/Parser/Parser.h"
#include "Common/constants.h"
#include "Common/util.h"
#include "Out2XML.h"

using namespace std;

#include "Common/tpp_hashmap.h"  // defines TPP_CONSTCHARP_HASHMAP

typedef TPP_CONSTCHARP_HASHMAP(SequestOut*) SequestOutPtrMap;

class CombineOut : public Out2XML {
  public:
  CombineOut();
  CombineOut(char* inpath1, char* inpath2, char* outpath);
  ~CombineOut();
  void combine();
  SequestOut* mergeOuts(SequestOut* out1, SequestOut* out2);
  void updateHits(SequestOut* out);
  void processData();
  void writeData();
  void writeParams();
  void writeOutFiles();
  void combineParams();
  //  void readOutFile(char* fileName, SequestOut* data, struct HeaderStruct * hdr);
  void readOutDir(char* path,  SequestOutPtrMap& outFiles, SequestParams*& seqParam);
  void readTgzFile(char* path, SequestOutPtrMap& outFiles, SequestParams*& seqParam);
 private:
  FILE* paramsFile_;
  char* paramsPath_;
  Boolean tgzfile1_;
  Boolean tgzfile2_;
  char inpath1_[4086];
  char inpath2_[4086];
  char outpath_[4086];
  SequestParams* seqParams1_;
  SequestParams* seqParams2_;
  SequestOutPtrMap outFiles1_; 
  SequestOutPtrMap outFiles2_;
  SequestOutPtrMap combOutFiles_;
  TPP_HASHMAP_T<char, char> newModsfromStatic_;
  TPP_HASHMAP_T<char, TPP_HASHMAP_T<char, char>* > newModsfromVariable_;

};


#endif
