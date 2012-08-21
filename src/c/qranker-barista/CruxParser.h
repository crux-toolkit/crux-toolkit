/*************************************************************************//**
 * \file CruxParser.cpp
 * AUTHOR: Manijeh Naseri and Sean McIlwain
 * CREATE DATE: 06/06/2012
 * \brief Parses tab-delimited search result files
 ****************************************************************************/
#ifndef CRUXPARSER_H
#define CRUXPARSER_H
#define CRUX

#include <sys/stat.h>
#ifndef _MSC_VER
#include <dirent.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <string>
#include <math.h>
#include <map>
#include <set>
#include <cstring>
#include "SpecFeatures.h"
#include "BipartiteGraph.h"
#include "SQTParser.h"
#ifdef CRUX
#include "CruxApplication.h"
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#endif
using namespace std;

class CruxParser:public SQTParser{
 public:

  CruxParser();
  virtual ~CruxParser();

 
/**
 * Parse tab delimited file.
 * Generates the same QRanker internal tables.
 * Set the matches in the tab delimited file. 
 */
  void readMatches(
    MatchFileReader& reader,///<Reader for the delimted file.
    int final_hits,///<Total number of matches
    enzyme enz, ///<Enzyme in used on search 
    bool decoy ///< Are all the matches decoy?
  );

  /*
  *gets the path of delimited file 
  *\returns true if it can open the file 
  */
  virtual bool read_search_results(
    string& cur_fname, ///< current delimited file path tp parse 
    bool decoy ///< Are all the matches decoy?
  ); 
  virtual string  get_parser_extension(); 


 protected:
  ofstream f_pepind2flanking_aa;
  ofstream f_psmind2xcorr_rank;
  ofstream f_psmind2match_spectrum;
  ofstream f_psmind2cleavage_type; 
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
