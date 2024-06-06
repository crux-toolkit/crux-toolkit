#ifndef XPRESS_GROUP_PEP_PARSER_H
#define XPRESS_GROUP_PEP_PARSER_H

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Quantitation/searchHitCache.h"

#include "Parsers/Parser/Parser.h"
#include "Quantitation/Option.h"
#include "Validation/MixtureModel/MixtureModel.h"
#include "Parsers/mzParser/mzParser.h"
#include "Parsers/Parser/TagFilter.h"

#define MAX_LOG 3

// things that are interesting to remember about a peptide search hit
typedef struct {
   double heavy_;
   double light_;
   double probability_;
} XPressRatioSearchHit;

typedef searchHitCache<XPressRatioSearchHit> XPressRatioSearchHitCache;

class XPressGroupPeptideParser {

 public:

  XPressGroupPeptideParser(const XPressRatioSearchHitCache &peptideParser, const peplist* peptides, double minprob/*, Boolean heavy2light*/);
  ~XPressGroupPeptideParser();
  void setFilter(Tag* tag);
  double getRatioLogSum();
  double getRatioSquareSum();
  int getRatioNum();
  RatioStruct getRatio();

 protected:

  void parse();
  Boolean peptideListMember(const char* pep);


  const peplist* peptides_;
  double min_probability_;
  //double ratio_sum_;

  double h2l_ratio_square_sum_;
  double h2l_ratio_log_sum_;
  
  double ratio_square_sum_;
  double ratio_log_sum_;

  //double inverse_ratio_sum_;
  //double inverse_ratio_square_sum_;
  int ratio_num_;
  //Boolean heavy2light_;

  const XPressRatioSearchHitCache &searchHits_;

};



#endif
