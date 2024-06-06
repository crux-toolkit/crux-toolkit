#ifndef Q3_GROUP_PEP_PARSER_H
#define Q3_GROUP_PEP_PARSER_H

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
   double q2_heavy_;
   double q2_light_;
   double probability_;
} Q3RatioSearchHit;

typedef searchHitCache<Q3RatioSearchHit> Q3RatioSearchHitCache;


class Q3GroupPeptideParser {

 public:

  Q3GroupPeptideParser(const Q3RatioSearchHitCache &searchHits, peplist* peptides, double minprob/*, Boolean heavy2light*/);
  ~Q3GroupPeptideParser();
  void setFilter(Tag* tag);
  double getRatioLogSum();
  double getRatioSquareSum();
  int getRatioNum();
  RatioStruct getRatio();

 protected:

  void parse();
  Boolean peptideListMember(const char* pep);

  peplist* peptides_;
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

  const Q3RatioSearchHitCache &searchHits_;

};











#endif
