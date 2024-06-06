#include "XPressGroupPeptideParser.h"


XPressGroupPeptideParser::XPressGroupPeptideParser(const XPressRatioSearchHitCache &searchHits, const peplist* peptides, double minprob/*, Boolean heavy2light*/) : searchHits_(searchHits) {
  peptides_ = peptides;
  min_probability_ = minprob;

  ratio_square_sum_ = 0.0;
  ratio_log_sum_ = 0.0;
  h2l_ratio_square_sum_ = 0.0;
  h2l_ratio_log_sum_ = 0.0;

  ratio_num_ = 0;

  parse();
}


XPressGroupPeptideParser::~XPressGroupPeptideParser() {
}

void XPressGroupPeptideParser::parse() {

   std::vector<const XPressRatioSearchHit *> hits = searchHits_.getMatches(*peptides_);

   // now we have a list of hits on the peptides of interest
   for (std::vector<const XPressRatioSearchHit *>::iterator iter = hits.begin();iter != hits.end();iter++) {
      double nextratio;
      
      if(((*iter)->light_ > 0.0 || (*iter)->heavy_ > 0.0) && 
         (min_probability_ == 0.0 || (*iter)->probability_ >= min_probability_)) {

         if((*iter)->heavy_ == 0.0)  // inf ratio
            nextratio = pow(10.0, MAX_LOG);
         else {
            nextratio = (*iter)->light_ / (*iter)->heavy_;
            
            if(nextratio > pow(10.0, MAX_LOG))
               nextratio = pow(10.0, MAX_LOG);
            else if(nextratio < pow(10.0, -1 * MAX_LOG))
               nextratio = pow(10.0, -1 * MAX_LOG);
         }
         double h2l_nextratio;
         if (nextratio <= pow(10.0, -1 * MAX_LOG)) {
            h2l_nextratio = pow(10.0, MAX_LOG);
         }
         else if (nextratio >= pow(10.0, MAX_LOG)) {
            h2l_nextratio = pow(10.0, -1 * MAX_LOG);
         }
         else {
            h2l_nextratio = 1 / nextratio;
         }

         h2l_ratio_square_sum_ += log10(h2l_nextratio) * log10(h2l_nextratio); //nextratio * nextratio;
         h2l_ratio_log_sum_ += log10(h2l_nextratio);

         ratio_square_sum_ += log10(nextratio) * log10(nextratio); //nextratio * nextratio;
         ratio_log_sum_ += log10(nextratio);
         ratio_num_++;    
      }
   }    
}


RatioStruct XPressGroupPeptideParser::getRatio() {
  RatioStruct ratio;

  if(ratio_num_ == 0) {
      ratio.iNumPeptides = 0;
      ratio.dRatio = -9.9;
      ratio.dStdDev = -9.9;
      ratio.dh2lRatio = -9.9;
      ratio.dh2lStdDev = -9.9;
  }
  else {

    double logmean = ratio_log_sum_ / ratio_num_;
    double h2l_logmean = h2l_ratio_log_sum_ / ratio_num_;

    if(logmean >= MAX_LOG) { // infinite ratio
      ratio.dRatio = 999.;
    }
    else if(logmean <= -1 * MAX_LOG) { // zero ratio
      ratio.dRatio = 0.0;
    }
    else {
      ratio.dRatio = pow((double)10, (double)logmean);
    }
    
    if(h2l_logmean >= MAX_LOG) { // infinite ratio
      ratio.dh2lRatio = 999.;
    }
    else if(h2l_logmean <= -1 * MAX_LOG) { // zero ratio
      ratio.dh2lRatio = 0.0;
    }
    else {
      ratio.dh2lRatio = pow((double)10, (double)h2l_logmean);
    }
    


    if(ratio_num_ == 1) {
      ratio.dStdDev = 0.0;
      ratio.dh2lStdDev = 0.0;
    }
    else {
      double variance = (ratio_square_sum_ / ratio_num_) - (logmean * logmean);
      double h2l_variance = ( h2l_ratio_square_sum_ / ratio_num_) - (h2l_logmean * h2l_logmean);
      if(variance < 0.0)
         variance = 0.0;
      if(h2l_variance < 0.0)
         h2l_variance = 0.0;

      ratio.dStdDev = sqrt(variance) * ratio.dRatio; //sqrt((ratio_square_sum_ / ratio_num_) - (logmean * logmean));
      ratio.dh2lStdDev = sqrt(h2l_variance) * ratio.dh2lRatio;

    }
    ratio.iNumPeptides = ratio_num_;
  }
  //cout << ratio.dRatio << " +- " << ratio.dStdDev << " (" << ratio_num_ << ")" << endl;
  return ratio;
}

double XPressGroupPeptideParser::getRatioLogSum() {
  return ratio_log_sum_;
}

double XPressGroupPeptideParser::getRatioSquareSum() {
  return ratio_square_sum_;
}

int XPressGroupPeptideParser::getRatioNum() {
  return ratio_num_;
}


Boolean XPressGroupPeptideParser::peptideListMember(const char* pep) {
  if(peptides_ == NULL)
    return False;
  return (peptides_->find(pep) != NULL);
}
