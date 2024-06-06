#ifndef scan_reader_class_h
#define scan_reader_class_h

#include "deep_class.h"
#include "Parsers/mzParser/mzParser.h"

#include <math.h>
#include <numeric>

class scan_reader {
public:
  bool getTIC(peptide_lists& my_peptide_lists, my_parameters& my_params);
  bool mzml(peptide_lists& my_peptide_lists, my_parameters& my_params);
   
private:

  //MH: Binary search function. Could and should be overloaded for PPM instead of MZ tolerance
  int findPeakMZ(mzParser::BasicSpectrum& spec, double mz, double tol);

  //MH: quick helper function for sorting from lowest to highes RT
  static bool compareRTime(const dsPSM& a, const dsPSM& b) { return a.xml_rtime < b.xml_rtime; }
  //static bool compareSeqScan(const dsPeptide& a, const dsPeptide& b) { 
  //  int i=a.pep_seq.compare(b.pep_seq);
  //  if(i==0) return (a.spec_sn<b.spec_sn);
  //  else return (i<0);
  //}
  

};

#endif
