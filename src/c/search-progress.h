/**
 * \file search-progress.h
 */

/*
 AUTHOR: Barbara Frewen
 CREATE DATE: October 5, 2009
 DESCRIPTION: A class to keep track of how many searches have been
 attempted and how many have been successfull.  Reports progress at
 appropriate intervals     
 REVISION: $Revision: 1.0 $

 */
#ifdef __cplusplus
extern "C" {
#endif

#include "carp.h"
#include "parameter.h"

#ifdef __cplusplus
}
#endif

using namespace std;

class SearchProgress{

 public:
  SearchProgress() :
    searches_attempted_(0), 
    searches_with_matches_(0), 
    progress_increment_(get_int_parameter("print-search-progress"))
      { };

    void report(int scan_num, int charge){
      if( ((searches_attempted_ + 1) % progress_increment_) == 0 ){
        carp(CARP_INFO, 
             "Searching spectrum number %i, charge %i, search number %i",
             scan_num, charge, searches_attempted_ + 1 );
      }
    };

    void increment(BOOLEAN_T generated_matches){
      searches_attempted_++;
      if( generated_matches ){
        searches_with_matches_++;
      }
    };

    int getNumSearchesWithMatches(){ return searches_with_matches_; };

 private:
  int searches_attempted_; ///< number of spec/charge's searched 
  int searches_with_matches_; ///< number of spec with results in .csm file
  int progress_increment_;  ///< how often to print progress

};







/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
