/**
 * \file SearchProgress.h
 * AUTHOR: Barbara Frewen
 *  CREATE DATE: October 5, 2009
 * \brief A class to keep track of how many searches have been
 * attempted and how many have been successfull.  Reports progress at
 * appropriate intervals
 */

#include "carp.h"
#include "parameter.h"


using namespace std;

/**
 * \class SearchProgress
 *
 * A class to keep track of how many searches have been performed and
 * output the information at appropriate intervals.  Instantiate a
 * class at the beginning of the search and call increment() after
 * each search.  Frequency of output messages is determined by the
 * parameter "print-search-progress".
 */
class SearchProgress{

 public:
  /**
   * Default constructor.
   */
  SearchProgress() :
    searches_attempted_(0), 
    searches_with_matches_(0), 
    progress_increment_(get_int_parameter("print-search-progress")),
    total_searches_(1)
    { 
      report_format_ = "Searching spectrum %i (%i+), search %i";
    };

  /**
   * Constructor that also sets the total number of searches expected
   * to be performed.
   */
  SearchProgress(int total_searches) :
    searches_attempted_(0), 
    searches_with_matches_(0), 
    progress_increment_(get_int_parameter("print-search-progress")),
    total_searches_(total_searches)
    { 
      report_format_ = "Searching spectrum %i (%i+), search %i"
        " of %i, %.0f%% complete";
    };

    void report(int scan_num, int charge){
      if( ((searches_attempted_ + 1) % progress_increment_) == 0 ){
        carp(CARP_INFO, report_format_, 
             scan_num, charge, searches_attempted_ + 1, total_searches_,
             (float)100 * (searches_attempted_ + 1) / (float)total_searches_ );
      }
    };

    void increment(bool generated_matches){
      searches_attempted_++;
      if( generated_matches ){
        searches_with_matches_++;
      }
    };

    int getNumSearchesWithMatches(){ return searches_with_matches_; };

 private:
  int searches_attempted_;    ///< number of spec/charge's searched 
  int searches_with_matches_; ///< number of spec with results in .txt file
  int progress_increment_;    ///< how often to print progress
  int total_searches_;        ///< expected number of searches to perform
  const char* report_format_; ///< format string to print
};







/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
