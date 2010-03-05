#include "check-match-collection.h"
#include "match_collection.h"
#include <vector>
#include <cstdlib>

// also "included" from match_collection.
void force_scored_by(MATCH_COLLECTION_T* match_collection, SCORER_TYPE_T type);
BOOLEAN_T calculate_delta_cn( MATCH_COLLECTION_T* mc, COMMAND_T search_type );

using namespace std;

// helper function
// return true if a == b +/- tol
bool is_close(double a, double b, double tol){
  double diff = fabs(a - b);
  return (diff < tol);
}

// declare things to set up
MATCH_COLLECTION_T* mc;
MATCH_T* m;
vector<MATCH_T*> match_list;
int num_matches;
// some unordered scores to use
float scores[10] = { 78.2, 50, 23.3, 109, 34.5, 50, 45.6, 50, 38, 64};

void match_collection_setup(){
  mc = new_empty_match_collection(FALSE); // not decoy

  // set up some matches with xcorrs to add
  num_matches = 8;
  for(int i=0; i<num_matches; i++){
    MATCH_T* m = new_match();
    set_match_score(m, XCORR, scores[i]);
    match_list.push_back(m);    
  }
}

void set_matches(MATCH_COLLECTION_T* mc, vector<MATCH_T*> matches){
  for(size_t i=0; i<matches.size(); i++){
    add_match_to_match_collection(mc, matches[i]);
  }
  force_scored_by(mc, XCORR);
}
void match_collection_teardown(){
  if( mc ){
    free(mc);
  }
  for(size_t i=0; i<match_list.size(); i++){
    free_match(match_list[i]);
    match_list[i] = NULL;
  }
}

START_TEST(test_create){
  // TODO add more create tests
  fail_unless( mc != NULL, "New empty match collection should not equal NULL.");
}
END_TEST

START_TEST(test_set){
  set_matches(mc, match_list);
  fail_unless( get_match_collection_match_total(mc) == num_matches,
               "Added %d matches but collection only has %d.", 
               num_matches, get_match_collection_match_total(mc));

  // check scores of all matches
  MATCH_ITERATOR_T* mi = new_match_iterator(mc, XCORR, FALSE); // don't sort
  int idx = 0;
  while( match_iterator_has_next(mi) ){
    MATCH_T* cur_match = match_iterator_next(mi);
    double xcorr = get_match_score(cur_match, XCORR);
    fail_unless( xcorr == scores[idx], 
                 "Xcorr of match %d is %.1f but should be %.1f.",
                 idx, xcorr, scores[idx]);
    idx++;
  }

  // also check that they sort correctly
  free_match_iterator(mi);
  mi = new_match_iterator(mc, XCORR, TRUE); // do sort
  idx = 0;
  double last_score = 10000;
  while( match_iterator_has_next(mi) ){
    MATCH_T* cur_match = match_iterator_next(mi);
    double xcorr = get_match_score(cur_match, XCORR);
    fail_unless( xcorr <= last_score, 
                 "Xcorr of match %d is %.1f and should be less than %.1f.",
                 idx, xcorr, last_score);
    idx++;
    last_score = xcorr;
  }
  free_match_iterator(mi);
}
END_TEST
// once these tests pass, we can use set_matches() and match_iterator in other tests

START_TEST(test_delta_cn){
  set_matches(mc, match_list);
  calculate_delta_cn(mc, SEQUEST_COMMAND);

  // check delta_cn scores
  MATCH_ITERATOR_T* mi = new_match_iterator(mc, XCORR, FALSE); // don't sort
  MATCH_T* cur_match = match_iterator_next(mi);
  double first_xcorr = get_match_score(cur_match, XCORR);

  // for sequest, delta cn of first should be 0
  double dcn = get_match_delta_cn(cur_match);
  fail_unless(dcn == 0, "First deltaCn should be 0 for SEQUEST but is %.1f", 
              dcn);

  // check remaining delta cns and save them for next test
  vector<double> dcns;
  while( match_iterator_has_next(mi) ){
    cur_match = match_iterator_next(mi);
    dcn = get_match_delta_cn(cur_match);
    if( dcns.empty() || dcn != dcns.back() ){ 
      dcns.push_back(dcn); 
    }
    double predicted_dcn = (first_xcorr - get_match_score(cur_match, XCORR))/first_xcorr;
    fail_unless( is_close(dcn, predicted_dcn, 0.001), 
                 "DeltaCn is %.10f and should be %.10f", dcn, predicted_dcn);
  }

  // now test with search command
  calculate_delta_cn(mc, SEARCH_COMMAND);

  // check delta_cn scores
  free_match_iterator(mi);
  mi = new_match_iterator(mc, XCORR, FALSE); // don't sort
  cur_match = match_iterator_next(mi);

  // for non sequest, delta cn of first should not be 0
  dcn = get_match_delta_cn(cur_match);
  fail_unless(dcn != 0, "First deltaCn should not be 0 for SEARCH but is %.1f", 
              dcn);

  // check deltaCn compared to SEQUEST values
  double last_xcorr = get_match_score(cur_match, XCORR);
  vector<double>::iterator previous_dcn = dcns.begin();

  while( match_iterator_has_next(mi) ){

    fail_unless( is_close(dcn, *previous_dcn, 0.001), 
                 "DeltaCn is %.10f and should be %.10f", dcn, *previous_dcn);
    
    last_xcorr = get_match_score(cur_match, XCORR);
    cur_match = match_iterator_next(mi);
    dcn = get_match_delta_cn(cur_match);

    if( last_xcorr != get_match_score(cur_match, XCORR) &&
        previous_dcn + 1 < dcns.end() ){
      previous_dcn++;
    }
  }
}
END_TEST
/* Boundry conditions test suite */
START_TEST(test_null){
  fail_unless( 1 == 1 );
}
END_TEST

Suite* match_collection_suite(){
  Suite* s = suite_create("match_collection");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_set);
  tcase_add_test(tc_core, test_delta_cn);
  //  tcase_add_test(tc_core, test_...);

  tcase_add_checked_fixture(tc_core, 
                            match_collection_setup, 
                            match_collection_teardown);
  suite_add_tcase(s, tc_core);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  tcase_add_test(tc_limits, test_null);
  tcase_add_checked_fixture(tc_limits, 
                            match_collection_setup, 
                            match_collection_teardown);
  suite_add_tcase(s, tc_limits);

  return s;
}













