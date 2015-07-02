#include <iostream>
#include <fstream>
#include "check-match-collection.h"
#include "MatchCollection.h"
#include "MatchFileWriter.h"
#include "MatchSearch.h"
#include <vector>
#include <cstdlib>
#include <string>

// also "included" from match_collection.
void force_scored_by(MatchCollection* match_collection, SCORER_TYPE_T type);
bool calculate_delta_cn( MatchCollection* mc, COMMAND_T search_type );

using namespace std;

// helper function
// return true if a == b +/- tol
bool is_close(double a, double b, double tol){
  double diff = fabs(a - b);
  return (diff < tol);
}

// declare things to set up
MatchCollection* mc;
Match* m;
vector<Match*> match_list;
int num_matches;
// some unordered scores to use
float scores[10] = { 78.2, 50, 23.3, 109, 34.5, 50, 45.6, 50, 38, 64};
Peptide* pep;
Protein* prot;
char protseq[] = "FAKESEQ";
SpectrumZState zstate;

void match_collection_setup(){
  mc = new MatchCollection(false); // not decoy

  // set up some matches with xcorrs and a peptide to add
  prot = new Protein("prot", protseq, strlen(protseq), NULL, 0, 0, NULL);
  pep = new Peptide((unsigned char)strlen(protseq), 7.77, prot, 1);
  num_matches = 8;
  for(int i=0; i<num_matches; i++){
    Match* m = new Match(pep, NULL, zstate, false);
    m->setScore(XCORR, scores[i]);
    match_list.push_back(m);    
  }
  set_verbosity_level(CARP_ERROR);
}

void set_matches(MatchCollection* mc, vector<Match*> matches){
  for(size_t i=0; i<matches.size(); i++){
    mc->addMatch(matches[i]);
  }
  mc->forceScoredBy(XCORR);
}
void match_collection_teardown(){
  if( mc ){
    delete mc; // frees the matches
  }
  if( prot ) { delete prot; }
  if( pep ) { delete pep; }

  for(size_t i=0; i<match_list.size(); i++){
    match_list[i] = NULL;
  }
}

// test a private helper function in match_collection.cpp
void get_target_decoy_filenames(vector<string>& target_decoy_names,
                                DIR* directory,
                                SET_TYPE_T type);

START_TEST(test_list_files){
  initialize_parameters();

  // CASE: search files, one decoy file
  create_output_directory("list-files-1", true);
  system("touch list-files-1/search.target.txt");
  system("touch list-files-1/somefileroot.search.target.txt");
  system("touch list-files-1/search.target.pep.xml");
  system("touch list-files-1/search.decoy.txt");
  system("touch list-files-1/somefileroot.search.decoy.txt");
  system("touch list-files-1/search.decoy.pep.xml");
  system("touch list-files-1/search.log.txt");

  DIR* dir = opendir("list-files-1");
  // for target,decoy
  for(int type=0; type < 2; type++){
    vector<string> found_names;

    get_target_decoy_filenames( found_names, dir, (SET_TYPE_T)type);
    fail_unless( found_names.size() == 2, "Should have found 2 filenames in "
                 "list-files-1 for type %d but only found %d.", 
                 type, found_names.size() );
    // check each name
    for(size_t i=0; i<found_names.size(); i++){
      fail_unless( found_names[i].find("search") != string::npos,
                   "Filename returned (%s) did not have 'search' in it.",
                   found_names[i].c_str());
    }
  }
  closedir( dir );

  // CASE: sequest files, one decoy file
  create_output_directory("list-files-2", true);
  system("touch list-files-2/sequest.decoy.pep.xml");
  system("touch list-files-2/sequest.decoy.sqt");
  system("touch list-files-2/sequest.decoy.txt");
  system("touch list-files-2/sequest.log.txt");
  system("touch list-files-2/sequest.target.pep.xml");
  system("touch list-files-2/sequest.target.sqt");
  system("touch list-files-2/sequest.target.txt");
  system("touch list-files-2/somefileroot.sequest.decoy.txt");
  system("touch list-files-2/somefileroot.sequest.target.txt");
  
  dir = opendir("list-files-2");
  // for target,decoy
  for(int type=0; type < 2; type++){
    vector<string> found_names;
    get_target_decoy_filenames( found_names, dir, (SET_TYPE_T)type);
    fail_unless( found_names.size() == 2, "Should have found 2 filenames in "
                 "list-files-2 for type %d but only found %d.", 
                 type, found_names.size() );
    // check each name
    for(size_t i=0; i<found_names.size(); i++){
      fail_unless( found_names[i].find("sequest") != string::npos,
                   "Filename returned (%s) did not have 'sequest' in it.",
                   found_names[i].c_str());
    }
  }
  closedir( dir );

  // CASE: search files, two decoys
  create_output_directory("list-files-3", true);
  system("touch list-files-3/search.decoy-1.pep.xml");
  system("touch list-files-3/search.decoy-1.txt");
  system("touch list-files-3/search.decoy-2.pep.xml");
  system("touch list-files-3/search.decoy-2.txt");
  system("touch list-files-3/search.log.txt");
  system("touch list-files-3/search.target.pep.xml");
  system("touch list-files-3/search.target.txt");
  system("touch list-files-3/somefileroot.search.decoy-1.txt");
  system("touch list-files-3/somefileroot.search.decoy-2.txt");
  system("touch list-files-3/somefileroot.search.target.txt");

  dir = opendir("list-files-3");
  // for target,decoy
  for(int type=0; type < 3; type++){
    vector<string> found_names;
    get_target_decoy_filenames( found_names, dir, (SET_TYPE_T)type);
    fail_unless( found_names.size() == 2, "Should have found 2 filenames in "
                 "list-files-3 for type %d but only found %d.", 
                 type, found_names.size() );
    // check each name
    for(size_t i=0; i<found_names.size(); i++){
      fail_unless( found_names[i].find("search") != string::npos,
                   "Filename returned (%s) did not have 'search' in it.",
                   found_names[i].c_str());
    }
  }
  closedir( dir );
}
END_TEST

START_TEST(test_create){
  // TODO add more create tests
  fail_unless( mc != NULL, "New empty match collection should not equal NULL.");
}
END_TEST

START_TEST(test_set){
  set_matches(mc, match_list);
  fail_unless( mc->getMatchTotal() == num_matches,
               "Added %d matches but collection only has %d.", 
               num_matches, mc->getMatchTotal());

  // check scores of all matches
  MatchIterator* mi = new MatchIterator(mc, XCORR, false); // don't sort
  int idx = 0;
  while( mi->hasNext() ){
    Match* cur_match = mi->next();
    double xcorr = cur_match->getScore(XCORR);
    fail_unless( xcorr == scores[idx], 
                 "Xcorr of match %d is %.1f but should be %.1f.",
                 idx, xcorr, scores[idx]);
    idx++;
  }

  // also check that they sort correctly
  delete mi;
  mi = new MatchIterator(mc, XCORR, true); // do sort
  idx = 0;
  double last_score = 10000;
  while( mi->hasNext() ){
    Match* cur_match = mi->next();
    double xcorr = cur_match->getScore(XCORR);
    fail_unless( xcorr <= last_score, 
                 "Xcorr of match %d is %.1f and should be less than %.1f.",
                 idx, xcorr, last_score);
    idx++;
    last_score = xcorr;
  }
  delete mi;
}
END_TEST
// once these tests pass, we can use set_matches() and match_iterator in other tests

START_TEST(test_print_rank){
  initialize_parameters();
  set_matches(mc, match_list);
  mc->populateMatchRank(XCORR);

  // other items for printing to tab file
  int z = 1;
  vector<int> possible_z(z);
  Spectrum* s = new Spectrum(7, 7, 7.77, possible_z, (char *)"fakename");
  const char* filename = "test-rank.txt";

  /* expected values
   rank 1 xcor 109
   rank 2 xcor 78.2
   rank 3 xcor 50
   rank 3 xcor 50
   rank 3 xcor 50
   rank 4 xcor 45.6
   rank 5 xcor 34.5
   rank 6 xcor 23.3
  */

  // try printing the top 3, should get 5 lines
  unlink(filename);
  MatchFileWriter* fout = new MatchFileWriter(filename);
  MatchSearch application;
  fout->addColumnNames(&application, false);
  fout->writeHeader();
  mc->printTabDelimited(fout, 3, s, XCORR);
  delete fout;
  sleep(2);  // wait to avoid NFS latency issues

  ifstream fin(filename, ifstream::in);
  string line;
  int count = -1; // it will count once after getting the eof
  while( ! fin.eof() ) {
    getline(fin, line);
    count++;
  }

  fail_unless( count == 6, 
               "For top-match=3, there should have been 6 lines printed "
               "but there were %d.", count);
  fin.close();

  // try printing top 5, should still get 5 lines
  unlink(filename);
  fout = new MatchFileWriter(filename);
  fout->addColumnNames(&application, false);
  fout->writeHeader();
  mc->printTabDelimited(fout, 5, s, XCORR);
  delete fout;

  fin.open(filename, ifstream::in);
  count = -1; // it will count once after getting the eof
  while( ! fin.eof() ) {
    getline(fin, line);
    count++;
  }

  fail_unless( count == 6, 
               "For top-match=3, there should have been 6 lines printed "
               "but there were %d.", count);
  fin.close();

  // try printing top 6, should get 6
  unlink(filename);
  fout = new MatchFileWriter(filename);
  fout->addColumnNames(&application, false);
  fout->writeHeader();
  mc->printTabDelimited(fout, 6, s, XCORR);
  delete fout;

  fin.open(filename, ifstream::in);
  count = -1; // it will count once after getting the eof
  while( ! fin.eof() ) {
    getline(fin, line);
    count++;
  }

  fail_unless( count == 7, 
               "For top-match=3, there should have been 7 lines printed "
               "but there were %d.", count);
  fin.close();

}
END_TEST

START_TEST(test_delta_cn){
  set_matches(mc, match_list);
  mc->calculateDeltaCn(SEQUEST_COMMAND);

  // check delta_cn scores
  MatchIterator* mi = new MatchIterator(mc, XCORR, false); // don't sort
  Match* cur_match = mi->next();
  double first_xcorr = cur_match->getScore(XCORR);

  // for sequest, delta cn of first should be 0
  double dcn = cur_match->getDeltaCn();
  fail_unless(dcn == 0, "First deltaCn should be 0 for SEQUEST but is %.1f", 
              dcn);

  // check remaining delta cns and save them for next test
  vector<double> dcns;
  while( mi->hasNext() ){
    cur_match = mi->next();
    dcn = cur_match->getDeltaCn();
    if( dcns.empty() || dcn != dcns.back() ){ 
      dcns.push_back(dcn); 
    }
    double predicted_dcn = (first_xcorr - cur_match->getScore(XCORR))/first_xcorr;
    fail_unless( is_close(dcn, predicted_dcn, 0.001), 
                 "DeltaCn is %.10f and should be %.10f", dcn, predicted_dcn);
  }

  // now test with search command
  mc->calculateDeltaCn(SEARCH_COMMAND);

  // check delta_cn scores
  delete mi;
  mi = new MatchIterator(mc, XCORR, false); // don't sort
  cur_match = mi->next();

  // for non sequest, delta cn of first should not be 0
  dcn = cur_match->getDeltaCn();
  fail_unless(dcn != 0, "First deltaCn should not be 0 for SEARCH but is %.1f", 
              dcn);

  // check deltaCn compared to SEQUEST values
  double last_xcorr = cur_match->getScore(XCORR);
  vector<double>::iterator previous_dcn = dcns.begin();

  while( mi->hasNext() ){

    fail_unless( is_close(dcn, *previous_dcn, 0.001), 
                 "DeltaCn is %.10f and should be %.10f", dcn, *previous_dcn);
    
    last_xcorr = cur_match->getScore(XCORR);
    cur_match = mi->next();
    dcn = cur_match->getDeltaCn();

    if( last_xcorr != cur_match->getScore(XCORR) &&
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
  tcase_add_test(tc_core, test_list_files);
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_set);
  tcase_add_test(tc_core, test_print_rank);
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
