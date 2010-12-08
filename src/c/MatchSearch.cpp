#include "MatchSearch.h"
#include "search.h"


using namespace std;

MatchSearch::MatchSearch() {
}

MatchSearch::~MatchSearch() {
}

int MatchSearch::main(int argc, char** argv){
  //TODO : Move all of match_search into this class.
  return search_main(argc, argv);
}

string MatchSearch::getName() {
  return "search-for-matches";
}

string MatchSearch::getDescription() {
  return 
  "Search a collection of spectra against a sequence "
  "database, returning a collection of peptide-spectrum "
  "matches (PSMs) scored by XCorr.";

}

string MatchSearch::getFileStem() {
  return "search";
}

