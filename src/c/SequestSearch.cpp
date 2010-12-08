#include "SequestSearch.h"

#include "sequest-search.h"


using namespace std;

SequestSearch::SequestSearch() {

}

SequestSearch::~SequestSearch() {
}

int SequestSearch::main(int argc,   ///< number of cmd line tokens
                        char** argv)///< array of cmd line tokens
{

  return sequest_search_main(argc, argv);

}// end main

string SequestSearch::getName() {
  return "sequest-search";
}

string SequestSearch::getDescription() {
  return 
  "Similar to search-for-matches but use Sp as a "
  "preliminary score followed by XCorr.";

}

string SequestSearch::getFileStem() {
  return "sequest";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
