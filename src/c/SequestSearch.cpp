/**
 * \file SequestSearch.cpp
 * \brief Object for running sequest-search
 *****************************************************************************/
#include "SequestSearch.h"
#include "sequest-search.h"

using namespace std;

/**
 * \returns a blank SequestSearch object
 */
SequestSearch::SequestSearch() {

}

/**
 * Destructor
 */
SequestSearch::~SequestSearch() {
}

/**
 * main method for SequestSearch
  */
int SequestSearch::main(int argc,   ///< number of cmd line tokens
                        char** argv)///< array of cmd line tokens
{

  return sequest_search_main(argc, argv);

}// end main

/**
 * \returns the command name for SequestSearch
 */
string SequestSearch::getName() {
  return "sequest-search";
}

/**
 * \returns the description for SequestSearch
 */
string SequestSearch::getDescription() {
  return 
  "Similar to search-for-matches but use Sp as a "
  "preliminary score followed by XCorr.";

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
