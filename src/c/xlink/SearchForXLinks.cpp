#include "SearchForXLinks.h"

#include "xlink_search.h"

using namespace std;

SearchForXLinks::SearchForXLinks() {

}

SearchForXLinks::~SearchForXLinks() {
}


int SearchForXLinks::main(int argc, char** argv) {
  return xlink_search_main(argc, argv);
}

string SearchForXLinks::getName() {
  return "search-for-xlinks";
}

string SearchForXLinks::getDescription() {
  return 
    "Search a collection of spectra against a sequence "
    "database returning a collection of matches "
    "corresponding to linear and cross-lnked peptides "
    "scored by XCorr.";
}
