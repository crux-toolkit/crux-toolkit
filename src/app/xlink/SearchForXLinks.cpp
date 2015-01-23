/**
 * \file SearchForXLinks.cpp
 * \brief Object for running search-for-xlinks
 *****************************************************************************/
#include "SearchForXLinks.h"

#include "xlink_search.h"

using namespace std;

/**
 * \returns a blank SearchForXLinks object
 */
SearchForXLinks::SearchForXLinks() {

}

/**
 * Destructor
 */
SearchForXLinks::~SearchForXLinks() {
}

/**
 * main method for SearchForXLinks
 */
int SearchForXLinks::main(int argc, char** argv) {

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);

  /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {
    "ms2 file", 
    "protein-database", 
    "link sites", 
    "link mass"
  };

  int num_arguments = sizeof(argument_list) / sizeof(char*);


  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);
  

  //The use-old-xlink parameter will determine
  //which codebase gets called.
  int ret;
  if (get_boolean_parameter("use-old-xlink")) {
    ret = xhhcSearchMain();
  } else {
    ret = xlinkSearchMain();
  }
  
  free_parameters();
  return ret;
}

/**
 * \returns the command name for SearchForXLinks
 */
string SearchForXLinks::getName() {
  return "search-for-xlinks";
}

/**
 * \returns the description for SearchForXLinks
 */
string SearchForXLinks::getDescription() {
  return 
    "Search a collection of spectra against a sequence "
    "database, returning a collection of matches "
    "corresponding to linear and cross-linked peptides "
    "scored by XCorr.";
}

/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T SearchForXLinks::getCommand() {
  return XLINK_SEARCH_COMMAND;
}

bool SearchForXLinks::needsOutputDirectory() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
