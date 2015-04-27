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

  initialize(argc, argv);

  //The use-old-xlink parameter will determine
  //which codebase gets called.
  int ret;
  if (get_boolean_parameter("use-old-xlink")) {
    ret = xhhcSearchMain();
  } else {
    ret = xlinkSearchMain();
  }
  
  return ret;
}

/**
 * \returns the command name for SearchForXLinks
 */
string SearchForXLinks::getName() const {
  return "search-for-xlinks";
}

/**
 * \returns the description for SearchForXLinks
 */
string SearchForXLinks::getDescription() const {
  return 
    "[[nohtml:Search a collection of spectra against a sequence database, "
    "returning a collection of matches corresponding to linear and "
    "cross-linked peptides scored by XCorr.]]"
    "[[html:<p>This command searches a protein database with a set of spectra. "
    "For each spectrum, the precursor mass is computed from either the "
    "precursor singly charged mass (m+h) or the mass-to-charge (m/z) and an "
    "assumed charge. Candidates molecules are linear peptides, dead-end "
    "products, self-loop products or cross-linked products whose mass lies "
    "within a specified range of the precursor mass. These candidate peptides "
    " are ranked using XCorr. The input protein database is in FASTA format."
    "</p><p>The algorithm is described in more detail in the following article:"
    "</p><blockquote>Sean McIlwain, Paul Draghicescu, Pragya Singh, David R. "
    "Goodlett and William Stafford Noble.<a href=\""
    "http://noble.gs.washington.edu/papers/mcilwain2010detecting.html\">"
    "&quot;Detecting cross-linked peptides by searching against a database of "
    "cross-linked peptide pairs.&quot;</a> <em>Journal of Proteome Research"
    "</em>. 2010.</blockquote>]]";
}

/**
 * \returns the command arguments
 */
vector<string> SearchForXLinks::getArgs() const {
  string arr[] = {
    "ms2 file", 
    "protein fasta file", 
    "link sites", 
    "link mass"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> SearchForXLinks::getOptions() const {
  string arr[] = {
    "use-old-xlink",
    "xlink-include-linears",
    "xlink-include-deadends",
    "xlink-include-selfloops",
    "xlink-prevents-cleavage",
    "max-xlink-mods",
    "min-mass",
    "max-mass",
    "min-length",
    "max-length",
    "mod",
    "max-mods",
    "enzyme",
    "custom-enzyme",
    "digestion",
    "missed-cleavages",
    "spectrum-min-mz",
    "spectrum-max-mz",
    "spectrum-charge",
    "compute-sp",
    "precursor-window",
    "precursor-window-type",
    "precursor-window-weibull",
    "min-weibull-points",
    "max-ion-charge",
    "scan-number",
    "mz-bin-width",
    "mz-bin-offset",
    "mod-mass-format",
    "use-flanking-peaks",
    "fragment-mass",
    "isotopic-mass",
    "isotope-windows",
    "compute-p-values",
    "seed",
    "xlink-print-db",
    "spectrum-parser",
    "use-z-line",
    "top-match",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
map<string, string> SearchForXLinks::getOutputs() const {
  map<string, string> outputs;
  outputs["search-for-xlinks.params.txt"] =
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs.";
  outputs["search-for-xlinks.target.txt"] =
    "a tab-delimited text file containing the peptide-spectrum matches (PSMs). "
    "See the <a href=\"txt-format.html\">txt file format</a> for a list of the fields.";
  outputs["search-for-xlinks.decoy.txt"] =
    "a tab-delimited text file containing the decoy PSMs. "
    "See the <a href=\"txt-format.html\">txt file format</a> for a list of the fields.";
  outputs["search-for-xlinks.qvalues.txt"] =
    "a tab-delimited text file containing the top ranked PSMs with calculated q-values. "
    "See the <a href=\"txt-format.html\">txt file format</a> for a list of the fields.";
  outputs["search-for-xlinks.log.txt"] =
    "a log file containing a copy of all messages that were printed to stderr.";
  return outputs;
}

/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T SearchForXLinks::getCommand() const {
  return XLINK_SEARCH_COMMAND;
}

bool SearchForXLinks::needsOutputDirectory() const {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
