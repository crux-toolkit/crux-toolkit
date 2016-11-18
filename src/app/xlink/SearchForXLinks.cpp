/**
 * \file SearchForXLinks.cpp
 * \brief Object for running search-for-xlinks
 *****************************************************************************/
#include "SearchForXLinks.h"

#include "xlink_search.h"
#include "util/mass.h"
#include "util/Params.h"

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

  //The use-old-xlink parameter will determine
  //which codebase gets called.
  int ret;
  if (Params::GetBool("use-old-xlink")) {
    ret = xhhcSearchMain();
  } else {
    ret = xlinkSearchMain();
  }
  
  return ret;
}

void SearchForXLinks::processParams() {
  for (char c = 'A'; c <= 'Z'; c++) {
    double deltaMass = Params::GetDouble(string(1, c));
    increase_amino_acid_mass(c, deltaMass);
  }
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
    "[[html:<p>This command compares a set of spectra to cross-linked "
    "peptides derived from a protein database in FASTA format. "
    "For each spectrum, the program generates a list of candidate "
    "molecules, including linear peptides, dead-end products, self-loop "
    "products and cross-linked products, with masses that lie within a "
    "specified range of the spectrum's precursor mass. These candidate "
    "molecules are ranked using XCorr, and the XCorr scores are assigned "
    "statistical confidence estimates using an empirical curve fitting procedure."
    "</p><p>The algorithm is described in more detail in the following article:"
    "</p><blockquote>Sean McIlwain, Paul Draghicescu, Pragya Singh, David R. "
    "Goodlett and William Stafford Noble. <a href=\""
    "http://pubs.acs.org/doi/abs/10.1021/pr901163d\">"
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
    "xlink-include-inter",
    "xlink-include-intra",
    "xlink-include-inter-intra",
    "xlink-prevents-cleavage",
    "require-xlink-candidate",
    "xlink-top-n",
    "max-xlink-mods",
    "min-mass",
    "max-mass",
    "min-length",
    "max-length",
    "mod",
    "cmod",
    "nmod",
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
    "precursor-window-type-weibull",
    "min-weibull-points",
    "use-a-ions",
    "use-b-ions",
    "use-c-ions",
    "use-x-ions",
    "use-y-ions",
    "use-z-ions",
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
    "verbosity",
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > SearchForXLinks::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("search-for-xlinks.target.txt",
    "a tab-delimited text file containing the peptide-spectrum matches (PSMs). "
    "See the <a href=\"../file-formats/txt-format.html\">txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("search-for-xlinks.decoy.txt",
    "a tab-delimited text file containing the decoy PSMs. "
    "See the <a href=\"../file-formats/txt-format.html\">txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("search-for-xlinks.qvalues.txt",
    "a tab-delimited text file containing the top ranked PSMs with calculated q-values. "
    "See the <a href=\"../file-formats/txt-format.html\">txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("search-for-xlinks.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("search-for-xlinks.log.txt",
    "a log file containing a copy of all messages that were printed to stderr."));
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
