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
    "</em>. 2010.</blockquote><p><strong>Modifications:</strong> Currently, "
    "<code>crux search-for-xlinks</code> supports static modifications (a "
    "change of mass applied to a given amino acid in every peptide in which it "
    "occurs). By default, a static modification of +57 Da to cysteine (C) is "
    "applied. Variable modifications (allowing peptides to be generated with "
    "and without a mass change to a given amino acid), are supported when <code>"
    "use-old-xlink=F</code>. Static and variable modifications can be "
    "specified in the parameter file, as described below.</p>]]";
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
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
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
    "See the txt file format for a list of the fields.";
  outputs["search-for-xlinks.decoy.txt"] =
    "a tab-delimited text file containing the decoy PSMs. See the txt file "
    "format for a list of the fields.";
  outputs["search-for-xlinks.qvalues.txt"] =
    "a tab-delimited text file containing the top ranked PSMs with calculated "
    "q-values. See the txt file format for a list of the fields.";
  outputs["search-for-xlinks.log.txt"] =
    "a log file containing a copy of all messages that were printed to "
    "stderr.";
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
