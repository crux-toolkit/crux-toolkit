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
  set_modspec();
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
    "</em>. 2010.</blockquote>"
    "<p>In search-for-xlinks, properties of the cross-linker are specified "
    "using the two required command line arguments, <link sites> and <link "
    "mass>. In addition, mass shifts associated with mono-links can be "
    "specified using the --mod option. Below are suggested parameter settings "
    "for some commonly used cross-linkers:</p>"
    "<table border=\"1\">"
    "  <tr>"
    "    <td><b>Linker</b></td>"
    "    <td><b>Link Mass</b></td>"
    "    <td><b>Link Sites</b></td>"
    "    <td><b>Mono Link</b></td>"
    "  </tr>"
    "  <tr>"
    "    <td>EDC</td>"
    "    <td>-18.0156</td>"
    "    <td>K:D,E,cterm</td>"
    "    <td>&nbsp;</td>"
    "  </tr>"
    "  <tr>"
    "    <td>BS2</td></td>"
    "    <td>96.0211296</td>"
    "    <td>K,nterm:K,nterm</td>"
    "    <td>--mono-link 1K+114.0316942:T:T,1K+113.0476524:T:T</td>"
    "  </tr>"
    "  <tr>"
    "    <td>BS3</td>"
    "    <td>138.0680742</td>"
    "    <td>K,nterm:K,nterm</td>"
    "    <td>--mono-link 1K+156.0786:T:T,1K+155.0946278:T:T</td>"
    "  </tr>"
    /* Need to find out the right parameters for SDA.
    "  <tr>"
    "    <td>SDA</td>"
    "    <td>82.0412979</td>"
    "    <td>K,nterm:X,nterm</td>"
    "    <td>--mod ????</td>"
    "  </tr>"
    */
    "  <tr>"
    "    <td>DSS</td>"
    "    <td>138.0680796</td>"
    "    <td>K,nterm:K,nterm</td>"
    "    <td>--mono-link 1K+156.0786:T:T,1K+155.0946278</td>"
    "  </tr>"
    "  <tr>"
    "    <td>AMAS</td>"
    "    <td>137.011</td>"
    "    <td>K,nterm:C</td>"
    "    <td>--mono-link 1KC+155.02156:T:T</td>"
    "  </tr>"
    "  <tr>"
    "    <td>GMBS</td>"
    "    <td>165.0422</td>"
    "    <td>K,nterm:C</td>"
    "    <td>--mono-link 1KC+183.05276:T:T</td>"
    "  </tr>"
    "  <tr>"
    "    <td>formaldehyde</td>"
    "    <td>9.98435</td>"
    "    <td>K,W,nterm:H,N,Y,K,W,R,nterm</td>"
    "    <td>--mono-link 1KW+12.0:T:T,1KW+30.010565:T:T</td>"
    "  </tr>"
    "</table>"

    "<ul>"
    "<li>"
    "Note that, unlike Tide, search-for-xlinks does not have a "
    "\"decoy-format\" option; instead, shuffled decoy peptides "
    "are always created.  In particular, the decoy database contains "
    "three copies of each target cross-linked species: one in which "
    "both peptides are shuffled, one in which only the first peptide "
    "is shuffled, and one in which only the second peptide is shuffled. "
    "In the tab-delimited output, these different types are indicated "
    "in the \"target/decoy\" column as \"target-target,\" "
    "\"target-decoy,\" \"decoy-target\" or \"decoy-decoy.\"</li>"

    "<li>"
    "In addition to the primary XCorr score, search-for-xlinks reports "
    "separate scores for the two cross-linked peptides. The way this "
    "calculation is done depends on whether the \"top-n\" parameter is "
    "set to 0 or not.  In the top-n=0 case, the XCorr scores of the "
    "two participating peptides are computed exactly.  When top-n is "
    "non-zero, then each peptide's score is calculated by using as "
    "the mass shift on the link site the remainder of the precursor "
    "mass. The crosslink peptide score is then calculated using the "
    "sum of the two peptide scores. In this approach, the true mass "
    "shift assigned to each peptide within a crosslinked pair can be "
    "different from that peptide's calculated mass shift. This mass "
    "shift affects the way the ions masses are calculated around the "
    "crosslink sites and will affect the final XCorr score. In "
    "particular, the smaller the precursor tolerance window, the closer "
    "the \"full, slower\" XCorr value calculated directly from the "
    "crosslinked peptides and the \"fast, inaccurate\" XCorr value "
    "calculated by summing the individual peptide scores using the "
    "remainder precursor mass shift will be. This is because the pairs "
    "of peptides chosen will have a \"precursor remainder\" mass shift "
    "closer to the true mass shift calculated directly from the "
    "crosslinked peptide.</li>"
    "</ul>"
    "]]";
}

/**
 * \returns the command arguments
 */
vector<string> SearchForXLinks::getArgs() const {
  string arr[] = {
    "ms2 file+", 
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
    "mods-spec",
    "cmod",
    "nmod",
    "max-mods",
    "mono-link",
    "mod-precision",
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
    "file-column",
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
    "concat",
    "xlink-print-db",
    "spectrum-parser",
    "use-z-line",
    "top-match",
    "print-search-progress",
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
