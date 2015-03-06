/**
 * \file ComputeQValues.cpp
 * \brief Object for calling compute-q-values
 *****************************************************************************/
#include "ComputeQValues.h"
#include "io/OutputFiles.h"
#include "analyze_psms.h"

using namespace std;


/**
 * \returns a blank ComputeQValues object
 */
ComputeQValues::ComputeQValues() {

}

/**
 * Destructor
 */
ComputeQValues::~ComputeQValues() {
}

/**
 * main method for ComputeQValues
 */
int ComputeQValues::main(int argc, char** argv) {
  //TODO : Figure out how to do this.
  analyze_matches_main(argc, argv);
  return 0;
}

/**
 * \returns the command name for ComputeQValues
 */
string ComputeQValues::getName() const {
  return "assign-confidence";
}

/**
 * \returns the description for ComputeQValues
 */
string ComputeQValues::getDescription() const {
  return
    "[[nohtml:Assign two types of statistical confidence measures (q-values "
    "and posterior error probabilities) to each PSM in a given set.]]"
    "[[html:<p>Given a collection of scored peptide-spectrum matches (PSMs), "
    "estimate two statistical confidence measures for each: a q-value and a "
    "posterior error probability (PEP).</p><h3>q-value</h3><p>The q-value is "
    "analogous to a p-value but incorporates false discovery rate multiple "
    "testing correction. The q-value associated with a score threshold T is "
    "defined as the minimal false discovery rate at which a score of T is "
    "deemed significant. In this setting, the q-value accounts for the fact "
    "that we are analyzing a large collection of PSMs.</p><p>To estimate "
    "q-values, <code>assign-confidence</code> searches the input directory for "
    "a corresponding set of decoy PSMs. The false discovery rate associated "
    "with a given score is estimated as the number of decoy scores above the "
    "threshold divided by the number of target scores above the threshold, "
    "multiplied by the ratio of the total number of targets to total number of "
    "decoys. This methodology is described in the following article:</p>"
    "<blockquote>Lukas K&auml;ll, John D. Storey, Michael J. MacCoss and "
    "William Stafford Noble. <a href=\"http://noble.gs.washington.edu/papers/"
    "kall2008assigning.html\">&quot;Assigning significance to peptides "
    "identified by tandem mass spectrometry using decoy databases&quot;</a>. "
    "<em>Journal of Proteome Research</em>. 7(1):29-34, 2008.</blockquote><p>"
    "Note that assign-confidence does not (yet) estimate the percentage of "
    "incorrect targets, as described in the above article. Hence, the method "
    "implemented here as &quot;decoy q-values&quot; is analogous to the &quot;"
    "Simple FDR&quot; procedure shown in Figure 4A of the above article.</p><p>"
    "In each case, the estimated FDRs are converted to q-values by ranking the "
    "PSMs by score and then taking, for each PSM, the minimum of the current "
    "FDR and all of the FDRs below it in the ranked list.</p><h3>Posterior "
    "error probability</h3><p>Unlike the q-value, which is calculated with "
    "respect to the collection of PSMs with scores above a specified "
    "threshold, the PEP (also known in the literature as the &quot;local "
    "FDR&quot;) is calculated with respect to a single score. The PEP is the "
    "probability that a particular PSM is incorrect. Crux's PEPs are estimated "
    "using the methodology described in this article:</p><blockquote>Lukas "
    "K&auml;ll, John Storey and William Stafford Noble. <a href=\""
    "http://noble.gs.washington.edu/papers/kall2008nonparametric.html\">"
    "&quot;Non-parametric estimation of posterior error probabilities "
    "associated with peptides identified by tandem mass spectrometry&quot;"
    "</a>. <em>Bioinformatics (Proceedings of the ECCB)</em>. 24(16):i42-i48, "
    "2008.</blockquote><p>A primer on multiple testing correction can be found "
    "here</p><blockquote>William Stafford Noble. <a href=\""
    "http://noble.gs.washington.edu/papers/noble2009how.html\">&quot;How does "
    "multiple testing correction work?&quot;</a> <em>Nature Biotechnology</em>. "
    "27(12):1135-1137, 2009.</blockquote><p>A discussion of q-values versus "
    "posterior error probabilities is provided in this article:</p><blockquote>"
    "Lukas K&auml;ll, John D. Storey, Michael J. MacCoss and William Stafford "
    "Noble. <a href=\"http://noble.gs.washington.edu/papers/"
    "kall2008posterior.html\">&quot;Posterior error probabilities and false "
    "discovery rates: two sides of the same coin&quot;</a>. <em>Journal of "
    "Proteome Research</em>. 7(1):40-44, 2008.</blockquote>]]";
}

/**
 * \returns the command arguments
 */
vector<string> ComputeQValues::getArgs() const {
  string arr[] = {
    "target input"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> ComputeQValues::getOptions() const {
  string arr[] = {
    "estimation-method",
    "decoy-prefix",
    "score",
    "smaller-is-better",
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
    "list-of-files",
    "fileroot"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
map<string, string> ComputeQValues::getOutputs() const {
  map<string, string> outputs;
  outputs["assign-confidence.target.txt"] =
    "a tab-delimited text file containing the PSMs. See txt file format for a "
    "list of the fields. The file will contain two additional columns, named "
    "\"<column name> q-value\" and \"<column name> PEP\" where "
    "\"<column name>\" is provided on the command line.";
  outputs["assign-confidence.log.txt"] =
    "a log file containing a copy of all messages that were printed to "
    "stderr.";
  outputs["assign-confidence.params.txt"] =
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs.";
  return outputs;
}

/**
 * \returns the filestem for ComputeQValues
 */
string ComputeQValues::getFileStem() const {
  return "assign-confidence";
}

COMMAND_T ComputeQValues::getCommand() const {
  return QVALUE_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not.
 */
bool ComputeQValues::needsOutputDirectory() const {
  return true;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
