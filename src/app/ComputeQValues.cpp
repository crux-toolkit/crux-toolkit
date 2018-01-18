/**
 * \file ComputeQValues.cpp
 * \brief Object for calling compute-q-values
 *****************************************************************************/
#include "ComputeQValues.h"
#include "io/OutputFiles.h"
#include "util/Params.h"
#include "PosteriorEstimator.h"


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
  vector<string> input_files = Params::GetStrings("target input");  

  if (Params::GetBool("list-of-files")) {
    get_files_from_list(input_files[0], input_files);
  }
  return main(input_files);
}

int ComputeQValues::main(const vector<string>& input_files) {
  // Prepare the output files.
 
  return 0;
}

/**
 * Compute posterior error probabilities (PEP) from the given target
 * and decoy scores.
 * \returns A newly allocated array of PEP for the target scores
 * sorted.
 */
vector<double> ComputeQValues::compute_PEP(
  const vector<FLOAT_T>& target_scores, ///< scores for target matches
  const vector<FLOAT_T>& decoy_scores,  ///< scores for decoy matches
  bool ascending ///< are the scores ascending or descending
) {
  if (target_scores.empty() || decoy_scores.empty()) {
    carp(CARP_FATAL, "Cannot compute PEP without target or decoy scores.");
  }

  vector< pair<double, bool> > score_labels =
    getScoreVector(target_scores, decoy_scores, ascending);
  double pi0 = estimatePi0(score_labels);

  // estimate PEPs
  vector<double> PEP_vector;

  PosteriorEstimator::estimatePEP(score_labels,
    true, // usePi0
    pi0, PEP_vector, 
    true);  // include decoy PEPs

  // now score_labels and PEPs are similarly sorted

  // pull out the PEPs in the order that the scores were given
  vector<double> pep(PEP_vector.size(), 0);
  for (int target_idx = 0; target_idx < target_scores.size(); target_idx++) {
    // the score to return next    
    double curr_target_score = target_scores[target_idx];

    // find its position in score_labels
    vector< pair<double, bool> >::iterator found_score_pos;
    if (ascending) {
      found_score_pos 
        = lower_bound(score_labels.begin(), score_labels.end(), 
                      make_pair(curr_target_score, true));
    } else {
      found_score_pos 
        = lower_bound(score_labels.begin(), score_labels.end(), 
                    make_pair(curr_target_score, true),
                    greater<pair<double, bool> >()); 
    }

    size_t found_index = distance(score_labels.begin(), found_score_pos);

    // pull out the PEP at the same position in PEP_vector
    pep[target_idx] = PEP_vector[found_index];
  }
  return pep;
}

// put all of the scores in a single vector of pairs: score, is_target
vector< pair<double, bool> > ComputeQValues::getScoreVector(
  const vector<FLOAT_T>& targetScores,
  const vector<FLOAT_T>& decoyScores,
  bool ascending
) {
  vector< pair<double, bool> > scores;
  scores.reserve(targetScores.size() + decoyScores.size());
  for (vector<FLOAT_T>::const_iterator i = targetScores.begin(); i != targetScores.end(); i++) {
    scores.push_back(make_pair(*i, true));
  }
  for (vector<FLOAT_T>::const_iterator i = decoyScores.begin(); i != decoyScores.end(); i++) {
    scores.push_back(make_pair(*i, false));
  }

  // sort them 
  if (ascending) {
    sort(scores.begin(), scores.end());
    PosteriorEstimator::setReversed(true);
  } else {
    sort(scores.begin(), scores.end(), greater< pair<double, bool> > ());
  }

  return scores;
}

double ComputeQValues::estimatePi0(
  const vector< pair<double, bool> >& scoreVector
) {
  vector<double> pVals;
  PosteriorEstimator::getPValues(scoreVector, pVals);
  return PosteriorEstimator::estimatePi0(pVals);
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
    "[[html:<p>Given target and decoy scores, estimate a q-value for each "
    "target score. The q-value is analogous to a p-value but incorporates "
    "false discovery rate multiple testing correction. The q-value associated "
    "with a score threshold T is defined as the minimal false discovery rate "
    "(FDR) at which a score of T is deemed significant. In this setting, the "
    "q-value accounts for the fact that we are analyzing a large collection of "
    "scores. For confidence estimation afficionados, please note that this "
    "definition of \"q-value\" is independent of the notion of \"positive FDR\" "
    "as defined in (Storey <em>Annals of Statistics</em> 31:2013-2015:2003).</p>"
    "<p>To estimate FDRs, <code>assign-confidence</code> uses one of two "
    "different procedures. Both require that the input contain both target and "
    "decoy scores. The default, target-decoy competition (TDC) procedure is "
    "described in this article:</p><blockquote>Josh E. Elias and Steve P. Gygi. "
    "\"Target-decoy search strategy for increased confidence in large-scale "
    "protein identifications by mass spectrometry.\" <em>Nature Methods</em>. "
    "4(3):207-14, 2007.</blockquote><p>Note that <code>assign-confidence</code> "
    "implements a variant of the protocol proposed by Elias and Gygi: rather "
    "than reporting a list that contains both targets and decoys, <code>"
    "assign-confidence</code> reports only the targets. The FDR estimate is "
    "adjusted accordingly (by dividing by 2).</p><p>The alternative, <em>"
    "mix-max</em> procedure is described in this article:</p><blockquote>Uri "
    "Keich, Attila Kertesz-Farkas and William Stafford Noble. \"An improved "
    "false discovery rate estimation procedure for shotgun proteomics.\" "
    "Submitted.</blockquote><p>Note that the mix-max procedure requires as "
    "input calibrated scores, such as Comet E-values or p-values produced "
    "using Tide-s <code>exact-p-value</code> option.</p>"
    "<p>The mix-max procedure requires that scores "
    "are reported from separate target and decoy searches. Thus, this approach "
    "is incompatible with a search that is run using the <code>--concat T"
    "</code> option to <code>tide-search</code> or the <code>--decoy_search 2"
    "</code> option to <code>comet</code>. On the other hand, the TDC "
    "procedure can take as input "
    "searches conducted in either mode (concatenated or separate). If given "
    "separate search results and asked to do TDC estimation, <code>"
    "assign-confidence</code> will carry out the target-decoy competition as "
    "part of the confidence estimation procedure.</p><p>In each case, the "
    "estimated FDRs are converted to q-values by sorting the scores then "
    "taking, for each score, the minimum of the current FDR and all of the FDRs "
    "below it in the ranked list.</p><p>A primer on multiple testing correction "
    "can be found here:</p><blockquote>William Stafford Noble. <a href=\""
    "http://www.nature.com/nbt/journal/v27/n12/full/nbt1209-1135.html\">\"How "
    "does multiple testing correction work?\"</a> <em>Nature Biotechnology</em>. "
    "27(12):1135-1137, 2009.</blockquote>]]";
}

/**
 * \returns the command arguments
 */
vector<string> ComputeQValues::getArgs() const {
  string arr[] = {
    "target input+"
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
    "sidak",
    "peptide-level",
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
vector< pair<string, string> > ComputeQValues::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("assign-confidence.target.txt",
    "a <a href=\"../file-formats/txt-format.html\">tab-delimited text file</a> that contains the "
    "targets, sorted by score. The file will contain one new column, named "
    "\"&lt;method&gt; q-value\", where &lt;method&gt; is either \"tdc\" or \"mix-max\"."));
  outputs.push_back(make_pair("assign-confidence.log.txt",
    "a log file containing a copy of all messages that were printed to stderr."));
  outputs.push_back(make_pair("assign-confidence.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
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
