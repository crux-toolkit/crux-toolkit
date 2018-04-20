/**
 *\file MakePinApplication.cpp 
 *****************************************************************************/
#include "MakePinApplication.h"
#include "io/PinWriter.h"
#include "parameter.h"
#include "util/Params.h"
#include "io/MatchCollectionParser.h"
#include "io/SQTReader.h"
#include "util/StringUtils.h"
#include <sstream>
#include <iomanip>
#include <ios>

using namespace std;

/**
 * \returns a blank PercolatorApplication object
 */
MakePinApplication::MakePinApplication() {
}

/**
 * Destructor
 */
MakePinApplication::~MakePinApplication() {}

/**
 * main method for MakePinApplication
 */
int MakePinApplication::main(int argc, char** argv) {
  vector<string> target_paths = Params::GetStrings("target input");
  vector<string> search_result_files;
  get_search_result_paths(target_paths, search_result_files);
  return main(search_result_files);
}

/**
 * \runs make-pin application
 */
int MakePinApplication::main(const vector<string>& paths) {
  //create MatchColletion 
  MatchCollectionParser parser;

  if (paths.empty()) {
    carp(CARP_FATAL, "No search paths found!");
  }

  MatchCollection* target_collection = new MatchCollection();
  MatchCollection* decoy_collection = new MatchCollection();

  int max_charge = 0;
  for (vector<string>::const_iterator iter = paths.begin(); iter != paths.end(); ++iter) {
    carp(CARP_INFO, "Parsing %s", iter->c_str());
    if (StringUtils::IEndsWith(*iter, ".sqt")) {
      SQTReader::readSymbols(*iter);
    }

    MatchCollection* current_collection = parser.create(iter->c_str(), "");
    if (!target_collection->getHasDistinctMatches() && current_collection->getHasDistinctMatches()) {
      target_collection->setHasDistinctMatches(true);
      decoy_collection->setHasDistinctMatches(true);
    }
    for (int scorer_idx = (int)SP; scorer_idx < (int)NUMBER_SCORER_TYPES; scorer_idx++) {
      SCORER_TYPE_T cur_type = (SCORER_TYPE_T)scorer_idx;
      bool scored = current_collection->getScoredType(cur_type);
      target_collection->setScoredType(cur_type, scored);
      decoy_collection->setScoredType(cur_type, scored);
    }
    MatchIterator match_iter(current_collection);
    while (match_iter.hasNext()) {
      Crux::Match* match = match_iter.next();
      if (match->getNullPeptide()) {
        decoy_collection->addMatch(match);
      } else {
        target_collection->addMatch(match);
      }
      int charge = match->getCharge();
      if (charge > max_charge) {
        max_charge = charge;
      }
    }
    delete current_collection;
  }

  carp(CARP_INFO, "There are %d target matches and %d decoys",
       target_collection->getMatchTotal(), decoy_collection->getMatchTotal());
  carp(CARP_INFO, "Maximum observed charge is %d.", max_charge);
  if (Params::GetInt("max-charge-feature") != 0) {
    max_charge = Params::GetInt("max-charge-feature");
    carp(CARP_INFO, "Maximum charge feature set to %d.", max_charge);
  }    
  if (target_collection->getMatchTotal() == 0) {
    carp(CARP_FATAL, "No target matches found!");
  } else if (decoy_collection->getMatchTotal() == 0) {
    carp(CARP_FATAL, "No decoy matches found!  Did you set 'decoy-prefix' properly?");
  }

  //prepare output file 
  string output_filename = Params::GetString("output-file");
  if (output_filename.empty()) {
    string fileroot = Params::GetString("fileroot");
    if (!fileroot.empty()) {
      fileroot += ".";
    }
    output_filename = fileroot + "make-pin.pin";
  }
  PinWriter writer;
  writer.openFile(output_filename, Params::GetString("output-dir"),
                  Params::GetBool("overwrite"));

  for (int i = 1; i <= max_charge; i++) {
    writer.setEnabledStatus("Charge" + StringUtils::ToString(i), true);
  }
  writer.setEnabledStatus("deltCn", target_collection->getScoredType(DELTA_CN));
  writer.setEnabledStatus("deltLCn", target_collection->getScoredType(DELTA_LCN));
  bool is_sp = target_collection->getScoredType(SP);
  writer.setEnabledStatus("lnrSp", is_sp);
  writer.setEnabledStatus("Sp", is_sp);
  writer.setEnabledStatus("IonFrac", is_sp);

  bool is_refactored_xcorr = target_collection->getScoredType(TIDE_SEARCH_REFACTORED_XCORR);
  writer.setEnabledStatus("XCorr", !is_refactored_xcorr);
  writer.setEnabledStatus("RefactoredXCorr", is_refactored_xcorr);
  writer.setEnabledStatus("NegLog10PValue",
                          target_collection->getScoredType(TIDE_SEARCH_EXACT_PVAL));
  writer.setEnabledStatus("NegLog10ResEvPValue",
                          target_collection->getScoredType(RESIDUE_EVIDENCE_PVAL));
  writer.setEnabledStatus("NegLog10CombinePValue",
                          target_collection->getScoredType(BOTH_PVALUE));

  if (writer.getEnabledStatus("lnNumSP") && target_collection->getHasDistinctMatches()) {
    writer.setEnabledStatus("lnNumSP", false);
    writer.setEnabledStatus("lnNumDSP", true);
  }

  //write .pin file 
  writer.printHeader();
  writer.write(target_collection, vector<MatchCollection*>(1, decoy_collection),
               Params::GetInt("top-match"));

  delete target_collection;
  delete decoy_collection;

  return 0;
}

/**
 * \returns the command name for PercolatorApplication
 */
string MakePinApplication::getName() const {
  return "make-pin";
}

/**
 * \returns the description for PercolatorApplication
 */
string MakePinApplication::getDescription() const {
  return
    "[[nohtml:Given a set of search results files, generate a pin file for "
    "input to crux percolator]]"
    "[[html:<p>Make-pin is a utility program that combines a collection of "
    "target and decoy peptide-spectrum matches (PSMs) into a single file in "
    "pin format, according to <a href=\"https://github.com/percolator/"
    "percolator/wiki/Interface\">this format</a>. The resulting file can be "
    "provided as input to <code><a href=\"percolator.html\">crux percolator</a>"
    "</code>.</p><p><code>make-pin</code> requires as input two sets of PSMs, "
    "one set derived from matching observed spectra against real "
    "(&quot;target&quot;) peptides and a second set derived from matching the "
    "same spectra against &quot;decoy&quot; peptides. The output file contains, "
    "for each PSM, a set of features for use by the Percolator algorithm. These "
    "features are summarized <a href=\"../file-formats/features.html\">here</a>.</p><p>Note "
    "that, in the stand-alone version of Percolator, the functionality provided "
    "by <code>crux make-pin</code> is incorporated into a program called "
    "<code>sqt2pin</code>. However, a significant difference between <code>crux "
    "percolator</code> and the stand-alone version of the program is that "
    "<code>crux percolator</code> does not require an explicit call to "
    "<code>crux make-pin</code>: if input is provided to <code>crux percolator"
    "</code> in a non-pin format, then the input will be automatically "
    "converted to pin format.</p>]]";
}

/**
 * \returns the command arguments
 */
vector<string> MakePinApplication::getArgs() const {
  string arr[] = {
    "target input+"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> MakePinApplication::getOptions() const {
  string arr[] = {
    "decoy-prefix",
    "fileroot",
    "filestem-prefixes",
    "max-charge-feature",
    "output-dir",
    "output-file",
    "overwrite",
    "parameter-file",
    "top-match",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > MakePinApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("make-pin.pin",
    "a tab-delimited file containing the input target and decoy PSMs in pin "
    "format. This file can be changed to an absolute path (see --output-file "
    "option)."));
  outputs.push_back(make_pair("make-pin.params.txt",
    "a file containing the name and value of all parameters for the current "
    "operation. Not all parameters in the file may have been used in the "
    "operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("make-pin.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "standard error."));
  return outputs;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool MakePinApplication::needsOutputDirectory() const {
  return true;
}


bool MakePinApplication:: hidden() const {
  return false;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
