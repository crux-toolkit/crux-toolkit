/**
 *\file MakePinApplication.cpp 
 *****************************************************************************/
#include "MakePinApplication.h"
#include "io/PinWriter.h"
#include "parameter.h"
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
  string target_path = get_string_parameter("target input");

  vector<string> search_result_files;
  get_search_result_paths(target_path, search_result_files);

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

  for (vector<string>::const_iterator iter = paths.begin(); iter != paths.end(); ++iter) {
    carp(CARP_INFO, "Parsing %s", iter->c_str());
    if (StringUtils::IEndsWith(*iter, ".sqt")) {
      SQTReader::readSymbols(*iter);
    }

    MatchCollection* current_collection = parser.create(iter->c_str(), "");
    for (int scorer_idx = (int)SP; scorer_idx < (int)NUMBER_SCORER_TYPES; scorer_idx++) {
      target_collection->setScoredType((SCORER_TYPE_T)scorer_idx, 
        current_collection->getScoredType((SCORER_TYPE_T)scorer_idx));
      decoy_collection->setScoredType((SCORER_TYPE_T)scorer_idx,
        current_collection->getScoredType((SCORER_TYPE_T)scorer_idx));
    } 
    MatchIterator* match_iter = new MatchIterator(current_collection);
    while(match_iter->hasNext()) {
      Crux::Match* match = match_iter->next();
      if (match->getNullPeptide()) {
        decoy_collection->addMatch(match);
      } else {
        target_collection->addMatch(match);
      }
    }
    delete match_iter;
    delete current_collection;
  }

  carp(CARP_INFO, "There are %d target matches and %d decoys",target_collection->getMatchTotal(), decoy_collection->getMatchTotal());
  if (target_collection->getMatchTotal() == 0) {
    carp(CARP_FATAL, "No target matches found!");
  }
  if (decoy_collection->getMatchTotal() == 0) {
    carp(CARP_FATAL, "No decoy matches found!  Did you set 'decoy-prefix' properly?");
  }

  PinWriter* writer = new PinWriter();

  string output_dir = get_string_parameter("output-dir");
 
  //perpare output file 
  string output_filename = get_string_parameter("output-file");
  if (output_filename.empty()) {
    string fileroot = get_string_parameter("fileroot");
    if (!fileroot.empty()) {
      fileroot += ".";
    }
    output_filename = fileroot + "make-pin.pin";
  }
  writer->openFile(output_filename.c_str(),output_dir.c_str(), get_boolean_parameter("overwrite"));

 //write .pin file 
  vector<MatchCollection*> decoys;
  decoys.push_back(decoy_collection);
  writer->write(target_collection, decoys, get_int_parameter("top-match"));

  //close file 
  writer->closeFile();

  delete target_collection;
  delete decoy_collection;
  delete writer;

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
    "features are summarized <a href=\"features.html\">here</a>.</p><p>Note "
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
    "target input"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> MakePinApplication::getOptions() const {
  string arr[] = {
    "top-match",
    "list-of-files",
    "decoy-prefix",
    "fileroot",
    "output-dir",
    "output-file",
    "overwrite",
    "parameter-file",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
map<string, string> MakePinApplication::getOutputs() const {
  map<string, string> outputs;
  outputs["make-pin.pin"] =
    "a tab-delimited file containing the input target and decoy PSMs in pin "
    "format. This file can be changed to an absolute path (see --output-file "
    "option).";
  outputs["make-pin.params.txt"] =
    "a file containing the name and value of all parameters for the current "
    "operation. Not all parameters in the file may have been used in the "
    "operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs.";
  outputs["make-pin.log.txt"] =
    "a log file containing a copy of all messages that were printed to "
    "standard error.";
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
