/**
 * \file CascadeSearchApplication.cpp
 * \brief Iterative PSM meta-search via Cascade protocol
 ************************************************************/
#include "CascadeSearchApplication.h"
#include "io/OutputFiles.h"
#include "AssignConfidenceApplication.h"
#include "TideSearchApplication.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "util/FileUtils.h"
#include "stdio.h"
#include "boost/filesystem.hpp"

using namespace std;

int const CascadeSearchApplication::CASCADE_TERMINATION_CONDITION = 20;

/**
 * \returns a blank CascadeSearchApplication object
 */
CascadeSearchApplication::CascadeSearchApplication() {

}

/**
 * Destructor
 */
CascadeSearchApplication::~CascadeSearchApplication() {
}

/**
 * main method for CascadeSearchApplication
 */
int CascadeSearchApplication::main(int argc, char** argv) {
  map<pair<string, unsigned int>, bool>* spectrum_flag = new map<pair<string, unsigned int>, bool>;

  carp(CARP_INFO, "Running cascade-search...");

  string database_string = Params::GetString("database-series");
  vector<string> database_indices = StringUtils::Split(database_string, ',');
  OutputFiles* output = new OutputFiles(this);

  int return_code;
  for (unsigned int cascade_cnt = 0; cascade_cnt < database_indices.size(); ++cascade_cnt) {

    //carry out tide-search
    TideSearchApplication TideSearchProgram;
    TideSearchProgram.setSpectrumFlag(spectrum_flag);
    return_code = TideSearchProgram.main(Params::GetStrings("tide spectra file"), database_indices[cascade_cnt]);
    if (return_code != 0) {
      return return_code;
    }

    //pass the output from Tide-Search to Assign-Confidence
    vector<string> bridge_file_name;
    bridge_file_name.push_back(TideSearchProgram.getOutputFileName());

    //carry out assign confidence
    AssignConfidenceApplication AssignConfidenceProgram;
    AssignConfidenceProgram.setSpectrumFlag(spectrum_flag);
    AssignConfidenceProgram.setIterationCnt(cascade_cnt);
    AssignConfidenceProgram.setOutput(output);
    AssignConfidenceProgram.setIndexName(database_indices[cascade_cnt]);
    AssignConfidenceProgram.setFinalIteration(cascade_cnt + 1 == database_indices.size());

    return_code = AssignConfidenceProgram.main(bridge_file_name);
    if (return_code != 0) {
      return return_code;
    }
    spectrum_flag = AssignConfidenceProgram.getSpectrumFlag();

    //remove tide-search and assign-confidence output files.
    string outputdir = Params::GetString("output-dir");
    RemoveTempFiles(outputdir, TideSearchProgram.getName());
    RemoveTempFiles(outputdir, AssignConfidenceProgram.getName());

    int numAccepted = AssignConfidenceProgram.getAcceptedPSMs();
    if (numAccepted < CASCADE_TERMINATION_CONDITION) {
      carp(CARP_INFO,
           "Terminating search early because only %d PSMs were accepted.",
           numAccepted);
      break;
    }
    carp(CARP_INFO, "Finished cascade-search of database %d.\n", cascade_cnt + 1);

  }
  delete output;

  return 0;
}

/**
 * \returns the command name for CascadeSearchApplication
 */
string CascadeSearchApplication::getName() const {
  return "cascade-search";
}

/**
 * \returns the description for CascadeSearchApplication
 */
string CascadeSearchApplication::getDescription() const {
  return
    "[[nohtml:An iterative procedure for incorporating information about "
    "peptide groups into the database search and confidence estimation "
    "procedure.]]"
    "[[html:<p>Cascade-search is a general procedure for incorporating information about "
    "peptide groups into the database search and confidence estimation procedure. Peptides "
    "may be grouped according to, for example, their enzymatic properties (zero, one, or "
    "two enzymatic termini) or the presence of different types of numbers of variable "
    "modifications. The algorithm works on a series of databases, each corresponding to a "
    "different peptide group. The database is searched in series, and after each search, any "
    "spectrum that is identified with a user-specified confidence threshold is sequestered "
    "from subsequent searches. The full cascade search procedure is described in this article: </p>"
    "<blockquote>Attila Kertesz-Farkas, Uri Keich and William Stafford Noble. "
    "<a href=\"http://pubs.acs.org/doi/abs/10.1021/pr501173s\">\"Tandem mass spectrum "
    "identification via cascaded search.\"</a> <i>Journal of Proteome Research</i>. "
    "14(8):3027-38, 2015. </blockquote>]]";
}

/**
 * \returns the command arguments
 */
vector<string> CascadeSearchApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+",
    "database-series"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}


/**
 * \returns the command options
 */
vector<string> CascadeSearchApplication::getOptions() const {

  string arr[] = {
    "q-value-threshold"
  };
  vector<string> options(arr, arr + sizeof(arr) / sizeof(string));
  addOptionsFrom<AssignConfidenceApplication>(&options);
  addOptionsFrom<TideSearchApplication>(&options);
  removeOptionFrom(&options, "top-match");

  return options;
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > CascadeSearchApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("cascade-search.target.txt",
    "a <a href=\"../file-formats/txt-format.html\">tab-delimited text file</a> containing the "
    "target PSMs accepted at a pre-defined fdr."));
  outputs.push_back(make_pair("cascade-search.log.txt",
    "a log file containing a copy of all messages that were printed to stderr."));
  outputs.push_back(make_pair("cascade-search.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  return outputs;
}

/**
 * \returns the filestem for CascadeSearchApplication
 */
string CascadeSearchApplication::getFileStem() const {
  return "cascade-search";
}

COMMAND_T CascadeSearchApplication::getCommand() const {
  return CASCADE_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not.
 */
bool CascadeSearchApplication::needsOutputDirectory() const {
  return true;
}

void CascadeSearchApplication::processParams() {

  if (Params::GetInt("top-match") != 1 && !Params::IsDefault("top-match")) {
    carp(CARP_WARNING, "Cascade-Search can work with top-match = 1 only.");
  }
  Params::Set("top-match", 1);

  if (Params::GetString("score-function") == "both") {
    Params::Set("score","combined p-value");
  } else if (Params::GetString("score-function") == "residue-evidence") {
    if (Params::GetBool("exact-p-value")) {
      Params::Set("score","res-ev score");
    } else {
      Params::Set("score","res-ev p-value");
    }
  } else if (Params::GetBool("exact-p-value")) {
    Params::Set("score", "exact p-value");
  }

  if (Params::GetBool("pin-output")) {
    carp(CARP_FATAL, "Cascade-Search cannot work with pinxml-output=T.");
  }
  if (Params::GetBool("pepxml-output")) {
    carp(CARP_FATAL, "Cascade-Search cannot work with pepxml-output=T.");
  }
  if (Params::GetBool("mzid-output")) {
    carp(CARP_FATAL, "Cascade-Search cannot work with mzid-output=T.");
  }
  if (Params::GetBool("sqt-output")) {
    carp(CARP_FATAL, "Cascade-Search cannot work with sqt-output=T.");
  }
}

void CascadeSearchApplication::RemoveTempFiles(const string& path, const string& prefix) {
  boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
  for (boost::filesystem::directory_iterator i(path); i != end_itr; ++i) {
    // Skip if not a file
    if (!boost::filesystem::is_regular_file(i->status())) {
      continue;
    }

    string filename = i->path().filename().generic_string();
    if (filename.compare(0, prefix.length(), prefix) == 0) {
      FileUtils::Remove(FileUtils::Join(path, filename));
    }
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
