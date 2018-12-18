/**
 * \file CruxApplication.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Abstract Object for a CruxApplication
 *****************************************************************************/
#include "CruxApplication.h"
#include "crux_version.h"
#include "io/carp.h"
#include "parameter.h"
#include "util/ArgParser.h"
#include "util/crux-utils.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "util/GlobalParams.h"
#include "util/WinCrux.h"

#include <iostream>
#define  NO_BOOST_DATE_TIME_INLINE
#include <boost/thread.hpp>

using namespace std;

/**
 * Frees an allocated CruxApplication
 */
CruxApplication::~CruxApplication() {
}

/**
 * \returns the arguments of the application
 */
vector<string> CruxApplication::getArgs() const {
  return vector<string>();
}

/**
 * \returns the options of the application
 */
vector<string> CruxApplication::getOptions() const {
  return vector<string>();
}

/**
 * \returns the outputs of the application as name -> description
 */
vector< pair<string, string> > CruxApplication::getOutputs() const {
  return vector< pair<string, string> >();
}

/**
 * \returns the file stem of the application, default getName.
 */
string CruxApplication::getFileStem() const {
  return getName();
}

/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T CruxApplication::getCommand() const {
  return MISC_COMMAND;
}


bool CruxApplication::needsOutputDirectory() const {
  return false;
}

void CruxApplication::initialize(int argc, char** argv) {
  initializeParams(getName(), getArgs(), getOptions(), argc, argv);

  set_verbosity_level(Params::GetInt("verbosity"));

  // Set seed for random number generation 
  if (Params::GetString("seed") == "time") {
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    mysrandom((unsigned)seconds); // Convert seconds to a unsigned int
  } else {
    mysrandom(StringUtils::FromString<unsigned>(Params::GetString("seed")));
  }

  // Create output directory if appliation needs it.
  if (needsOutputDirectory()) {
    // Create output directory 
    string output_folder = Params::GetString("output-dir");
    if (create_output_directory(output_folder, Params::GetBool("overwrite")) == -1) {
      carp(CARP_FATAL, "Unable to create output directory %s.", output_folder.c_str());
    }

    // Open the log file to record carp messages 
    open_log_file(getFileStem() + ".log.txt");
  
    // Store the host name, start date and time, version number, and command line.
    carp(CARP_INFO, "CPU: %s", hostname());
    carp(CARP_INFO, "Crux version: %s", CRUX_VERSION);
    carp(CARP_INFO, date_and_time());
    log_command_line(argc, argv);
  }

  // Start the timer.
  carp(CARP_INFO, "Beginning %s.", getName().c_str());
  wall_clock();

  processParams();
  Params::Finalize();
  GlobalParams::set();
  if (!Params::GetBool("no-analytics")) {
    // Post data to Google Analytics using a separate thread
    boost::thread analytics_thread(postToAnalytics, getName());
  }

  if (needsOutputDirectory()) {
    // Write the parameter file
    string paramFile = make_file_path(getFileStem() + ".params.txt");
    ofstream* file = FileUtils::GetWriteStream(paramFile, Params::GetBool("overwrite"));
    if (file == NULL) {
      throw runtime_error("Could not open " + paramFile + " for writing");
    }
    Params::Write(file);
    delete file;
  }
}

/**
 * Should this application be kept from the usage statement?
 */
bool CruxApplication::hidden() const {
  return false;
}

/**
 * Read in all parameters from command line and parameter file
 */
void CruxApplication::initializeParams(
  const string& appName,
  const vector<string>& appArgs,
  const vector<string>& appOptions,
  int argc,
  char** argv
) {
  initialize_parameters();
  set_verbosity_level(Params::GetInt("verbosity"));

  // Parse command line
  ArgParser argParser;
  try {
    argParser.Parse(argc, argv, appArgs);

    // Read parameter file if specified
    string parameter_file = argParser.GetOption("parameter-file");
    if (!parameter_file.empty()) {
      parse_parameter_file(parameter_file.c_str());
    }
    // Process command line options
    const map<string, string>& options = argParser.GetOptions();
    for (map<string, string>::const_iterator i = options.begin(); i != options.end(); i++) {
      try {
        Params::Set(i->first, i->second);
      } catch (const runtime_error& e) {
        throw ArgParserException(e.what());
      }
    }
    // Process command line arguments
    const map< string, vector<string> >& args = argParser.GetArgs();
    for (map< string, vector<string> >::const_iterator i = args.begin(); i != args.end(); i++) {
      for (vector<string>::const_iterator j = i->second.begin(); j != i->second.end(); j++) {
        try {
          Params::AddArgValue(i->first, *j);
        } catch (const runtime_error& e) {
          throw ArgParserException(e.what());
        }
      }
    }
  } catch (const ArgParserException& e) {
    carp(CARP_FATAL, "%s\n\n%s\n", e.what(),
         getUsage(appName, appArgs, appOptions, e.ShowFullUsage()).c_str());
  } catch (const runtime_error& e) {
    carp(CARP_FATAL, "%s", e.what());
  }
}

/**
 * Process parameters after they have been set up, but before they have been
 * finalized
 */
void CruxApplication::processParams() {
}

string CruxApplication::getUsage(
  const string& appName,
  const vector<string>& args,
  const vector<string>& options,
  bool full
) {
  vector<string> argDisplay;
  for (vector<string>::const_iterator i = args.begin(); i != args.end(); i++) {
    argDisplay.push_back(StringUtils::EndsWith(*i, "+") ?
      "<" + i->substr(0, i->length() - 1) + ">+" :"<" + *i + ">");
  }

  stringstream usage;
  usage << "USAGE:" << endl
        << endl
        << "  crux " << appName << " [options]";
  for (vector<string>::const_iterator i = argDisplay.begin(); i != argDisplay.end(); i++) {
    usage << ' ' << *i;
  }

  if (full) {
    usage << endl << endl
          << "REQUIRED ARGUMENTS:";
    for (vector<string>::const_iterator i = argDisplay.begin(); i != argDisplay.end(); i++) {
      stringstream line;
      string argName = i->substr(1, i->length() - (StringUtils::EndsWith(*i, "+") ? 3 : 2));
      line << *i << ' ' << Params::ProcessHtmlDocTags(Params::GetUsage(argName));
      usage << endl << endl << StringUtils::LineFormat(line.str(), 80, 2);
    }
    usage << endl << endl
          << "OPTIONAL ARGUMENTS:" << endl;
    for (vector<string>::const_iterator i = options.begin(); i != options.end(); i++) {
      if (!Params::IsVisible(*i)) {
        continue;
      }
      string defaultString = Params::GetStringDefault(*i);
      if (defaultString.empty()) {
        defaultString = "<empty>";
      }
      usage << endl
            << "  [--" << *i << " " << Params::GetAcceptedValues(*i) << "]" << endl
            << StringUtils::LineFormat(Params::ProcessHtmlDocTags(Params::GetUsage(*i)) +
                                       " Default = " + defaultString + ".", 80, 5);
    }
    if (options.empty()) {
      usage << endl
            << "  This command does not support any optional parameters.";
    }
  } else {
    usage << endl << endl
          << "Run \"crux " << appName << "\" with no arguments to see a full "
             "list of options.";
  }
  return usage.str();
}

