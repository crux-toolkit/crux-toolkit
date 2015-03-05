/**
 * \file CruxApplication.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Abstract Object for a CruxApplication
 *****************************************************************************/
#include "CruxApplication.h"

#include "io/carp.h"
#include "parameter.h"
#include "util/ArgParser.h"
#include "util/crux-file-utils.h"
#include "util/Params.h"
#include "util/WinCrux.h"

#include <iostream>

using namespace std;

/**
 * Frees an allocated CruxApplication
 */
CruxApplication::~CruxApplication() {
}

/**
 * \returns the file stem of the application, default getName.
 */
string CruxApplication::getFileStem() {
  return getName();
}

/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T CruxApplication::getCommand() {
  return MISC_COMMAND;
}


bool CruxApplication::needsOutputDirectory() {
  return false;
}

void CruxApplication::initialize(
  const char** argument_list, ///< list of required arguments
  int num_arguments,          ///< number of elements in arguments_list
  const char** option_list,   ///< list of optional flags
  int num_options,            ///< number of elements in options_list
  int argc,                   ///< number of tokens on cmd line
  char** argv                 ///< array of command line tokens
  ) {
  initializeParams(getName(), argument_list, num_arguments,
                   option_list, num_options, argc, argv);
  processParams();
  Params::Finalize();

  set_verbosity_level(Params::GetInt("verbosity"));

  carp(CARP_INFO, "Beginning %s.", getName().c_str());

  // Set seed for random number generation 
  if (Params::GetString("seed") == "time") {
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    mysrandom((unsigned)seconds); // Convert seconds to a unsigned int
  } else {
    mysrandom(StringUtils::FromString<unsigned>(Params::GetString("seed")));
  }
  
  // Start the timer.
  wall_clock();

  // Create output directory if appliation needs it.
  if (needsOutputDirectory()) {
    // Create output directory 
    string output_folder = Params::GetString("output-dir");
    if (create_output_directory(output_folder, Params::GetBool("overwrite")) == -1) {
      carp(CARP_FATAL, "Unable to create output directory %s.", output_folder.c_str());
    }

    // Open the log file to record carp messages 
    open_log_file(getFileStem() + ".log.txt");
  
    // Store the host name, start date and time, and command line.
    carp(CARP_INFO, "CPU: %s", hostname());
    carp(CARP_INFO, date_and_time());
    log_command_line(argc, argv);

    // Write the parameter file
    string paramFile = make_file_path(getFileStem() + ".params.txt");
    ofstream* file = create_file(paramFile.c_str(), Params::GetBool("overwrite"));
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
bool CruxApplication::hidden() {
  return false;
}

/**
 * Read in all parameters from command line and parameter file
 */
void CruxApplication::initializeParams(
  const string& appName,
  const char** argument_list, ///< list of required arguments
  int num_arguments,          ///< number of elements in arguments_list
  const char** option_list,   ///< list of optional flags
  int num_options,            ///< number of elements in options_list
  int argc,                   ///< number of tokens on cmd line
  char** argv                 ///< array of command line tokens
) {
  initialize_parameters();
  set_verbosity_level(Params::GetInt("verbosity"));

  // Parse command line
  vector<string> argSpecs;
  for (int i = 0; i < num_arguments; i++) {
    argSpecs.push_back(argument_list[i]);
  }
  ArgParser argParser;
  try {
    argParser.Parse(argc, argv, argSpecs);

    // Read parameter file if specified
    string parameter_file = argParser.GetOption("parameter-file");
    if (!parameter_file.empty()) {
      parse_parameter_file(parameter_file.c_str());
      read_mods_from_file(parameter_file.c_str());
    }
    // Process command line options
    const map<string, string>& options = argParser.GetOptions();
    for (map<string, string>::const_iterator i = options.begin(); i != options.end(); i++) {
      Params::Set(i->first, i->second);
    }
    // Process command line arguments
    const map< string, vector<string> >& args = argParser.GetArgs();
    for (map< string, vector<string> >::const_iterator i = args.begin(); i != args.end(); i++) {
      for (vector<string>::const_iterator j = i->second.begin(); j != i->second.end(); j++) {
        Params::AddArgValue(i->first, *j);
      }
    }
  } catch (const runtime_error& e) {
    carp(CARP_FATAL, "%s\n\n%s\n", e.what(),
         getUsage(appName, argument_list, num_arguments, option_list, num_options).c_str());
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
  const char** argument_list, ///< list of required arguments
  int num_arguments,          ///< number of elements in arguments_list
  const char** option_list,   ///< list of optional flags
  int num_options             ///< number of elements in options_list
) {
  vector<string> argDisplay;
  for (int i = 0; i < num_arguments; i++) {
    string argName = argument_list[i];
    if (StringUtils::EndsWith(argName, "+")) {
      argName = argName.substr(0, argName.length() - 1);
      argDisplay.push_back("<" + argName + ">+");
    } else {
      argDisplay.push_back("<" + argName + ">");
    }
  }

  stringstream usage;
  usage << "USAGE:" << endl
        << endl
        << "  crux " << appName << " [options]";
  for (vector<string>::const_iterator i = argDisplay.begin(); i != argDisplay.end(); i++) {
    usage << ' ' << *i;
  }
  usage << endl << endl
        << "REQUIRED ARGUMENTS:";
  for (vector<string>::const_iterator i = argDisplay.begin(); i != argDisplay.end(); i++) {
    stringstream line;
    string argName = *i;
    argName = argName.substr(
      1, argName.length() - (StringUtils::EndsWith(argName, "+") ? 3 : 2));
    line << *i << ' ' << Params::GetUsage(argName);
    usage << endl << endl << StringUtils::LineFormat(line.str(), 80, 2);
  }
  usage << endl << endl
        << "OPTIONAL ARGUMENTS:" << endl;
  for (int i = 0; i < num_options; i++) {
    usage << endl
          << "  [--" << option_list[i] << " <" << Params::GetType(option_list[i])
          << ">]" << endl
          << StringUtils::LineFormat(Params::GetUsage(option_list[i]), 80, 5);
  }
  usage << endl << endl
        << "Additional parameters are documented in the online documentation.";
  return usage.str();
}

