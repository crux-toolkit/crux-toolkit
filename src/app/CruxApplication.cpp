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
#define  BOOST_DATE_TIME_NO_LIB
#include <boost/thread.hpp>

using namespace std;

/* Code added by Rufino*/

/* Iterators to the first and last pointers of the main appications list*/
std::vector<CruxApplication *>::const_iterator firstInd, lastInd;


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
    // Post usage data to Google Analytics 4 using async i/o
    postToGA4(getName());
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

  /* The following code is intended to warn the user for a bad usage of command line options on applications
    Author: Rufino Haroldo Locon
    HSE University, Moscow, Russia
    July 2023
  */
  
    bool no_mach = false; // -> this flag is used to check if a parameter in the command line is accepted by the app pointed by appName
    std::string rejected_options = ""; // -> string used to show the final list of parameters thar are not permitted in the current application
    std::string apps_allowed = ""; // -> string used to show the final list of apps where a parameter rejected is permitted
    std::string separator = ", ";
    std::vector<std::string> next_app_options; // -> string vector to keep the option list from every app in applications list
    bool is_tide_search = false; // -> flag used to check if the current application is tide-search
    
    carp(CARP_DEBUG, "Total Options in command line: %s", std::to_string(options.size()).c_str()); // message used to know how many options are in the command line

    if( options.size() != 0 ) { // if there are more than 0 options in command line

      carp(CARP_DEBUG, "Default values will be overwritten for:"); // message used to indicate in wich options the default value will be overwrite 

      for (map<string, string>::const_iterator options_index = options.begin(); options_index != options.end(); options_index++) {
        carp(CARP_DEBUG, "->: %s", options_index->first.c_str());
      }
      for (vector<string>::const_iterator it = appOptions.begin(); it != appOptions.end(); ++it) {
        carp(CARP_DEBUG, "%s options->: %s", appName.c_str(), (*it).c_str());
      }

      /* This loop is intended to check the parameters of the current application against all the appplications in application list*/
      for (map<string, string>::const_iterator options_index = options.begin(); options_index != options.end(); options_index++) {

        if( !(std::find(appOptions.begin(), appOptions.end(), options_index->first) != appOptions.end()) ) { // if one parameter is not allowed for the current app,
                                                                                                            // then check in the list of applications

          // The no-analytics and seed options are available to all
          // applications and handled somewhat differently,
          // so don't check.
          if( options_index->first.compare("no-analytics") == 0 
              || options_index->first.compare("seed") == 0 ) {
            break;
          }

          for ( vector<CruxApplication *>::const_iterator currentApplication = firstInd; currentApplication != lastInd; currentApplication++ ) {
            //in this loop, a search of the parameter is made on all options for all apps

            next_app_options = (*currentApplication)->getOptions();

            if( (std::find(next_app_options.begin(), next_app_options.end(), options_index->first) != next_app_options.end()) ){

              if (appName.compare("tide-search") == 0 && (*currentApplication)->getName().compare("tide-index") == 0) {
                // if the current app is tide-search, and the index in pointed to tide-index, the parameter is then allowed for use in the current app.
                is_tide_search = true;
              } else {

                if(apps_allowed.empty()) {

                  apps_allowed = apps_allowed + (*currentApplication)->getName();

                } else {

                  apps_allowed = apps_allowed + separator + (*currentApplication)->getName();
                }
              }
            }
          }

          rejected_options = rejected_options + "\t-> " + options_index->first +" used in " + apps_allowed + "\n";
          apps_allowed = "";

          if (!is_tide_search) {
            no_mach = true;
          }
        }
      }

      if ( no_mach ){
        carp(CARP_FATAL, "Options found that don't match %s application:\n%s", appName.c_str(), rejected_options.c_str());
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
      "<" + i->substr(0, i->length() - 1) + ">" :"<" + *i + ">");
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

/* Code added by Rufino*/

/* setListApp method receives two parameters, iterator f (first) and interator l (last), those iterators
    come from the crux-main.cpp entry point after application list is initialized. Those parameters are
    assigned to the local interators firstIndex (pointer to the first entry of the list) and q (pointer to the first entry of the list)*/

void CruxApplication::setApplicationsList(std::vector<CruxApplication *>::const_iterator firstIndex, std::vector<CruxApplication *>::const_iterator lastIndex){

  firstInd = firstIndex;
  lastInd = lastIndex;

}
