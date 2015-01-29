/**
 * \file CruxApplication.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Abstract Object for a CruxApplication
 *****************************************************************************/
#include "CruxApplication.h"

#include "io/carp.h"
#include "parameter.h"
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

  // Verbosity level for set-up/command line reading 
  set_verbosity_level(CARP_WARNING);

  // Initialize parameter.c and set default values
  initialize_parameters();

  // Define optional and required arguments
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);

  // Parse the command line, including optional params file
  // Includes syntax, type, and bounds checking, dies on error 
  string cmd_name = this->getName();
  char* full_cmd = cat_string("crux ", cmd_name.c_str());

  parse_cmd_line_into_params_hash(argc, argv, cmd_name.c_str());

  free(full_cmd);

  carp(CARP_INFO, "Beginning %s.", 
       this->getName().c_str());

  
  // Set seed for random number generation 
  if (get_string_parameter("seed") == "time") {
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    mysrandom((unsigned)seconds); // Convert seconds to a unsigned int
  }
  else{
    mysrandom((unsigned)atoi(get_string_parameter("seed").c_str()));
  }
  
  // Start the timer.
  wall_clock();

  // Create output directory if appliation needs it.
  if (this->needsOutputDirectory()) {

    // Create output directory 
    string output_folder = get_string_parameter("output-dir");
    bool overwrite = get_boolean_parameter("overwrite");
    int result = create_output_directory(output_folder, overwrite);
    if( result == -1 ){
      carp(CARP_FATAL, "Unable to create output directory %s.", output_folder.c_str());
    }
  
    string cmd_file_name = this->getFileStem();
  
    // Open the log file to record carp messages 
    open_log_file(this->getFileStem() + ".log.txt");
  
    // Store the host name, start date and time, and command line.
    carp(CARP_INFO, "CPU: %s", hostname());
    carp(CARP_INFO, date_and_time());
    log_command_line(argc, argv);
 
    // Write the parameter file
    writeParamFile();
  }
}


/**
 * Should this application be kept from the usage statement?
 */
bool CruxApplication::hidden(){
  return false;
}

/**
 * Writes the parameter file
 */
void CruxApplication::writeParamFile() {
  print_parameter_file(getFileStem() + ".params.txt");
}

