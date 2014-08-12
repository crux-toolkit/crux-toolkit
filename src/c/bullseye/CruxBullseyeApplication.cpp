/**
 * \file CruxBullseyeApplication.cpp 
 * \brief Given a ms1 and ms2 file, run hardklor followed by the bullseye algorithm.
 *****************************************************************************/
#include "CruxBullseyeApplication.h"
#include "CruxHardklorApplication.h"
#include "DelimitedFileWriter.h"

#include "crux-utils.h"

using namespace std;

/**
 * \returns a blank CruxBullseyeApplication object
 */
CruxBullseyeApplication::CruxBullseyeApplication() {

}

/**
 * Destructor
 */
CruxBullseyeApplication::~CruxBullseyeApplication() {
}

/**
 * main method for CruxBullseyeApplication
 */
int CruxBullseyeApplication::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "max-persist",
    "exact-match",
    "exact-tolerance",
    "persist-tolerance",
    "gap-tolerance",
    "scan-tolerance",
    "bullseye-max-mass",
    "bullseye-min-mass",
    "retention-tolerance",
    "spectrum-format",
    "parameter-file",
    "verbosity"
  };

  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"MS1 spectra", "MS2 spectra"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  /* Get parameters. */
  string hardklor_output;
  string input_ms1 = get_string_parameter("MS1 spectra");
  string input_ms2 = get_string_parameter("MS2 spectra");

  string out_format = get_string_parameter_pointer("spectrum-format");
  if (out_format.empty() || out_format == "__NULL_STR") {
    out_format = "ms2";
  } else if (out_format != "ms2" && out_format != "bms2" &&
             out_format != "cms2" && out_format != "mgf") {
    carp(CARP_FATAL, "spectrum-format must be ms2, bms2, cms2, or mgf, but was %s",
         out_format.c_str());
  }
  string match_ms2 = make_file_path("bullseye.pid." + out_format);
  string nomatch_ms2 = make_file_path("bullseye.no-pid." + out_format);
  bool overwrite = get_boolean_parameter("overwrite");

  hardklor_output = string(get_string_parameter_pointer("hardklor-file"));

  if (hardklor_output == "__NULL_STR") {
    hardklor_output = make_file_path("hardklor.mono.txt");
    if ((overwrite) || (!file_exists(hardklor_output))) {
      carp(CARP_DEBUG,"Calling hardklor");
      bool ret = CruxHardklorApplication::main(input_ms1);
      if (ret != 0) {
        carp(CARP_WARNING, "Hardklor failed:%d", ret);
        return ret;
      }
    }
  }

  /* build argument list */
  vector<string> be_args_vec;
  be_args_vec.push_back("bullseye");
  
  /* add flags */
  
  be_args_vec.push_back("-c");
  be_args_vec.push_back(
    DelimitedFileWriter::to_string(get_double_parameter("max-persist"))
    );
  
  if (get_boolean_parameter("exact-match")) {
    be_args_vec.push_back("-e");
    be_args_vec.push_back("-p");
    be_args_vec.push_back(
      DelimitedFileWriter::to_string(get_double_parameter("exact-tolerance"))
      );
  }
  
  be_args_vec.push_back("-g");
  be_args_vec.push_back(
    DelimitedFileWriter::to_string(get_int_parameter("gap-tolerance"))
    );
  
  be_args_vec.push_back("-r");
  be_args_vec.push_back(
    DelimitedFileWriter::to_string(get_double_parameter("persist-tolerance"))
    );
  
  be_args_vec.push_back("-n");
  be_args_vec.push_back(
    DelimitedFileWriter::to_string(get_double_parameter("bullseye-min-mass"))
    );

  be_args_vec.push_back("-m");
  be_args_vec.push_back(
    DelimitedFileWriter::to_string(get_double_parameter("bullseye-max-mass")));

  be_args_vec.push_back("-s"); 
  be_args_vec.push_back(
    //TODO- I don't know why bullseye.cpp adds 1 to the value passed in...
    DelimitedFileWriter::to_string<int>(get_int_parameter("scan-tolerance")-1)
    );
  
  be_args_vec.push_back("-t");
  be_args_vec.push_back(
    DelimitedFileWriter::to_string(get_double_parameter("retention-tolerance"))
    );
  


  /* add arguments */
  be_args_vec.push_back(hardklor_output);
  be_args_vec.push_back(input_ms2);
  be_args_vec.push_back(match_ms2);
  be_args_vec.push_back(nomatch_ms2);


  /* build argv line */
  int be_argc = be_args_vec.size();

  char** be_argv = new char*[be_argc];

  be_argv[0] = (char*)be_args_vec[0].c_str();
  for (int idx = 1;idx < be_argc ; idx++) {
    be_argv[idx] = (char*)be_args_vec[idx].c_str();
    carp(CARP_DEBUG, "be_argv[%d]=%s", idx, be_argv[idx]);
  }

  /* Call bullseyeMain */
  int ret = bullseyeMain(be_argc, be_argv);

  delete []be_argv;

  return ret;
}

/**
 * \returns the command name for CruxBullseyeApplication
 */
string CruxBullseyeApplication::getName() {
  return "bullseye";
}

/**
 * \returns the description for CruxBullseyeApplication
 */
string CruxBullseyeApplication::getDescription() {

  return "Assign high resolution precursor m/z values to MS/MS data "
    "using the Hardklör algorithm.";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CruxBullseyeApplication::needsOutputDirectory() {
  return true;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
