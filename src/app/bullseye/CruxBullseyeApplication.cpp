/**
 * \file CruxBullseyeApplication.cpp 
 * \brief Given a ms1 and ms2 file, run hardklor followed by the bullseye algorithm.
 *****************************************************************************/
#include "CruxBullseyeApplication.h"
#include "app/hardklor/CruxHardklorApplication.h"
#include "util/CarpStreamBuf.h"
#include "io/DelimitedFileWriter.h"

#include "util/crux-utils.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"

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

int CruxBullseyeApplication::main(int argc, char** argv) {
  string out_format = Params::GetString("spectrum-format");
  if (out_format.empty()) {
    out_format = "ms2";
  } else if (out_format != "ms2" && out_format != "bms2" &&
             out_format != "cms2" && out_format != "mgf") {
    carp(CARP_FATAL, "spectrum-format must be ms2, bms2, cms2, or mgf, but was %s",
         out_format.c_str());
  }
  return main(
    Params::GetString("MS1 spectra"),
    Params::GetString("MS2 spectra"),
    make_file_path("bullseye.pid." + out_format),
    make_file_path("bullseye.no-pid." + out_format)
  );
}

/**
 * main method for CruxBullseyeApplication
 */
int CruxBullseyeApplication::main(
  const string& input_ms1,
  const string& input_ms2,
  const string& match_ms2,
  const string& nomatch_ms2
) {
  /* Get parameters. */
  string hardklor_output = Params::GetString("hardklor-file");
  if (hardklor_output.empty()) {
    hardklor_output = make_file_path("hardklor.mono.txt");
    if (Params::GetBool("overwrite") || (!FileUtils::Exists(hardklor_output))) {
      carp(CARP_DEBUG, "Calling hardklor");
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
  be_args_vec.push_back(Params::GetString("max-persist"));
  
  if (Params::GetBool("exact-match")) {
    be_args_vec.push_back("-e");
    be_args_vec.push_back("-p");
    be_args_vec.push_back(Params::GetString("exact-tolerance"));
  }
  
  be_args_vec.push_back("-g");
  be_args_vec.push_back(Params::GetString("gap-tolerance"));
  
  be_args_vec.push_back("-r");
  be_args_vec.push_back(Params::GetString("persist-tolerance"));
  
  be_args_vec.push_back("-n");
  be_args_vec.push_back(Params::GetString("bullseye-min-mass"));

  be_args_vec.push_back("-m");
  be_args_vec.push_back(Params::GetString("bullseye-max-mass"));

  be_args_vec.push_back("-s"); 
  be_args_vec.push_back(
    //TODO- I don't know why bullseye.cpp adds 1 to the value passed in...
    StringUtils::ToString(Params::GetInt("scan-tolerance") - 1));
  
  be_args_vec.push_back("-t");
  be_args_vec.push_back(Params::GetString("retention-tolerance"));
  


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

  // Re-route stream to log file
  CarpStreamBuf buffer;
  streambuf* old = cout.rdbuf();
  cout.rdbuf(&buffer);

  /* Call bullseyeMain */
  int ret = bullseyeMain(be_argc, be_argv);

  // Recover stream
  cout.rdbuf(old);

  delete []be_argv;

  return ret;
}

/**
 * \returns the command name for CruxBullseyeApplication
 */
string CruxBullseyeApplication::getName() const {
  return "bullseye";
}

/**
 * \returns the description for CruxBullseyeApplication
 */
string CruxBullseyeApplication::getDescription() const {
  return
    "[[nohtml:Assign high resolution precursor m/z values to MS/MS data using "
    "the Hardklör algorithm.]]"
    "[[html:<p>Bullseye assigns high resolution precursor m/z values to "
    "fragmentation (MS2) spectra. Bullseye uses the Hardkl&ouml;r algorithm to "
    "identify persistent isotope distributions (PPIDs) in precursor (MS1) "
    "scans. For each PPID, MS2 scans that occur within a specified time and "
    "m/z range are assigned the average monoisotopic m/z from the PPID "
    "assigned as the precursor m/z. A detailed description of the Bullseye "
    "algorithm is given in </p><quote>Hsieh EJ, Hoopmann MR, Maclean B, "
    "MacCoss MJ. <a href=\"http://pubs.acs.org/doi/abs/10.1021/pr900816a\">"
    "&quot;Comparison of Database Search Strategies for High Precursor Mass "
    "Accuracy MS/MS Data&quot;</a>. <em>Journal of Proteome Research</em>. "
    "9(2):1138-43, 2010.</quote><p>Note that, in complex samples, it is not "
    "unusual for multiple PPIDs to be found near an MS2 spectrum. In those "
    "cases, Bullseye will assign both mass measurements to the spectrum. In a "
    ".ms2 file, multiple Z line entries will be made for the scan number.</p>"
    "<p>It is possible to reduce the number of scans that receive multiple "
    "PPIDs by adjusting Bullseye's parameters. For example, reducing the "
    "retention time tolerance (&quot;--retention-tolerance&quot;) or reducing "
    "the tolerance for persistent peptides (&quot;--persist-tolerance&quot;) "
    "will reduce the chances of multiple PPIDs being assigned.</p><p>Bullseye "
    "uses Hardkl&ouml;r, so all of the <a href=\"hardklor.html\">Hardkl&ouml;r "
    "parameters</a> may also be used with Bullseye. For users familiar with "
    "the standalone version of Bullseye, the parameter mapping is "
    "<a href=\"bullseye_standalone_to_crux.html\">here</a>.</p>]]";
}

/**
 * \returns the command arguments
 */
vector<string> CruxBullseyeApplication::getArgs() const {
  string arr[] = {
    "MS1 spectra",
    "MS2 spectra"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> CruxBullseyeApplication::getOptions() const {
  string arr[] = {
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
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > CruxBullseyeApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("bullseye.pid.<format>",
    "a file containing the fragmentation spectra for which accurate masses "
    "were successfully inferred. Unless otherwise specified (with the "
    "--spectrum-format option), the output file format is \".ms2\". Note that "
    "if the output format is \".ms2,\" then a single spectrum may have "
    "multiple \"Z\" lines, each indicating a charge state and accurate mass. "
    "In addition, Bullseye inserts an \"I\" line (for charge-dependent "
    "analysis) corresponding to each \"Z\" line. The \"I\" line contains "
    "\"EZ\" in the second column, the charge and mass from the associated "
    "\"Z\" line in the third and fourth colummns, followed by the "
    "chromatographic apex and the intensity at the chromatographic apex."));
  outputs.push_back(make_pair("bullseye.no-pid.<format>",
    "a file containing the fragmentation spectra for which accurate masses "
    "were not inferred."));
  outputs.push_back(make_pair("hardklor.mono.txt",
    "a tab-delimited text file containing one line for each isotope "
    "distribution, as described <a href=\"hardklor.html\">here</a>."));
  outputs.push_back(make_pair("bullseye.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("bullseye.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "standard error."));
  return outputs;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CruxBullseyeApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T CruxBullseyeApplication::getCommand() const {
  return BULLSEYE_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
