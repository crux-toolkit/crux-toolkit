/**
 * \file PrintProcessedSpectra.cpp
 *
 * AUTHOR: Barbara Frewen
 * CREATE DATE: September 18, 2009
 *
 * DESCRIPTION: Main method for the print-processed-spectra command.
 *              For every spectrum in an ms2 file, process as for
 *              xcorr and print peaks in ms2 format to new file.
 */

#include "PrintProcessedSpectra.h"
#include "io/SpectrumCollectionFactory.h"
#include "util/Params.h"
#include "util/GlobalParams.h"

using namespace std;

/**
 * \returns a blank PrintProcessedSpectra object
 */
PrintProcessedSpectra::PrintProcessedSpectra() {

}

/**
 * Destructor
 */
PrintProcessedSpectra::~PrintProcessedSpectra() {
}

/**
 * main method for PrintProcessedSpectra
 */
int PrintProcessedSpectra::main(int argc, char** argv) {
  // Get arguments and options
  string input_ms2_name  = Params::GetString("ms2 file");
  string output_ms2_name = Params::GetString("output file");
  output_ms2_name = prefix_fileroot_to_name(output_ms2_name);
  string output_dir = Params::GetString("output-dir");
  bool overwrite = Params::GetBool("overwrite");
  OBSERVED_PREPROCESS_STEP_T stop_after = GlobalParams::getStopAfter();

  // open output file
  create_output_directory(output_dir, overwrite);
  FILE* output_ms2 = create_file_in_path(output_ms2_name, output_dir, overwrite);
  // open input file
  Crux::SpectrumCollection* spectra =
    SpectrumCollectionFactory::create(input_ms2_name);
  if (spectra == NULL) {
    carp(CARP_FATAL, "Could not read spectra from %s.", input_ms2_name.c_str());
  }

  spectra->parse();
  carp(CARP_DEBUG, "Found %d spectra in file.", 
       spectra->getNumSpectra());

  // write header to output file
  //const char* header = spectra->getComment();
  //fprintf(output_ms2, "%s", header);
  fprintf(output_ms2, "H\tComment\tSpectra processed as for Xcorr\n");

  // create iterator for getting spectra
  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);

  if (spectrum_iterator == NULL) {
    carp(CARP_FATAL, "Could create spectrum iterator");
  }

  // loop over all spectra, process, print
  while (spectrum_iterator->hasNext()) {
    SpectrumZState cur_zstate;
    int cur_charge = 0;
    Crux::Spectrum* cur_spectrum = 
      spectrum_iterator->next(cur_zstate);

    cur_charge = cur_zstate.getCharge();
    carp(CARP_DETAILED_INFO, "Processing spectrum %d charge %d.",
         cur_spectrum->getFirstScan(), cur_charge);

    // change the peak values
    FLOAT_T* intensities = NULL;
    int max_mz_bin = 0;
    Scorer::getProcessedPeaks(cur_spectrum, cur_charge, XCORR,
                        &intensities, &max_mz_bin, stop_after);

    // print processed spectrum
    cur_spectrum->printProcessedPeaks(cur_zstate,
                                      intensities, max_mz_bin, output_ms2);
  }

  // close output file
  delete spectra;
  fclose(output_ms2);

  return(0);
}

/**
 * \returns the command name for PrintProcessedSpectra
 */
string PrintProcessedSpectra::getName() const {
  return "print-processed-spectra";
}

/**
 * \returns the description for PrintProcessedSpectra
 */
string PrintProcessedSpectra::getDescription() const {
  return 
    "[[nohtml:Process spectra as for scoring xcorr and print the results to a "
    "file.]]"
    "[[html:<p>Pre-process each spectrum in a given file in preparation for "
    "computing XCorr. The pre-processing steps are described in detail in this "
    "paper:</p><blockquote>J. K. Eng, B. Fischer, J. Grossman and M. J. "
    "MacCoss. <a href=\"http://pubs.acs.org/doi/abs/10.1021/pr800420s\">&quot;A "
    "fast SEQUEST cross correlation algorithm"
    ".&quot;</a> <em>Journal of Proteome Research</em>. "
    "7(10):4598-4602, 2008.</blockquote><p>The output of this program is "
    "equivalent to the spectrum shown in Figure 1D of the above paper.</p>]]";
}

/**
 * \returns the command arguments
 */
vector<string> PrintProcessedSpectra::getArgs() const {
  string arr[] = {
    "ms2 file",
    "output file"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> PrintProcessedSpectra::getOptions() const {
  string arr[] = {
    "stop-after",
    "output-units",
    "spectrum-parser",
    "use-z-line",
    "verbosity",
    "parameter-file", 
    "overwrite"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > PrintProcessedSpectra::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("output file",
    "The name of the file in which the processed spectra will be printed in "
    "MS2 format."));
  return outputs;
}

/**
 * \returns the file stem of the application, default getName.
 */
string PrintProcessedSpectra::getFileStem() const {
  return "processed-spectra";
}

COMMAND_T PrintProcessedSpectra::getCommand() const {
  return PROCESS_SPEC_COMMAND;
}

bool PrintProcessedSpectra::needsOutputDirectory() const {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

