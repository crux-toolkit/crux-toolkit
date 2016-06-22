/**
 * \file CruxHardklorApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "CruxHardklorApplication.h"
#include "CAveragine.h"
#include "CHardklor.h"
#include "CHardklor2.h"
#include "CMercury8.h"
#include "CModelLibrary.h"
#include "CHardklorParser.h"
#include "util/CarpStreamBuf.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "io/DelimitedFileWriter.h"

using namespace std;

CruxHardklorApplication::CruxHardklorApplication() {
}

CruxHardklorApplication::~CruxHardklorApplication() {
}

int CruxHardklorApplication::main(int argc, char** argv) {
  return main(Params::GetString("spectra"));
}

int CruxHardklorApplication::main(const string& ms1) {
  carp(CARP_INFO, "Hardklor v2.19, April 10 2015");
  carp(CARP_INFO, "Mike Hoopmann, Mike MacCoss");
  carp(CARP_INFO, "Copyright 2007-2015");
  carp(CARP_INFO, "University of Washington");

  string cdm = Params::GetString("cdm");
  if (cdm == "B") {
    cdm = "none";
  } else if (cdm == "F") {
    cdm = "fft";
  } else if (cdm == "P") {
    cdm = "patterson";
  } else if (cdm == "Q") {
    cdm = "quick";
  } else if (cdm == "s") {
    cdm = "senko";
  }
  bool xmlOutput = Params::GetBool("hardklor-xml-output");

  vector<char*> hardklorArgs;
  addArg(&hardklorArgs, "hardklor");
  addArg(&hardklorArgs, "-cmd");
  addArg(&hardklorArgs, "algorithm", Params::GetString("hardklor-algorithm"));
  addArg(&hardklorArgs, "averagine_mod", Params::GetString("averagine-mod"));
  addArg(&hardklorArgs, "boxcar_averaging", Params::GetString("boxcar-averaging"));
  addArg(&hardklorArgs, "boxcar_filter", Params::GetString("boxcar-filter"));
  addArg(&hardklorArgs, "boxcar_filter_ppm", Params::GetString("boxcar-filter-ppm"));
  addArg(&hardklorArgs, "centroided", Params::GetBool("centroided"));
  addArg(&hardklorArgs, "charge_algorithm", cdm);
  addArg(&hardklorArgs, "charge_max", Params::GetString("max-charge"));
  addArg(&hardklorArgs, "charge_min", Params::GetString("min-charge"));
  addArg(&hardklorArgs, "correlation", Params::GetString("corr"));
  addArg(&hardklorArgs, "depth", Params::GetString("depth"));
  addArg(&hardklorArgs, "distribution_area", Params::GetBool("distribution-area"));
  addArg(&hardklorArgs, "hardklor_data", Params::GetString("hardklor-data-file"));
  addArg(&hardklorArgs, "instrument", Params::GetString("instrument"));
  addArg(&hardklorArgs, "isotope_data", Params::GetString("isotope-data-file"));
  addArg(&hardklorArgs, "max_features", Params::GetString("max-features"));
  addArg(&hardklorArgs, "ms_level", Params::GetString("mzxml-filter"));
  addArg(&hardklorArgs, "mz_max", Params::GetString("mz-max"));
  addArg(&hardklorArgs, "mz_min", Params::GetString("mz-min"));
  addArg(&hardklorArgs, "mz_window", Params::GetString("mz-window"));
  addArg(&hardklorArgs, "resolution", Params::GetString("resolution"));
  addArg(&hardklorArgs, "scan_range_max", Params::GetString("scan-range-max"));
  addArg(&hardklorArgs, "scan_range_min", Params::GetString("scan-range-min"));
  addArg(&hardklorArgs, "sensitivity", Params::GetString("sensitivity"));
  addArg(&hardklorArgs, "signal_to_noise", Params::GetString("signal-to-noise"));
  addArg(&hardklorArgs, "smooth", Params::GetString("smooth"));
  addArg(&hardklorArgs, "sn_window", Params::GetString("sn-window"));
  addArg(&hardklorArgs, "static_sn", Params::GetBool("static-sn"));
  addArg(&hardklorArgs, "xml", xmlOutput);

  addArg(&hardklorArgs, ms1);

  string outputFile = "hardklor.mono.";
  outputFile += xmlOutput ? "xml" : "txt";
  outputFile = make_file_path(outputFile);
  addArg(&hardklorArgs, outputFile);

  CHardklorParser hp;
  hp.parseCMD(hardklorArgs.size(), &hardklorArgs[0]);

  for (vector<char*>::iterator i = hardklorArgs.begin(); i != hardklorArgs.end(); i++) {
    delete *i;
  }

  // Create all the output files that will be used
  for (int i = 0; i< hp.size(); i++) {
    const char* out = &hp.queue(i).outFile[0];
    if (FileUtils::Exists(out) && !Params::GetBool("overwrite")) {
      carp(CARP_FATAL, "The file '%s' already exists and cannot be overwritten. "
           "Use --overwrite T to replace or choose a different output file name",
           out);
    }
    fstream fptr(out, ios::out);
  }

  // Re-route stdout/stderr
  CarpStreamBuf buffer;
  streambuf* oldCout = cout.rdbuf();
  streambuf* oldCerr = cerr.rdbuf();
  cout.rdbuf(&buffer);
  cerr.rdbuf(&buffer);

  CAveragine* averagine = new CAveragine(hp.queue(0).MercuryFile, hp.queue(0).HardklorFile);
  CMercury8* mercury = new CMercury8(hp.queue(0).MercuryFile);
  CModelLibrary* models = new CModelLibrary(averagine, mercury);

  CHardklor h(averagine, mercury);
  CHardklor2 h2(averagine, mercury, models);
  vector<CHardklorVariant> pepVariants;
  CHardklorVariant hkv;

  for (int i = 0; i < hp.size(); i++) {
    if (hp.queue(i).algorithm == Version2) {
      pepVariants.clear();
      if (!hp.queue(i).noBase) {
        pepVariants.push_back(hkv);
      }
      for (unsigned j = 0; j < hp.queue(i).variant->size(); j++) {
        pepVariants.push_back(hp.queue(i).variant->at(j));
      }
      models->eraseLibrary();
      models->buildLibrary(hp.queue(i).minCharge, hp.queue(i).maxCharge, pepVariants);
      h2.GoHardklor(hp.queue(i));
    } else {
      h.GoHardklor(hp.queue(i));
    }
  }

  // Recover stdout/stderr
  cout.rdbuf(oldCout);
  cerr.rdbuf(oldCerr);

  delete models;
  delete averagine;
  delete mercury;

  return 0;
}

void CruxHardklorApplication::addArg(
  vector<char*>* args,
  const string& arg
) {
  char* toAdd = new char[arg.length() + 1];
  strcpy(toAdd, arg.c_str());
  args->push_back(toAdd);
}

void CruxHardklorApplication::addArg(
  vector<char*>* args,
  const string& name,
  const string& value) {
  if (!value.empty()) {
    addArg(args, "-" + name);
    addArg(args, value);
  }
}

void CruxHardklorApplication::addArg(
  vector<char*>* args,
  const string& name,
  bool value) {
  addArg(args, name, string(value ? "1" : "0"));
}

/**
 * \returns the command name for CruxHardklorApplication
 */
string CruxHardklorApplication::getName() const {
  return "hardklor";
}

/**
 * \returns the description for CruxHardklorApplication
 */
string CruxHardklorApplication::getDescription() const {
  return
    "[[nohtml:Identify isotopic distributions from high-resolution mass "
    "spectra.]]"
    "[[html:<p>Hardkl&ouml;r analyzes high-resolution mass spectra, "
    "identifying protein or peptide isotope distributions and determining the "
    "corresponding monoisotopic masses and charge states. The algorithm aims "
    "to identify persistence peptide isotope distribution (PPIDs), i.e., "
    "isotope distributions that recur over multiple scans. Hardkl&ouml;r is "
    "specifically designed to handle overlapping isotope distributions in a "
    "single spectrum. A detailed description of the Hardkl&ouml;r algorithm is "
    "given in</p><blockquote>Hoopmann MR, Finney GL and MacCoss MJ. <a href=\""
    "http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2556510/\">&quot;High speed "
    "data reduction, feature selection, and MS/MS spectrum quality assessment "
    "of shotgun proteomics datasets using high resolution mass spectrometry."
    "&quot;</a> <em>Analytical Chemistry</em>. 79:5630-5632 (2007)."
    "</blockquote>]]";
}

/**
 * \returns the command arguments
 */
vector<string> CruxHardklorApplication::getArgs() const {
  string arr[] = {
    "spectra"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> CruxHardklorApplication::getOptions() const {
  string arr[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "hardklor-algorithm",
    "averagine-mod",
    "boxcar-averaging",
    "boxcar-filter",
    "boxcar-filter-ppm",
    "centroided",
    "cdm",
    "min-charge",
    "max-charge",
    "corr",
    "depth",
    "distribution-area",
    "hardklor-data-file",
    "instrument",
    "isotope-data-file",
    "max-features",
    "mzxml-filter",
    "mz-max",
    "mz-min",
    "mz-window",
    "resolution",
    "scan-range-max",
    "scan-range-min",
    "sensitivity",
    "signal-to-noise",
    "smooth",
    "sn-window",
    "static-sn",
    "parameter-file",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > CruxHardklorApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("hardklor.mono.txt",
    "a tab-delimited text file containing one line for each isotope "
    "distribution. The columns appear in the following order:<ol><li><strong>"
    "scan</strong>: The scan number assigned to this spectrum in the input file."
    "</li><li><strong>retention time</strong>: The time (in seconds) at which the "
    "spectrum was collected.</li><li><strong>mass</strong>: The uncharged "
    "monoisotopic mass of the protein or peptide.</li><li><strong>charge</strong>: "
    "The inferred charge state of the protein or peptide.</li><li><strong>intensity"
    "</strong>: The intensity of the base isotope peak of the model used to predict "
    "the protein or peptide.</li><li><strong>m/z</strong>: The m/z of the base peak."
    "</li><li><strong>s/n</strong>: The signal-to-noise threshold, i.e., the relative "
    "abundance a peak must exceed in the spectrum window to be considered in the "
    "scoring algorithm. Note that this is a local noise threshold for the area of "
    "the spectrum that the peptide was identified in.</li><li><strong>modifications"
    "</strong>: Deviations to the averagine model. Only modifications specified by "
    "the user are considered. If no modifications are found in a particular PPID, "
    "then the column is marked with an underscore.</li><li><strong>dotp</strong>: "
    "The dot product score applies to all predictions in a given spectrum window. "
    "Thus, if two protein or peptide predictions share the same spectrum window, "
    "then they have a single dot product score that is the score of their combined "
    "peaks.</li></ol>"));
  outputs.push_back(make_pair("hardklor.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("hardklor.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "stderr."));
  return outputs;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CruxHardklorApplication::needsOutputDirectory() const {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
