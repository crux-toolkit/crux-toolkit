#include "io/carp.h"
#include "bullseye/CruxBullseyeApplication.h"
#include "util/FileUtils.h"
#include "MakePinApplication.h"
#include "PercolatorApplication.h"
#include "Pipeline.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "TideSearchApplication.h"
#include "TideIndexApplication.h"
#include "CometApplication.h"
#include "CruxLFQApplication.h"
#include "SpectralCounts.h"

using namespace std;

PipelineApplication::PipelineApplication() {
}

PipelineApplication::~PipelineApplication() {
}

int PipelineApplication::main(int argc, char** argv) {
  checkParams();

  carp(CARP_INFO, "Running pipeline with the following steps:");
  for (vector<CruxApplication*>::iterator i = apps_.begin(); i != apps_.end(); i++) {
    carp(CARP_INFO, "--> %s", (*i)->getName().c_str());
  }

  vector<string> spectra = Params::GetStrings("mass spectra");
  string database = Params::GetString("peptide source");

  vector<string> resultsFiles;
  while (!apps_.empty()) {
    CruxApplication* cur = apps_.front();
    carp(CARP_INFO, "Running %s...", cur->getName().c_str());
    int ret;
    switch (cur->getCommand()) {
      case BULLSEYE_COMMAND:
        ret = runBullseye(cur, &spectra);
        break;
      case COMET_COMMAND:
      case TIDE_SEARCH_COMMAND:
        ret = runSearch(cur, spectra, database, &resultsFiles);
        break;
      case PERCOLATOR_COMMAND:
        ret = runPostProcessor(cur, resultsFiles);
        break;
      case SPECTRAL_COUNTS_COMMAND:
        ret = runSpectralCounts(cur);
        break;
      case CRUX_LFQ_COMMAND:
        ret = runLFQ(cur, spectra);
        break;
      default:
        carp(CARP_FATAL, "Pipeline is not set up to run command '%s'",
                         cur->getName().c_str());
        break;
    }
    delete cur;
    apps_.erase(apps_.begin());

    if (ret != 0) {
      carp(CARP_FATAL, "Error running %s", cur->getName().c_str());
    }
  }

  return 0;
}

void PipelineApplication::checkParams() {
  string search = Params::GetString("search-engine");
  if (search == "comet") {
    if (Params::GetInt("decoy_search") == 0) {
      carp(CARP_FATAL, "Cannot perform post-processing without decoys. "
                       "Set decoy_search to 1 or 2.");
    }
  }
}

vector<string> PipelineApplication::getExpectedResultsFiles(
  CruxApplication* app,
  const vector<string>& spectra
) {
  string outputBase = make_file_path(app->getName());
  if (app->getCommand() == COMET_COMMAND) {
    if (Params::GetInt("decoy_search") < 1) {
      carp(CARP_WARNING, "Searching without decoys (see decoy_search option)");
    }

    bool multiFile = spectra.size() > 1;
    string outputExt;
    if (Params::GetBool("output_txtfile")) {
      outputExt =".txt";
    } else if (Params::GetBool("output_pepxmlfile")) {
      outputExt =".pep.xml";
    } else if (Params::GetBool("output_sqtfile")) {
      outputExt =".sqt";
    } else if (Params::GetBool("output_percolatorfile") && !multiFile) {
      outputExt =".pin";
    } else {
      carp(CARP_FATAL, "No valid Comet output options enabled");
    }

    vector<string> outputBases;
    if (!multiFile) {
      outputBases.push_back(outputBase);
    } else {
      // TODO Same file stems?
      for (vector<string>::const_iterator i = spectra.begin(); i != spectra.end(); i++) {
        outputBases.push_back(outputBase + "." + FileUtils::Stem(*i));
      }
    }

    vector<string> resultsFiles;
    for (vector<string>::const_iterator i = outputBases.begin(); i != outputBases.end(); i++) {
        if (app->getCommand() == COMET_COMMAND && Params::GetInt("decoy_search") < 2) {
            // Comet doesn't append '.target' to its output files for concatenated search
            resultsFiles.push_back(*i + outputExt);
        } else {
            resultsFiles.push_back(*i + ".target" + outputExt);
        }
      if (Params::GetInt("decoy_search") == 2) {
        resultsFiles.push_back(*i + ".decoy" + outputExt);
      }
    }
    return resultsFiles;
  }

  bool concat = Params::GetBool("concat");
  string outputExt;
  if (Params::GetBool("txt-output")) {
    outputExt = ".txt";
  } else if (Params::GetBool("pepxml-output")) {
    outputExt = ".pep.xml";
  } else if (Params::GetBool("sqt-output")) {
    outputExt = ".sqt";
  } else if (Params::GetBool("mzid-output")) {
    outputExt = ".mzid";
  } else if (Params::GetBool("pin-output")) {
    outputExt = ".pin";
  } else {
    carp(CARP_FATAL, "No valid Tide output options enabled");
  }
  vector<string> resultsFiles;
  if (concat || outputExt == ".mzid" || outputExt == ".pin") {
    resultsFiles.push_back(outputBase + outputExt);
  } else {
    resultsFiles.push_back(outputBase + ".target" + outputExt);
    resultsFiles.push_back(outputBase + ".decoy" + outputExt);
  }
  return resultsFiles;
}

int PipelineApplication::runBullseye(CruxApplication* app, vector<string>* spectra) {
  if (app->getCommand() != BULLSEYE_COMMAND) {
    carp(CARP_FATAL, "Something went wrong.");
  }

  string outFormat = Params::GetString("spectrum-format");
  if (outFormat.empty()) {
    outFormat = "ms2";
  }
  for (vector<string>::iterator i = spectra->begin(); i != spectra->end(); i++) {
    string ms1 = *i;
    if (StringUtils::IEndsWith(ms1, ".ms2") || StringUtils::IEndsWith(ms1, ".cms2")) {
      string ms1Check = ms1.substr(0, ms1.length() - 1) + '1';
      if (FileUtils::Exists(ms1Check)) {
        ms1 = ms1Check;
      }
    }
    string outBase = make_file_path(app->getFileStem() + "." + FileUtils::BaseName(*i));
    string outMatch = outBase + ".pid." + outFormat;
    string outNoMatch = outBase + ".nopid." + outFormat;
    int ret = ((CruxBullseyeApplication*)app)->main(ms1, *i, outMatch, outNoMatch);
    if (ret != 0) {
      carp(CARP_ERROR, "Error running Bullseye on '%s'", i->c_str());
      return ret;
    }
    *i = outMatch;
  }
  return 0;
}

int PipelineApplication::runSearch(
  CruxApplication* app,
  const vector<string>& spectra,
  const string& database,
  vector<string>* resultsFiles
) {
  bool comet = app->getCommand() == COMET_COMMAND;
  bool tide = app->getCommand() == TIDE_SEARCH_COMMAND;
  if (!comet && !tide) {
    carp(CARP_FATAL, "Something went wrong.");
  }

  carp(CARP_INFO, "Search will be run with the following files against database '%s':",
                  database.c_str());
  for (vector<string>::const_iterator i = spectra.begin(); i != spectra.end(); i++) {
    carp(CARP_INFO, "--> %s", i->c_str());
  }

  *resultsFiles = getExpectedResultsFiles(app, spectra);

  if (comet) {
    return ((CometApplication*)app)->main(spectra);
  }
  return ((TideSearchApplication*)app)->main(spectra);
}

int PipelineApplication::runPostProcessor(
  CruxApplication* app,
  const vector<string>& resultsFiles
) {
  if (app->getCommand() != PERCOLATOR_COMMAND) {
    carp(CARP_FATAL, "Pipeline post-processing only supports Percolator, "
                     "but got '%s'.", app->getName().c_str());
  }

  carp(CARP_INFO, "Post-processing will be run using the following files:");
  for (vector<string>::const_iterator i = resultsFiles.begin(); i != resultsFiles.end(); i++) {
    carp(CARP_INFO, "--> %s", i->c_str());
  }

  string pin;
  if (resultsFiles.size() == 1 && StringUtils::IEndsWith(resultsFiles.front(), ".pin")) {
    pin = resultsFiles.front();
  } else {
    // If passed anything but a single pin file, run make-pin
    pin = make_file_path("make-pin.pin");
    carp(CARP_INFO, "Running make-pin");
    if (MakePinApplication::main(resultsFiles) != 0) {
      carp(CARP_FATAL, "make-pin failed. Not running Percolator.");
    }
    carp(CARP_INFO, "Finished make-pin.");
  }
  return ((PercolatorApplication*)app)->main(pin, "", "pipeline");
}

int PipelineApplication::runSpectralCounts(
  CruxApplication* app
) {
  string fileroot = Params::GetString("fileroot");
  if (!fileroot.empty()) {
    fileroot.append(".");
  }
  string filename = FileUtils::Join(Params::GetString("output-dir"), fileroot + "percolator.target.psms.txt");
  return ((SpectralCounts*)app)->main(filename);
}

int PipelineApplication::runLFQ(
  CruxApplication* app,
  const vector<string>& spectra
) {
  if (app->getCommand() != CRUX_LFQ_COMMAND) {
    carp(CARP_FATAL, "Something went wrong.");
  }

  string fileroot = Params::GetString("fileroot");
  if (!fileroot.empty()) {
    fileroot.append(".");
  }
  string psmFile = FileUtils::Join(Params::GetString("output-dir"), fileroot + "percolator.target.psms.txt");
  string specfileReplicates = Params::GetString("specfile-replicates");

  carp(CARP_INFO, "LFQ will be run using PSM file '%s' and the following spectra files:",
                  psmFile.c_str());
  for (vector<string>::const_iterator i = spectra.begin(); i != spectra.end(); i++) {
    carp(CARP_INFO, "--> %s", i->c_str());
  }

  return ((CruxLFQApplication*)app)->main(psmFile, spectra, specfileReplicates);
}

string PipelineApplication::getName() const {
  return "pipeline";
}

string PipelineApplication::getDescription() const {
  return
    "[[nohtml:Runs a series of Crux tools on a protein database and one or more "
    "sets of tandem mass spectra.]]"
    "[[html:<p>Given one or more sets of tandem mass spectra as well as a "
    "protein database, this command runs a series of Crux tools and reports all "
    "of the results in a single output directory. There are four steps in the "
    "pipeline:</p><ol><li><a href=\"bullseye.html\">Bullseye</a> to assign high-"
    "resolution precursor m/z values to MS/MS data. This step is optional.</li>"
    "<li>Database searching using either <a href=\"tide-search.html\">"
    "Tide-search</a> or <a href=\"comet.html\">Comet</a>. The database can be "
    "provided as a file in FASTA format, or additionally, an index as produced "
    "by <a href=\"tide-index.html\">tide-index</a>.</li><li>Post-processing "
    "using <a href=\"percolator.html\">Percolator</a>.</li><li>Pseudo quantitation "
    "using <a href=\"spectral-counts.html\">spectral-counts</a></li><li>Label-free "
    "quantitation using <a href=\"lfq.html\">lfq</a></li></ol>"
    "<p>All of the command line options associated with the individual tools "
    "in the pipeline can be used with the <code>pipeline</code> command.</p>]]";
}

vector<string> PipelineApplication::getArgs() const {
  string arr[] = {
    "mass spectra+",
    "peptide source"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> PipelineApplication::getOptions() const {
  string arr[] = {
    "bullseye",
    "search-engine"
  };
  vector<string> options(arr, arr + sizeof(arr) / sizeof(string));

  addOptionsFrom<CruxBullseyeApplication>(&options);
  addOptionsFrom<TideIndexApplication>(&options);
  addOptionsFrom<TideSearchApplication>(&options);
  addOptionsFrom<CometApplication>(&options);
  addOptionsFrom<PercolatorApplication>(&options);
  addOptionsFrom<SpectralCounts>(&options);
  addOptionsFrom<CruxLFQApplication>(&options);

  // Pipeline always feeds LFQ the Percolator PSMs file (see runLFQ()), so
  // psm-file-format isn't actually a user choice here.
  removeOptionFrom(&options, "psm-file-format");

  return options;
}

vector< pair<string, string> > PipelineApplication::getOutputs() const {
  vector< pair<string, string> > outputs;

  addOutputsFrom<CruxBullseyeApplication>(&outputs);
  addOutputsFrom<TideSearchApplication>(&outputs);
  addOutputsFrom<CometApplication>(&outputs);
  addOutputsFrom<PercolatorApplication>(&outputs);
  addOutputsFrom<SpectralCounts>(&outputs);
  addOutputsFrom<CruxLFQApplication>(&outputs);

  return outputs;
}

string PipelineApplication::getFileStem() const {
  return getName();
}

COMMAND_T PipelineApplication::getCommand() const {
  return PIPELINE_COMMAND;
}

bool PipelineApplication::needsOutputDirectory() const {
  return true;
}

bool PipelineApplication::hidden() const {
  return false;
}

void PipelineApplication::processParams() {
  if (Params::GetBool("bullseye")) {
    apps_.push_back(new CruxBullseyeApplication());
  }

  const vector<string>& spectra = Params::GetStrings("mass spectra");
  if (Params::GetString("search-engine") == "comet") {
    Params::Set("database_name", Params::GetString("peptide source"));
    for (vector<string>::const_iterator i = spectra.begin(); i != spectra.end(); i++) {
      Params::AddArgValue("input spectra", *i);
    }
    apps_.push_back(new CometApplication());
  } else {
    for (vector<string>::const_iterator i = spectra.begin(); i != spectra.end(); i++) {
      Params::AddArgValue("tide spectra file", *i);
    }
    Params::Set("tide database", Params::GetString("peptide source"));
    apps_.push_back(new TideSearchApplication());
  }

  apps_.push_back(new PercolatorApplication());

  if (Params::GetString("protein-database").empty()) {
    Params::Set("protein-database", Params::GetString("peptide source"));
  }
  apps_.push_back(new SpectralCounts());

  // Pipeline always post-processes with Percolator, so LFQ must parse its
  // percolator.target.psms.txt output rather than an assign-confidence file.
  Params::Set("psm-file-format", "percolator");
  apps_.push_back(new CruxLFQApplication());

  for (vector<CruxApplication*>::iterator i = apps_.begin(); i != apps_.end(); i++) {
    (*i)->processParams();
  }
}

