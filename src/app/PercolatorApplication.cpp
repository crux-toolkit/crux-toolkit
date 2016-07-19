/**
 * \file PercolatorApplication.cpp 
 * \brief Runs Percolator
 *****************************************************************************/
#include "MakePinApplication.h"
#include "PercolatorApplication.h"
#include "PercolatorAdapter.h"
#include "Caller.h"
#include "util/Params.h"
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ios>
#include "util/CarpStreamBuf.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "io/MzIdentMLWriter.h"
#include "model/ProteinMatchCollection.h"
#include "io/PMCDelimitedFileWriter.h"
#include "io/PMCPepXMLWriter.h"
#include "io/PMCSQTWriter.h"


using namespace std;

/**
 * \returns a blank PercolatorApplication object
 */
PercolatorApplication::PercolatorApplication() {
}

/**
 * Destructor
 */
PercolatorApplication::~PercolatorApplication() {
}

/**
 * main method for PercolatorApplication
 */
int PercolatorApplication::main(int argc, char** argv) {
  string input_pin = Params::GetString("pin");
  carp(CARP_INFO, "Reading file %s", input_pin.c_str());

  // Check if we need to run make-pin first
  if (!Params::GetBool("feature-file-in") &&
      (Params::GetBool("list-of-files") ||
       StringUtils::IEndsWith(input_pin, ".txt") ||
       StringUtils::IEndsWith(input_pin, ".sqt") ||
       StringUtils::IEndsWith(input_pin, ".pep.xml") ||
       StringUtils::IEndsWith(input_pin, ".mzid"))) {
    vector<string> result_files;
    get_search_result_paths(input_pin, result_files);

    input_pin = make_file_path("make-pin.pin");

    vector<string>::const_iterator fileIter = result_files.begin();
    if (StringUtils::IEndsWith(*fileIter, ".pin")) {
      if (FileUtils::Exists(input_pin)) {
        if (Params::GetBool("overwrite")) {
          FileUtils::Remove(input_pin);
        } else {
          carp(CARP_FATAL, "The file '%s' already exists and cannot be overwritten. "
               "Use --overwrite T to replace or choose a different output file name",
               input_pin.c_str());
        }
      }
      FileUtils::Copy(*fileIter, input_pin);
      string headers;
      fstream out(input_pin.c_str());
      if (!out.good()) {
        carp(CARP_FATAL, "Filestream error '%s'", input_pin.c_str());
      }
      getline(out, headers);
      out.seekp(0, ios_base::end);
      for (fileIter++; fileIter != result_files.end(); fileIter++) {
        if (!StringUtils::IEndsWith(*fileIter, ".pin")) {
          FileUtils::Remove(input_pin);
          carp(CARP_FATAL, "Cannot mix .pin with non-pin files");
        }
        ifstream in(fileIter->c_str());
        if (!in.good()) {
          FileUtils::Remove(input_pin);
          carp(CARP_FATAL, "Error opening file '%s' for reading", fileIter->c_str());
        }
        string inLine;
        getline(in, inLine);
        if (headers != inLine) {
          FileUtils::Remove(input_pin);
          carp(CARP_FATAL, "Headers in pin file '%s' were '%s', but expected '%s'",
               fileIter->c_str(), inLine.c_str(), headers.c_str());
        }
        while (!in.eof()) {
          getline(in, inLine);
          out << inLine;
          if (in.peek() != EOF) {
            out << endl;
          }
        }
      }
    } else {
      carp(CARP_INFO, "Running make-pin");
      if (MakePinApplication::main(result_files) != 0 || !FileUtils::Exists(input_pin)) {
        carp(CARP_FATAL, "make-pin failed. Not running Percolator.");
      }
      carp(CARP_INFO, "Finished make-pin.");
    }
  }
  return main(input_pin);
}

/**
 * \brief runs percolator on the input pin
 * \returns whether percolator was successful or not
 */
int PercolatorApplication::main(
  const string& input_pin ///< file path of pin to process.
  ) {
  /* build argument list */
  vector<string> perc_args_vec;
  perc_args_vec.push_back("percolator");

  string output_target_peptides = make_file_path(getFileStem() + ".target.peptides.txt");
  string output_target_psms = make_file_path(getFileStem() + ".target.psms.txt");
  string output_target_proteins = make_file_path(getFileStem() + ".target.proteins.txt");
  string output_decoy_peptides = make_file_path(getFileStem() + ".decoy.peptides.txt");
  string output_decoy_psms = make_file_path(getFileStem() + ".decoy.psms.txt");
  string output_decoy_proteins = make_file_path(getFileStem() + ".decoy.proteins.txt");
  // Target peptides file is written to prevent writing to stdout
  perc_args_vec.push_back("-r");
  perc_args_vec.push_back(output_target_peptides);
  if (Params::GetBool("original-output")) {
    perc_args_vec.push_back("-B");
    perc_args_vec.push_back(output_decoy_peptides);
    perc_args_vec.push_back("-m");
    perc_args_vec.push_back(output_target_psms);
    perc_args_vec.push_back("-M");
    perc_args_vec.push_back(output_decoy_psms);
  }

  //add verbosity
  perc_args_vec.push_back("-v");
  int verbosity = get_verbosity_level();
  if (verbosity <= CARP_FATAL) {
    perc_args_vec.push_back("0");
  } else if (verbosity <= CARP_ERROR) {
    perc_args_vec.push_back("1");
  } else if (verbosity <= CARP_WARNING) {
    perc_args_vec.push_back("1");
  } else if (verbosity <= CARP_INFO) {
    perc_args_vec.push_back("2");
  } else if (verbosity <= CARP_DETAILED_INFO) {
    perc_args_vec.push_back("3");
  } else if (verbosity <= CARP_DEBUG) {
    perc_args_vec.push_back("4");
  } else if (verbosity <= CARP_DETAILED_DEBUG) {
    perc_args_vec.push_back("5");
  } else if (verbosity <= CARP_MAX) {
    perc_args_vec.push_back("5");
  }

  perc_args_vec.push_back("-P");
  string decoy_pre = Params::GetString("decoy-prefix");
  perc_args_vec.push_back(!decoy_pre.empty() ? decoy_pre : "random_");

  string seed_parameter = Params::GetString("percolator-seed");
  unsigned int seed_value;
  if (seed_parameter == "time") {
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    // percolator accepts seed values 1-20000
    seed_value = (unsigned int)seconds % 20000 + 1;
  } else if ((seed_value = StringUtils::FromString<unsigned>(seed_parameter)) == 0) {
    // seed 0 causes segfault in percolator 
    ++seed_value;
  }
  perc_args_vec.push_back("--seed");
  perc_args_vec.push_back(StringUtils::ToString(seed_value));

  if (Params::GetBool("pout-output")) {
    perc_args_vec.push_back("-X");
    perc_args_vec.push_back(make_file_path(getFileStem() + ".pout.xml"));
    if (Params::GetBool("decoy-xml-output")) {
      perc_args_vec.push_back("-Z");
    }
  }

  perc_args_vec.push_back("-p");
  perc_args_vec.push_back(Params::GetString("c-pos"));
 
  perc_args_vec.push_back("-n");
  perc_args_vec.push_back(Params::GetString("c-neg"));
 
  perc_args_vec.push_back("--trainFDR");
  perc_args_vec.push_back(Params::GetString("train-fdr"));
 
  perc_args_vec.push_back("--testFDR");
  perc_args_vec.push_back(Params::GetString("test-fdr"));

  perc_args_vec.push_back("--maxiter");
  perc_args_vec.push_back(Params::GetString("maxiter"));

  if (Params::GetBool("quick-validation")) {
    perc_args_vec.push_back("--quick-validation");
  }


  if (Params::GetBool("feature-file-out")) {
    perc_args_vec.push_back("--tab-out");
    perc_args_vec.push_back(make_file_path(getFileStem() + ".feature.txt"));
  }

  if (Params::GetBool("output-weights")) {
    perc_args_vec.push_back("--weights");
    perc_args_vec.push_back(make_file_path(getFileStem() + ".weights.txt"));
  }
  
  if (!Params::GetString("init-weights").empty()) {
    perc_args_vec.push_back("--init-weights");
    perc_args_vec.push_back(Params::GetString("init-weights"));
  }

  if (!Params::GetString("default-direction").empty()) {  
    perc_args_vec.push_back("--default-direction");
    perc_args_vec.push_back(Params::GetString("default-direction"));
  }

  if (Params::GetBool("unitnorm")) {
    perc_args_vec.push_back("-u");
  }

  if (Params::GetBool("test-each-iteration")) {
    perc_args_vec.push_back("--test-each-iteration");
  }

  if (Params::GetBool("override")) {
    perc_args_vec.push_back("--override");
  }
 
  if (Params::GetBool("klammer")) {
    perc_args_vec.push_back("--klammer");
  }

  /* --doc option disabled, need retention times in pin file
  int doc_parameter = Params::GetInt("doc");
  if(doc_parameter >= 0) {
    perc_args_vec.push_back("--doc");
    perc_args_vec.push_back(to_string(doc_parameter));
  }
  */

  // FIXME include schema as part of distribution and add option to turn on validation
  perc_args_vec.push_back("-s");

  bool set_protein = Params::GetBool("protein");
  if (set_protein) {
    perc_args_vec.push_back("-A");

    if (Params::GetDouble("fido-alpha") > 0) {
      perc_args_vec.push_back("--fido-alpha");
      perc_args_vec.push_back(Params::GetString("fido-alpha"));
    }
    if (Params::GetDouble("fido-beta") > 0) {
      perc_args_vec.push_back("--fido-beta");
      perc_args_vec.push_back(Params::GetString("fido-beta"));
    }
    if (Params::GetDouble("fido-gamma") > 0) {
      perc_args_vec.push_back("--fido-gamma");
      perc_args_vec.push_back(Params::GetString("fido-gamma"));
    }

    if (Params::GetBool("fido-protein-level-pi0")) {
      perc_args_vec.push_back("--fido-protein-level-pi0");
    }
    if (Params::GetBool("fido-empirical-protein-q")) {
       perc_args_vec.push_back("--fido-empirical-protein-q");
    }

    perc_args_vec.push_back("--fido-gridsearch-depth");
    perc_args_vec.push_back(Params::GetString("fido-gridsearch-depth"));


    perc_args_vec.push_back("--fido-fast-gridsearch");
    perc_args_vec.push_back(Params::GetString("fido-fast-gridsearch"));

    perc_args_vec.push_back("--fido-protein-truncation-threshold");
    perc_args_vec.push_back(Params::GetString("fido-protein-truncation-threshold"));

    if (Params::GetBool("fido-split-large-components")) {
      perc_args_vec.push_back("--fido-split-large-components");
    }

    if (Params::GetBool("post-processing-qvality")) {
      perc_args_vec.push_back("--post-processing-qvality");
    }

    if (Params::GetBool("post-processing-tdc")) {
      perc_args_vec.push_back("--post-processing-tdc");
    }

    perc_args_vec.push_back("--fido-gridsearch-mse-threshold");
    perc_args_vec.push_back(Params::GetString("fido-gridsearch-mse-threshold"));

    // Target proteins file is written to prevent writing to stdout
    perc_args_vec.push_back("-l");
    perc_args_vec.push_back(output_target_proteins);
    if (Params::GetBool("original-output")) {
      perc_args_vec.push_back("-L");
      perc_args_vec.push_back(output_decoy_proteins);
    }
  }

  perc_args_vec.push_back(input_pin);

  /* build argv line */

  string perc_cmd;
  vector<const char*> perc_argv;
  for (vector<string>::const_iterator i = perc_args_vec.begin();
       i != perc_args_vec.end();
       i++) {
    perc_argv.push_back(i->c_str());
    perc_cmd += " " + *i;
  }

  carp(CARP_DEBUG, "cmd:%s", perc_cmd.c_str());
  
  /* Re-route stdeer to log file. */
  CarpStreamBuf buffer;
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  /* Call percolatorMain */
  PercolatorAdapter pCaller;
  int retVal = -1;
  if (pCaller.parseOptions(perc_args_vec.size(), (char**)&perc_argv.front())) {
    retVal = pCaller.run();
  } 
  
  if (retVal != 0) {
    carp(CARP_FATAL, "Error running percolator:%d", retVal);
  }

  // get percolator score information into crux objects
  ProteinMatchCollection* target_pmc = pCaller.getProteinMatchCollection();
  ProteinMatchCollection* decoy_pmc = pCaller.getDecoyProteinMatchCollection();
  if (target_pmc == NULL || decoy_pmc == NULL) {
    carp(CARP_WARNING, "Failed translating Percolator objects into Crux objects");
  }

  string output_dir = Params::GetString("output-dir");

  // write txt
  if (!Params::GetBool("original-output")) {
    FileUtils::Remove(output_target_peptides);
    if (set_protein) {
      FileUtils::Remove(output_target_proteins);
    }

    if (Params::GetBool("txt-output")) {
      PMCDelimitedFileWriter txt_writer;
      txt_writer.writeFile(this, output_target_psms,
                           PMCDelimitedFileWriter::PSMS, target_pmc);
      txt_writer.writeFile(this, output_decoy_psms,
                           PMCDelimitedFileWriter::PSMS, decoy_pmc);
      txt_writer.writeFile(this, output_target_peptides,
                           PMCDelimitedFileWriter::PEPTIDES, target_pmc);
      txt_writer.writeFile(this, output_decoy_peptides,
                           PMCDelimitedFileWriter::PEPTIDES, decoy_pmc);

      if (set_protein) {
        txt_writer.writeFile(this, output_target_proteins,
                             PMCDelimitedFileWriter::PROTEINS, target_pmc);
        txt_writer.writeFile(this, output_decoy_proteins,
                             PMCDelimitedFileWriter::PROTEINS, decoy_pmc);
      }
    }
  }

  // write mzid
  if (Params::GetBool("mzid-output")) {
    MzIdentMLWriter mzid_writer, decoy_mzid_writer;
    string mzid_path = make_file_path(getFileStem() + ".target.mzid");
    mzid_writer.openFile(mzid_path, Params::GetBool("overwrite"));
    mzid_writer.addProteinMatches(target_pmc);
    mzid_writer.closeFile();
    mzid_path = make_file_path(getFileStem() + ".decoy.mzid");
    decoy_mzid_writer.openFile(mzid_path, Params::GetBool("overwrite"));
    decoy_mzid_writer.addProteinMatches(decoy_pmc);
    decoy_mzid_writer.closeFile();
  }
  
  // write pepxml
  if (Params::GetBool("pepxml-output")) {
    PMCPepXMLWriter pep_writer;
    string pep_path = make_file_path(getFileStem() + ".target.pep.xml");
    pep_writer.openFile(pep_path.c_str(), Params::GetBool("overwrite"));
    pep_writer.write(target_pmc);
    pep_writer.closeFile();
    pep_path = make_file_path(getFileStem() + ".decoy.pep.xml");
    pep_writer.openFile(pep_path.c_str(), Params::GetBool("overwrite"));
    pep_writer.write(decoy_pmc);
    pep_writer.closeFile();
  }

  Globals::clean();

  /* Recover stderr */
  std::cerr.rdbuf(old);

  return retVal;
}

COMMAND_T PercolatorApplication::getCommand() const {
  return PERCOLATOR_COMMAND;

}

/**
 * \returns the command name for PercolatorApplication
 */
string PercolatorApplication::getName() const {
  return "percolator";
}

/**
 * \returns the description for PercolatorApplication
 */
string PercolatorApplication::getDescription() const {
  return
    "[[nohtml:Re-rank a collection of PSMs using the Percolator algorithm. "
    "Optionally, also produce protein rankings using the Fido algorithm.]]"
    "[[html:<p>Percolator is a semi-supervised learning algorithm that "
    "dynamically learns to separate target from decoy peptide-spectrum matches "
    "(PSMs). The algorithm is described in this article:</p><blockquote> Lukas "
    "K&auml;ll, Jesse Canterbury, Jason Weston, William Stafford Noble and "
    "Michael J. MacCoss. <a href=\"http://noble.gs.washington.edu/papers/"
    "kall2007semi-supervised.html\">&quot;Semi-supervised learning for peptide "
    "identification from shotgun proteomics datasets.&quot;</a> <em>Nature "
    "Methods</em>. 4(11):923-925, 2007.</blockquote><p>Percolator requires as "
    "input two collections of PSMs, one set derived from matching observed "
    "spectra against real (&quot;target&quot;) peptides, and a second derived "
    "from matching the same spectra against &quot;decoy&quot; peptides. The "
    "output consists of ranked lists of PSMs, peptides and proteins. Peptides "
    "and proteins are assigned two types of statistical confidence estimates: "
    "q-values and posterior error probabilities.</p><p>The features used by "
    "Percolator to represent each PSM are summarized <a href=\"../file-formats/features.html\">"
    "here</a>.</p><p>Percolator also includes code from <a href=\""
    "http://noble.gs.washington.edu/proj/fido/\">Fido</a>, whch performs "
    "protein-level inference. The Fido algorithm is described in this article:"
    "</p><blockquote>Oliver Serang, Michael J. MacCoss and William Stafford "
    "Noble. <a href=\"http://pubs.acs.org/doi/abs/10.1021/pr100594k\">"
    "&quot;Efficient marginalization to compute protein posterior probabilities "
    "from shotgun mass spectrometry data.&quot;</a> <em>Journal of Proteome "
    "Research</em>. 9(10):5346-5357, 2010.</blockquote><p>Crux includes code "
    "from <a href=\"http://per-colator.com/\">Percolator</a>. Crux Percolator "
    "differs from the stand-alone version of Percolator in the following "
    "respects:</p><ul><li>In addition to the native Percolator XML file "
    "format, Crux Percolator supports additional input file formats (SQT, "
    "PepXML, tab-delimited text) and output file formats (PepXML, mzIdentML, "
    "tab-delimited text).</li><li>To maintain consistency with the rest of the "
    "Crux commands, Crux Percolator uses different parameter syntax than the "
    "stand-alone version of Percolator.</li><li>Like the rest of the Crux "
    "commands, Crux Percolator writes its files to an output directory, logs "
    "all standard error messages to a log file, and is capable of reading "
    "parameters from a parameter file.</li></ul>]]";
}

/**
 * \returns the command arguments
 */
vector<string> PercolatorApplication::getArgs() const {
  string arr[] = {
    "pin"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> PercolatorApplication::getOptions() const {
  string arr[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "txt-output",
    "pout-output",
    "mzid-output",
    "pepxml-output",
    "feature-file-out",
    "list-of-files",
    "parameter-file",
    "feature-file-in",
    "protein",
    "decoy-xml-output",
    "decoy-prefix",
    "c-pos",
    "c-neg",
    "train-fdr",
    "test-fdr",
    "maxiter",
    "quick-validation",
    "output-weights",
    "init-weights",
    "default-direction",
    "unitnorm",
    "fido-alpha",
    "fido-beta",
    "fido-gamma",
    "test-each-iteration",
    "override",
    "percolator-seed",
    "klammer",
    "only-psms",
    //"doc",
    "fido-protein-level-pi0",
    "fido-empirical-protein-q",
    "fido-gridsearch-depth",
    "post-processing-tdc",
    "fido-gridsearch-mse-threshold",
    "fido-fast-gridsearch",
    "fido-protein-truncation-threshold",
    "fido-split-large-components",
    "post-processing-qvality",
    "verbosity",
    "top-match"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > PercolatorApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("percolator.target.proteins.txt",
    "a tab-delimited file containing the target protein matches. See "
    "<a href=\"../file-formats/txt-format.html\">here</a> for a list of the fields."));
  outputs.push_back(make_pair("percolator.decoy.proteins.txt",
    "a tab-delimited file containing the decoy protein matches. See "
    "<a href=\"../file-formats/txt-format.html\">here</a> for a list of the fields."));
  outputs.push_back(make_pair("percolator.target.peptides.txt",
    "a tab-delimited file containing the target peptide matches. See "
    "<a href=\"../file-formats/txt-format.html\">here</a> for a list of the fields."));
  outputs.push_back(make_pair("percolator.decoy.peptides.txt",
    "a tab-delimited file containing the decoy peptide matches. See "
    "<a href=\"../file-formats/txt-format.html\">here</a> for a list of the fields."));
  outputs.push_back(make_pair("percolator.target.psms.txt",
    "a tab-delimited file containing the target PSMs. See "
    "<a href=\"../file-formats/txt-format.html\">here</a> for a list of the fields."));
  outputs.push_back(make_pair("percolator.decoy.psms.txt",
    "a tab-delimited file containing the decoy PSMs. See "
    "<a href=\"../file-formats/txt-format.html\">here</a> for a list of the fields."));
  outputs.push_back(make_pair("percolator.params.txt",
    "a file containing the name and value of all parameters for the current "
    "operation. Not all parameters in the file may have been used in the "
    "operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("percolator.pep.xml",
    "a file containing the PSMs in "
    "<a href=\"http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML\">"
    "pepXML format</a>. This file can be used as input to some of the tools in the "
    "<a href=\"http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP\">"
    "Transproteomic Pipeline</a>."));
  outputs.push_back(make_pair("percolator.mzid",
    "a file containing the protein, peptide, and spectrum matches in <a href=\""
    "http://www.psidev.info/mzidentml\">mzIdentML format</a>."));
  outputs.push_back(make_pair("percolator.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "standard error."));
  return outputs;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool PercolatorApplication::needsOutputDirectory() const {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
