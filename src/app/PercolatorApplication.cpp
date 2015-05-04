/**
 * \file PercolatorApplication.cpp 
 * \brief Runs Percolator
 *****************************************************************************/
#include "MakePinApplication.h"
#include "PercolatorApplication.h"
#include "PercolatorAdapter.h"
#include "Caller.h"
#include "parameter.h"
#include "util/Params.h"
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ios>
#include "util/CarpStreamBuf.h"
#include "util/FileUtils.h"
#include "io/MzIdentMLWriter.h"
#include "model/ProteinMatchCollection.h"
#include "io/PMCDelimitedFileWriter.h"
#include "io/PMCPepXMLWriter.h"
#include "io/PMCSQTWriter.h"


using namespace std;
  /**
   * Turn the given value into a string.  For floating point numbers,
   * use the --precision option value.
   */
  template<typename TValue>
  static string to_string(TValue value) {

    ostringstream oss;
    oss << setprecision(get_int_parameter("precision")) << fixed;
    oss << value;
    string out_string = oss.str();
    return out_string;
  }

  /**
   * Turn the given value into a string.  For floating point numbers,
   * use the given precision.
   */
  template<typename TValue>
  static string to_string
    (TValue& value,
     int precision,
     bool fixed_float = true) {

    ostringstream oss;
    oss << setprecision(precision);
    if (fixed_float) {
      oss << fixed;
    } else {
      oss.unsetf(ios_base::floatfield);
    }
    oss << value;
    string out_string = oss.str();
    return out_string;
  }


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
  initialize(argc, argv);

  string input_pin = Params::GetString("pin");

  // Check if we need to run make-pin first
  if (Params::GetBool("list-of-files") || 
      has_extension(input_pin.c_str(), "txt") ||
      has_extension(input_pin.c_str(), "sqt") ||
      has_extension(input_pin.c_str(), "pep.xml") ||
      has_extension(input_pin.c_str(), "mzid")) {

    vector<string> result_files;
    get_search_result_paths(input_pin, result_files);

    string pin_location = make_file_path("make-pin.pin");

    const char* make_pin_file = pin_location.c_str();

    carp(CARP_INFO, "Running make-pin");
    int ret = MakePinApplication::main(result_files);

    if (ret != 0 || !FileUtils::Exists(make_pin_file)) {
      carp(CARP_FATAL, "make-pin failed. Not running Percolator.");
    }
    carp(CARP_INFO, "Finished make-pin.");
    input_pin = string(make_pin_file);
  } else if (!has_extension(input_pin.c_str(), "pin")) {
      carp(CARP_FATAL, "input file %s is not recognized.", input_pin.c_str() );
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

  if(Params::GetBool("decoy-xml-output")){
    perc_args_vec.push_back("-Z");
  }

  perc_args_vec.push_back("-P");
  string decoy_pre = Params::GetString("decoy-prefix");
  if(decoy_pre.length()){
     perc_args_vec.push_back(decoy_pre);
  } else {
     perc_args_vec.push_back("random_");
  }

  string seed_parameter = Params::GetString("percolator-seed");
  unsigned int seed_value;
  if (seed_parameter == "time") {
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    // percolator accepts seed values 1-20000
    seed_value = (unsigned int)seconds % 20000 + 1;
  } else {
    // seed 0 causes segfault in percolator 
    stringstream seed_extractor(seed_parameter);
    seed_extractor >> seed_value;
    if (seed_value == 0) {
      ++seed_value;
    }
  }
  stringstream seed_stream;
  seed_stream << seed_value;
  perc_args_vec.push_back("--seed");
  perc_args_vec.push_back(seed_stream.str());

  if (Params::GetBool("pout-output")) {
    perc_args_vec.push_back("-X");
    perc_args_vec.push_back(make_file_path(getFileStem() + ".pout.xml"));
  }

  perc_args_vec.push_back("-p");
  perc_args_vec.push_back(to_string(Params::GetDouble("c-pos")));
 
  perc_args_vec.push_back("-n");
  perc_args_vec.push_back(to_string(Params::GetDouble("c-neg")));
 
  perc_args_vec.push_back("--trainFDR");
  perc_args_vec.push_back(to_string(Params::GetDouble("train-fdr")));
 
  perc_args_vec.push_back("--testFDR");
  perc_args_vec.push_back(to_string(Params::GetDouble("test-fdr")));

  perc_args_vec.push_back("--maxiter");
  perc_args_vec.push_back(to_string(Params::GetInt("maxiter")));

  if (Params::GetBool("quick-validation")) {
    perc_args_vec.push_back("--quick-validation");
  }

  perc_args_vec.push_back("--train-ratio");
  perc_args_vec.push_back(to_string(Params::GetDouble("train-ratio")));

  if(Params::GetBool("feature-file")){ 
    perc_args_vec.push_back("--tab-out");
    string feature_output=make_file_path(getFileStem() + ".feature.txt");
    perc_args_vec.push_back(feature_output);
  }

  if(Params::GetBool("output-weights")){
    perc_args_vec.push_back("--weights");
    string output_wght=make_file_path(getFileStem() + ".weights.txt");
    perc_args_vec.push_back(output_wght);
  }
  
  if (!Params::GetString("input-weights").empty()) {
    perc_args_vec.push_back("--init-weights");
    perc_args_vec.push_back(Params::GetString("input-weights"));
  }

  if (!Params::GetString("default-direction").empty()) {  
    perc_args_vec.push_back("--default-direction");
    perc_args_vec.push_back(Params::GetString("default-direction"));
  }

  if(Params::GetBool("unitnorm"))
    perc_args_vec.push_back("-u");

  if(Params::GetBool("test-each-iteration"))
    perc_args_vec.push_back("--test-each-iteration");

  if(Params::GetBool("static-override")){
    perc_args_vec.push_back("--override");
  }
 
  if(Params::GetBool("klammer"))
      perc_args_vec.push_back("--klammer");

  /* --doc option disabled, need retention times in pin file
  int doc_parameter = get_int_parameter("doc");
  if(doc_parameter >= 0) {
    perc_args_vec.push_back("--doc");
    perc_args_vec.push_back(to_string(doc_parameter));
  }
  */

  // FIXME include schema as part of distribution and add option to turn on validation
  perc_args_vec.push_back("-s");

  if(Params::GetBool("allow-protein-group"))
    perc_args_vec.push_back("--allow-protein-group");
 
  bool set_protein = Params::GetBool("protein");
  if(set_protein){ 
    perc_args_vec.push_back("-A");

    if (Params::GetDouble("alpha") > 0) {
      perc_args_vec.push_back("--fido-alpha");
      perc_args_vec.push_back(to_string(Params::GetDouble("alpha")));
    }
    if (Params::GetDouble("beta") > 0) {
      perc_args_vec.push_back("--fido-beta");
      perc_args_vec.push_back(to_string(Params::GetDouble("beta")));
    }
    if (Params::GetDouble("gamma") > 0) {
      perc_args_vec.push_back("--fido-gamma");
      perc_args_vec.push_back(to_string(Params::GetDouble("gamma")));
    }

    if(Params::GetBool("protein-level-pi0"))
      perc_args_vec.push_back("-I");
   
    if(Params::GetBool("empirical-protein-q"))
       perc_args_vec.push_back("--empirical-protein-q");

    if(!Params::GetBool("group-proteins"))
       perc_args_vec.push_back("--fido-no-group-proteins");

    if (Params::GetBool("no-separate-proteins")) {
      perc_args_vec.push_back("--fido-no-separate-proteins");
    }
    
    if(Params::GetBool("no-prune-proteins"))
      perc_args_vec.push_back("--fido-no-prune-proteins"); 

    perc_args_vec.push_back("--fido-gridsearch-depth");
    perc_args_vec.push_back(to_string(Params::GetInt("deepness")));

    if (Params::GetBool("reduce-tree-in-gridsearch")) {
      perc_args_vec.push_back("--fido-reduce-tree-in-gridsearch");
    }

    if (Params::GetBool("post-processing-tdcn")) {
      perc_args_vec.push_back("--post-processing-tdcn");
    }

    perc_args_vec.push_back("--grid-search-mse-threshold");
    perc_args_vec.push_back(Params::GetString("grid-search-mse-threshold"));

    if (Params::GetBool("truncation")) {
      perc_args_vec.push_back("--fido-truncation");
    }

    if (Params::GetBool("protein-group-level-inference")) {
      perc_args_vec.push_back("--fido-protein-group-level-inference");
    }

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
  int perc_argc = perc_args_vec.size();

  string perc_cmd = perc_args_vec[0];
  
  char** perc_argv = new char*[perc_argc];

  perc_argv[0] = (char*)perc_args_vec[0].c_str();
  //cerr<<"perc_argv["<<0<<"]= "<<perc_argv[0]<<endl;
  for (int idx = 1;idx < perc_argc ; idx++) {
    perc_argv[idx] = (char*)perc_args_vec[idx].c_str();
    perc_cmd = perc_cmd + " " + perc_argv[idx];
    carp(CARP_DEBUG, "perc_argv[%d]=%s", idx, perc_argv[idx]);
    //cerr<<"perc_argv["<<idx<<"]= "<<perc_argv[idx]<<endl;
  }

  carp(CARP_DEBUG, "cmd:%s", perc_cmd.c_str());
  
  /* Re-route stdeer to log file. */
  CarpStreamBuf buffer;
  streambuf* old = std::cerr.rdbuf();
  std::cerr.rdbuf(&buffer);

  /* Call percolatorMain */
  PercolatorAdapter* pCaller = new PercolatorAdapter();
  int retVal = -1;
  if (pCaller->parseOptions(perc_argc, perc_argv)) {
    retVal = pCaller->run();
  } 
  
  carp(CARP_DEBUG, "Percolator retval:%d", retVal);
  if (retVal != 0) {
    carp(CARP_FATAL, "Error running percolator:%d", retVal);
  }

  // get percolator score information into crux objects
  ProteinMatchCollection* protein_match_collection =
    pCaller->getProteinMatchCollection();
  ProteinMatchCollection* decoy_protein_match_collection =
    pCaller->getDecoyProteinMatchCollection();
  string output_dir = Params::GetString("output-dir");

  // write txt
  if (!Params::GetBool("original-output")) {
    FileUtils::Remove(output_target_peptides);
    FileUtils::Remove(output_target_proteins);
  }
  if (Params::GetBool("txt-output") && !Params::GetBool("original-output")) {
    PMCDelimitedFileWriter txt_writer;
    txt_writer.writeFile(this, output_target_psms,
                         PMCDelimitedFileWriter::PSMS, protein_match_collection);
    txt_writer.writeFile(this, output_decoy_psms,
                         PMCDelimitedFileWriter::PSMS, decoy_protein_match_collection);
    txt_writer.writeFile(this, output_target_peptides,
                         PMCDelimitedFileWriter::PEPTIDES, protein_match_collection);
    txt_writer.writeFile(this, output_decoy_peptides,
                         PMCDelimitedFileWriter::PEPTIDES, decoy_protein_match_collection);

    if (set_protein) {
      txt_writer.writeFile(this, output_target_proteins,
                           PMCDelimitedFileWriter::PROTEINS, protein_match_collection);
      txt_writer.writeFile(this, output_decoy_proteins,
                           PMCDelimitedFileWriter::PROTEINS, decoy_protein_match_collection);
    }
  }

  // write mzid
  if (Params::GetBool("mzid-output")) {
    MzIdentMLWriter mzid_writer, decoy_mzid_writer;
    string mzid_path = make_file_path(getFileStem() + ".target.mzid");
    mzid_writer.openFile(mzid_path, Params::GetBool("overwrite"));
    mzid_writer.addProteinMatches(protein_match_collection);
    mzid_writer.closeFile();
    mzid_path = make_file_path(getFileStem() + ".decoy.mzid");
    decoy_mzid_writer.openFile(mzid_path, Params::GetBool("overwrite"));
    decoy_mzid_writer.addProteinMatches(decoy_protein_match_collection);
    decoy_mzid_writer.closeFile();
  }
  
  // write pepxml
  if (Params::GetBool("pepxml-output")) {
    PMCPepXMLWriter pep_writer;
    string pep_path = make_file_path(getFileStem() + ".target.pep.xml");
    pep_writer.openFile(pep_path.c_str(), Params::GetBool("overwrite"));
    pep_writer.write(protein_match_collection);
    pep_writer.closeFile();
    pep_path = make_file_path(getFileStem() + ".decoy.pep.xml");
    pep_writer.openFile(pep_path.c_str(), Params::GetBool("overwrite"));
    pep_writer.write(decoy_protein_match_collection);
    pep_writer.closeFile();
  }

  delete protein_match_collection;
  
  delete pCaller;
  Globals::clean();
  delete []perc_argv;

  /* Recover stderr */
  std::cerr.rdbuf( old );



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
    "Percolator to represent each PSM are summarized <a href=\"features.html\">"
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
    "respects:</p><ul><li>In addition to the native Percolator tab-delimited "
    "file format (.pin), Crux Percolator supports additional input file "
    "formats (SQT, PepXML, Crux tab-delimited text) and output file formats "
    "(PepXML, mzIdentML, Crux tab-delimited text).</li>"
    "<li>To maintain consistency with the rest of the "
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
    "feature-file",
    "list-of-files",
    "parameter-file",
    "protein",
    "decoy-xml-output",
    "decoy-prefix",
    "c-pos",
    "c-neg",
    "train-fdr",
    "test-fdr",
    "maxiter",
    "quick-validation",
    "train-ratio",
    "output-weights",
    "input-weights",
    "default-direction",
    "unitnorm",
    "alpha",
    "beta",
    "gamma",
    "test-each-iteration",
    "static-override",
    "percolator-seed",
    "klammer",
    "only-psms",
    //"doc",
    "allow-protein-group",
    "protein-level-pi0",
    "empirical-protein-q",
    "group-proteins",
    "no-separate-proteins",
    "no-prune-proteins",
    "deepness",
    "reduce-tree-in-gridsearch",
    "post-processing-tdcn",
    "grid-search-mse-threshold",
    "truncation",
    "protein-group-level-inference",
    "verbosity",
    "top-match"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
map<string, string> PercolatorApplication::getOutputs() const {
  map<string, string> outputs;
  outputs["percolator.target.proteins.txt"] =
    "a tab-delimited file containing the target protein matches. See "
    "<a href=\"txt-format.html\">here</a> for a list of the fields.";
  outputs["percolator.decoy.proteins.txt"] =
    "a tab-delimited file containing the decoy protein matches. See "
    "<a href=\"txt-format.html\">here</a> for a list of the fields.";
  outputs["percolator.target.peptides.txt"] =
    "a tab-delimited file containing the target peptide matches. See "
    "<a href=\"txt-format.html\">here</a> for a list of the fields.";
  outputs["percolator.decoy.peptides.txt"] =
    "a tab-delimited file containing the decoy peptide matches. See "
    "<a href=\"txt-format.html\">here</a> for a list of the fields.";
  outputs["percolator.target.psms.txt"] =
    "a tab-delimited file containing the target PSMs. See "
    "<a href=\"txt-format.html\">here</a> for a list of the fields.";
  outputs["percolator.decoy.psms.txt"] =
    "a tab-delimited file containing the decoy PSMs. See "
    "<a href=\"txt-format.html\">here</a> for a list of the fields.";
  outputs["percolator.params.txt"] =
    "a file containing the name and value of all parameters for the current "
    "operation. Not all parameters in the file may have been used in the "
    "operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs.";
  outputs["percolator.pep.xml"] =
    "a file containing the PSMs in "
    "<a href=\"http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML\">"
    "pepXML format</a>. This file can be used as input to some of the tools in the "
    "<a href=\"http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP\">"
    "Transproteomic Pipeline</a>.";
  outputs["percolator.mzid"] =
    "a file containing the protein, peptide, and spectrum matches in <a href=\""
    "http://www.psidev.info/mzidentml\">mzIdentML format</a>.";
  outputs["percolator.log.txt"] =
    "a log file containing a copy of all messages that were printed to "
    "standard error.";
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
