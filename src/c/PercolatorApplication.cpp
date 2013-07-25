/**
 * \file PercolatorApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "MakePinApplication.h"
#include "PercolatorApplication.h"
#include "PercolatorAdapter.h"
#include "Caller.h"
#include "parameter.h"
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ios>
#include "CarpStreamBuf.h"
#include "MzIdentMLWriter.h"
#include "ProteinMatchCollection.h"
#include "PMCDelimitedFileWriter.h"
#include "PMCPepXMLWriter.h"
#include "PMCSQTWriter.h"


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

   /* Define optional command line arguments */

  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "txt-output",
    "mzid-output",
    "pepxml-output",
    "feature-file",
    "feature-in-file",
    "parameter-file",
    "protein",
    "decoy-xml-output",
    "decoy-prefix",
    "c-pos",
    "c-neg",
    "train-fdr",
    "test-fdr",
    "maxiter",
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
    "seed",
    "klammer",
    //"doc",
    "only-psms",
    "allow-protein-group",
    "protein-level-pi0",
    "empirical-protein-q",
    "group-proteins",
    "no-prune-proteins",
    "deepness",
    "verbosity"
  };

  int num_options = sizeof(option_list) / sizeof(char*);
  
  /* Define required command line arguments */
  const char* argument_list[] = {"pin.xml"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, 
    num_arguments,
    option_list, 
    num_options, 
    argc, 
    argv)
  ;

  string input_pinxml(get_string_parameter_pointer("pin.xml"));
  
  if (!get_boolean_parameter("feature-in-file")) {
    // Check if we need to run make-pin first
    if (has_extension(input_pinxml.c_str(), "txt") ||
        has_extension(input_pinxml.c_str(), "sqt") ||
        has_extension(input_pinxml.c_str(), "pep.xml")) {
      string input_decoy(input_pinxml);
      int target_pos = input_pinxml.find("target");
      if (target_pos < 0) {
        int decoy_pos = input_pinxml.find("decoy");
        if (decoy_pos < 0) {
          carp(CARP_FATAL, "Not a PIN, target, or decoy file: %s",
                           input_pinxml.c_str());
        }
        // user gave decoy results file
        input_pinxml.replace(decoy_pos, 5, "target");
      } else {
        // user gave target results file
        input_decoy.replace(target_pos, 6, "decoy");
      }
  
      // check that files exist
      if (!file_exists(input_pinxml)) {
        carp(CARP_FATAL, "Target file %s not found", input_pinxml.c_str());
      } else if (!file_exists(input_decoy)) {
        carp(CARP_FATAL, "Decoy file %s not found", input_decoy.c_str());
      }
  
      string pin_location = string(get_string_parameter_pointer("output-dir")) +
                            "/make-pin.pin.xml";
      const char* make_pin_file = pin_location.c_str();
  
      carp(CARP_INFO, "Running make-pin with '%s' and decoy file '%s'.",
           input_pinxml.c_str(), input_decoy.c_str());
      int ret = MakePinApplication::main(input_pinxml, input_decoy);
      if (ret != 0 || !file_exists(make_pin_file)) {
        carp(CARP_FATAL, "make-pin failed. Not running Percolator.");
      }
      carp(CARP_INFO, "Finished make-pin.");
      input_pinxml = string(make_pin_file);
    }
  }
  return main(input_pinxml);

}

/**
 * \brief runs percolator on the input pin.xml
 * \returns whether percolator was successful or not
 */
int PercolatorApplication::main(
  const string& input_pinxml ///< file path of pin.xml to process.
  ) {

  /* Write the data files */

  string output_xml = make_file_path("percolator.target.pout.xml");
  string output_target_tab = make_file_path("percolator.target.txt");
  string output_decoy_tab = make_file_path("percolator.decoy.txt");



  /* build argument list */
  vector<string> perc_args_vec;
  perc_args_vec.push_back("percolator");


  perc_args_vec.push_back("-X");
  perc_args_vec.push_back(output_xml);

  // These files will be removed later
  // They are only set so that the tab-delimited info is not written to stdout
  perc_args_vec.push_back("-r");
  perc_args_vec.push_back(output_target_tab);
  perc_args_vec.push_back("-B");
  perc_args_vec.push_back(output_decoy_tab);


  //Add options

  //TODO remove this dependency.
  //check that the -X is set or  not     
    
  bool set_protein = get_boolean_parameter("protein");

  if(get_boolean_parameter("protein")){
     perc_args_vec.push_back("-A");
  }
  
  if(get_boolean_parameter("decoy-xml-output")){
    perc_args_vec.push_back("-Z");
  }

  perc_args_vec.push_back("-P");
  string decoy_pre=get_string_parameter_pointer("decoy-prefix");
  if(decoy_pre.length()){
     perc_args_vec.push_back(decoy_pre);
  }else{
     perc_args_vec.push_back("random_");
  }

  string seed_parameter = get_string_parameter_pointer("seed");
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

  perc_args_vec.push_back("-p");
  perc_args_vec.push_back(to_string(get_double_parameter("c-pos")));
 
  perc_args_vec.push_back("-n");
  perc_args_vec.push_back(to_string(get_double_parameter("c-neg")));
 
  perc_args_vec.push_back("--trainFDR");
  perc_args_vec.push_back(to_string(get_double_parameter("train-fdr")));
 
  perc_args_vec.push_back("--testFDR");
  perc_args_vec.push_back(to_string(get_double_parameter("test-fdr")));

  perc_args_vec.push_back("--maxiter");
  perc_args_vec.push_back(to_string(get_int_parameter("maxiter")));

  perc_args_vec.push_back("--train-ratio");
  perc_args_vec.push_back(to_string(get_double_parameter("train-ratio")));

  if(get_boolean_parameter("feature-file")){ 
    perc_args_vec.push_back("--tab-out");
    string feature_output=make_file_path("percolator.feature.txt");
    perc_args_vec.push_back(feature_output);
  }

  if(get_boolean_parameter("output-weights")){
    perc_args_vec.push_back("--weights");
    string output_wght=make_file_path("percolator.weights.txt");
    perc_args_vec.push_back(output_wght);
  }
  
  if(string(get_string_parameter_pointer("input-weights"))!="__NULL_STR"){
    perc_args_vec.push_back("--init-weights");
    perc_args_vec.push_back(get_string_parameter_pointer("input-weights"));
  }

  perc_args_vec.push_back("--default-direction");
  perc_args_vec.push_back(to_string(get_int_parameter("default-direction")));


  if(get_boolean_parameter("unitnorm"))
    perc_args_vec.push_back("-u");
  

  if(set_protein){
    perc_args_vec.push_back("--alpha");
    perc_args_vec.push_back(to_string(get_double_parameter("alpha")));
 
    perc_args_vec.push_back("--beta");
    perc_args_vec.push_back(to_string(get_double_parameter("beta")));

    perc_args_vec.push_back("--gamma");
    perc_args_vec.push_back(to_string(get_double_parameter("gamma")));
  }
  
  
  if(get_boolean_parameter("test-each-iteration"))
    perc_args_vec.push_back("--test-each-iteration");

  if(get_boolean_parameter("static-override")){
    perc_args_vec.push_back("--override");
  }

 
  if(get_boolean_parameter("klammer"))
      perc_args_vec.push_back("--klammer");


  /* --doc option disabled, need retention times in pinxml file
  int doc_parameter = get_int_parameter("doc");
  if(doc_parameter >= 0) {
    perc_args_vec.push_back("--doc");
    perc_args_vec.push_back(to_string(doc_parameter));
  }
  */

  if (get_boolean_parameter("only-psms") && !set_protein) { 
    perc_args_vec.push_back("--unique-peptides");
  }
  
  // FIXME include schema as part of distribution and add option to turn on validation
    perc_args_vec.push_back("-s");
  
  if(get_boolean_parameter("allow-protein-group"))
    perc_args_vec.push_back("allow-protein-group");
 
  if(set_protein){ 
    if(get_boolean_parameter("protein-level-pi0"))
      perc_args_vec.push_back("-I");
   
    if(get_boolean_parameter("empirical-protein-q"))
       perc_args_vec.push_back("--empirical-protein-q");

    if(get_boolean_parameter("group-proteins"))
       perc_args_vec.push_back("--group-proteins");
    
    if(get_boolean_parameter("no-prune-proteins"))
      perc_args_vec.push_back("--no-prune-proteins"); 
   
    perc_args_vec.push_back("--deepness");
    perc_args_vec.push_back(to_string(get_int_parameter("deepness")));
  }
  
  //This has to be a the end in order for it to work.
  if(get_boolean_parameter("feature-in-file")) {
    perc_args_vec.push_back("-j");
  }

  
   perc_args_vec.push_back(input_pinxml);

   
  /* build argv line */
  int perc_argc = perc_args_vec.size();

  char** perc_argv = new char*[perc_argc];

  perc_argv[0] = (char*)perc_args_vec[0].c_str();
  //cerr<<"perc_argv["<<0<<"]= "<<perc_argv[0]<<endl;
  for (int idx = 1;idx < perc_argc ; idx++) {
    perc_argv[idx] = (char*)perc_args_vec[idx].c_str();
    carp(CARP_DEBUG, "perc_argv[%d]=%s", idx, perc_argv[idx]);
    //cerr<<"perc_argv["<<idx<<"]= "<<perc_argv[idx]<<endl;
  }

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

  // remove tab files
  remove(output_target_tab.c_str());
  remove(output_decoy_tab.c_str());

  // get percolator score information into crux objects
  
  ProteinMatchCollection* protein_match_collection = pCaller->getProteinMatchCollection();
  ProteinMatchCollection* decoy_protein_match_collection =
    pCaller->getDecoyProteinMatchCollection();
  string output_dir = string(get_string_parameter_pointer("output-dir"));

  // write txt
  if (get_boolean_parameter("txt-output")) {
    PMCDelimitedFileWriter txt_writer;
    string txt_path = make_file_path("percolator.target");
    string decoy_txt_path = make_file_path("percolator.decoy");

    txt_writer.writeFile(this, txt_path + ".psms.txt",
                         PMCDelimitedFileWriter::PSMS, protein_match_collection);
    txt_writer.writeFile(this, decoy_txt_path + ".psms.txt",
                         PMCDelimitedFileWriter::PSMS, decoy_protein_match_collection);
    txt_writer.writeFile(this, txt_path + ".peptides.txt",
                         PMCDelimitedFileWriter::PEPTIDES, protein_match_collection);
    txt_writer.writeFile(this, decoy_txt_path + ".peptides.txt",
                         PMCDelimitedFileWriter::PEPTIDES, decoy_protein_match_collection);

    if (set_protein) {
      txt_writer.writeFile(this, txt_path + ".proteins.txt",
                           PMCDelimitedFileWriter::PROTEINS, protein_match_collection);
      txt_writer.writeFile(this, decoy_txt_path + ".proteins.txt",
                           PMCDelimitedFileWriter::PROTEINS, decoy_protein_match_collection);
    }
  }

  // write mzid
  if (get_boolean_parameter("mzid-output")) {
    MzIdentMLWriter mzid_writer, decoy_mzid_writer;
    string mzid_path = make_file_path("percolator.target.mzid");
    mzid_writer.openFile(mzid_path, get_boolean_parameter("overwrite"));
    mzid_writer.addProteinMatches(protein_match_collection);
    mzid_writer.closeFile();
    mzid_path = make_file_path("percolator.decoy.mzid");
    decoy_mzid_writer.openFile(mzid_path, get_boolean_parameter("overwrite"));
    decoy_mzid_writer.addProteinMatches(decoy_protein_match_collection);
    decoy_mzid_writer.closeFile();
  }
  
  // write pepxml
  if (get_boolean_parameter("pepxml-output")) {
    PMCPepXMLWriter pep_writer;
    string pep_path = make_file_path("percolator.target.pep.xml");
    pep_writer.openFile(pep_path.c_str(), get_boolean_parameter("overwrite"));
    pep_writer.write(protein_match_collection);
    pep_writer.closeFile();
    pep_path = make_file_path("percolator.decoy.pep.xml");
    pep_writer.openFile(pep_path.c_str(), get_boolean_parameter("overwrite"));
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

COMMAND_T PercolatorApplication::getCommand() {
  return PERCOLATOR_COMMAND;

}

/**
 * \returns the command name for PercolatorApplication
 */
string PercolatorApplication::getName() {
  return "percolator";
}

/**
 * \returns the description for PercolatorApplication
 */
string PercolatorApplication::getDescription() {

  return "Runs percolator";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool PercolatorApplication::needsOutputDirectory() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
