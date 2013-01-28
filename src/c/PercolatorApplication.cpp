/**
 * \file PercolatorApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "PercolatorApplication.h"
#include "PercolatorAdapter.h"
#include "src/external/percolator/src/Caller.h"
#include "parameter.h"
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ios>
#include "CarpStreamBuf.h"
#include "MzIdentMLWriter.h"
#include "ProteinMatchCollection.h"


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
    "feature-file",
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
    "doc",
    "unique-peptides",
    "no-schema-validation",
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

  string input_pinxml(get_string_parameter("pin.xml"));

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
  string decoy_pre=get_string_parameter("decoy-prefix");
  if(decoy_pre.length()){
     perc_args_vec.push_back(decoy_pre);
  }else{
     perc_args_vec.push_back("random_");
  }

  perc_args_vec.push_back("-p");
  perc_args_vec.push_back(to_string(get_double_parameter("c-pos")));
 
  perc_args_vec.push_back("-n");
  perc_args_vec.push_back(to_string(get_double_parameter("c-neg")));
 
  perc_args_vec.push_back("--trainFDR");
  perc_args_vec.push_back(to_string(get_double_parameter("trian-fdr")));
 
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


  if(get_boolean_parameter("doc"))
    perc_args_vec.push_back("--doc");


  if(get_boolean_parameter("unique-peptides")&& set_protein == false){ 
    perc_args_vec.push_back("--unique-peptides");
   }
  
  if(get_boolean_parameter("no-schema-validation"))
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
  Caller* pCaller = new Caller();
  int retVal = -1;
  if (pCaller->parseOptions(perc_argc, perc_argv)) {
    retVal = pCaller->run();
  }

  // FIXME ==========================================
  //MatchCollection* res = PercolatorAdapter::psmScoresToMatchCollection(&(pCaller->fullset));
  //ProteinMatchCollection* protein_match_collection = new ProteinMatchCollection(res);
  ProteinMatchCollection* protein_match_collection = new ProteinMatchCollection();
  PercolatorAdapter::addPsmScores(protein_match_collection, &(pCaller->fullset));
  PercolatorAdapter::addPeptideScores(protein_match_collection, &(pCaller->fullset));
  PercolatorAdapter::addProteinScores(protein_match_collection, &(pCaller->fullset));
  remove("adapter_test.mzid");
  MzIdentMLWriter writer;
  writer.openFile("adapter_test.mzid", true);
  //writer.addMatches(res);
  writer.addProteinMatches(protein_match_collection);
  writer.closeFile();
  
  remove("adapter_test.txt");
  MatchFileWriter writer2("adapter_test.txt");
  vector<bool> to_print;
  to_print.assign(NUMBER_MATCH_COLUMNS, false);
  to_print[SCAN_COL] = true;
  to_print[CHARGE_COL] = true;
  to_print[SEQUENCE_COL] = true;
  to_print[PROTEIN_ID_COL] = true;
  writer2.addColumnNames(this, true, to_print);
  writer2.writeHeader();
  vector<int> zStates;
  Crux::Spectrum spectrum(1, -1, 500.0, zStates, "");  // only first scan + precursor m/z matters for this
  //res->printTabDelimited(&writer2, 100, &spectrum, PERCOLATOR_SCORE);
  
  // FIXME ==========================================


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
