/**
 *\file MakePinApplication.cpp 
 *****************************************************************************/
#include "MakePinApplication.h"
#include "PinWriter.h"
#include "parameter.h"
#include "MatchCollectionParser.h"
#include <sstream>
#include <iomanip>
#include <ios>
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
MakePinApplication::MakePinApplication() {
}

/**
 * Destructor
 */
MakePinApplication::~MakePinApplication() {}

/**
 * main method for MakePinApplication
 */
int MakePinApplication::main(int argc, char** argv) {

   /* Define optional command line arguments */

  const char* option_list[] = {
    "top-match",
    "fileroot",
    "output-dir",
    "overwrite",
    "output-file",
    "verbosity",
    "parameter-file",
    "list-of-files",
    "filestem-prefixes"
  };

  int num_options = sizeof(option_list) / sizeof(char*);
  
  /* Define required command line arguments */
  const char* argument_list[] = {
    "target input", 
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, 
    num_arguments,
    option_list, 
    num_options, 
    argc, 
    argv)
  ;

  string target_path = string(get_string_parameter_pointer("target input"));

  vector<string> search_result_files;
  get_search_result_paths(target_path, search_result_files);

  return main(search_result_files);
}

/**
 * \runs make-pin application
 */
int MakePinApplication::main(vector<string>& paths) {
  //create MatchColletion 
  MatchCollectionParser parser;

  if (paths.size() == 0) {
    carp(CARP_FATAL, "No search paths found!");
  }

  MatchCollection* target_collection = new MatchCollection();
  MatchCollection* decoy_collection = new MatchCollection();

  for (vector<string>::iterator iter = paths.begin(); iter != paths.end(); ++iter) {
    carp(CARP_INFO, "Parsing %s", iter->c_str());
    MatchCollection* current_collection = parser.create(iter->c_str(), "__NULL_STR");
    for (int scorer_idx = (int)SP; scorer_idx < (int)NUMBER_SCORER_TYPES; scorer_idx++) {
      target_collection->setScoredType((SCORER_TYPE_T)scorer_idx, 
        current_collection->getScoredType((SCORER_TYPE_T)scorer_idx));
      decoy_collection->setScoredType((SCORER_TYPE_T)scorer_idx,
        current_collection->getScoredType((SCORER_TYPE_T)scorer_idx));
    } 
    MatchIterator* match_iter = new MatchIterator(current_collection);
    while(match_iter->hasNext()) {
      Crux::Match* match = match_iter->next();
      if (match->getNullPeptide()) {
        decoy_collection->addMatch(match);
      } else {
        target_collection->addMatch(match);
      }
    }
    delete match_iter;
    delete current_collection;
  }

  carp(CARP_INFO, "There are %d target matches and %d decoys",target_collection->getMatchTotal(), decoy_collection->getMatchTotal());
  if (target_collection->getMatchTotal() == 0) {
    carp(CARP_FATAL, "No target matches found!");
  }
  if (decoy_collection->getMatchTotal() == 0) {
    carp(CARP_FATAL, "No decoy matches found!  Did you set 'decoy-prefix' properly?");
  }

  PinWriter* writer = new PinWriter();

  string output_dir=get_string_parameter_pointer("output-dir");
 
  //perpare output file 
  string output_filename = get_string_parameter_pointer("output-file");
  if(output_filename=="" || output_filename=="__NULL_STR"){
    string fileroot = get_string_parameter_pointer("fileroot");
    if (fileroot == "" || fileroot == "__NULL_STR")
      fileroot = "";
    else
      fileroot += ".";
    output_filename = fileroot + "make-pin.pin";
  }
  writer->openFile(output_filename.c_str(),output_dir.c_str(), get_boolean_parameter("overwrite"));

 //write .pin file 
  vector<MatchCollection*> decoys;
  decoys.push_back(decoy_collection);
  writer->write(target_collection, decoys, get_int_parameter("top-match"));

  //close file 
  writer->closeFile();

  delete target_collection;
  delete decoy_collection;
  delete writer;

  return 0;
}


/**
 * \returns the command name for PercolatorApplication
 */
string MakePinApplication::getName() {
  return "make-pin";
}

/**
 * \returns the description for PercolatorApplication
 */
string MakePinApplication::getDescription() {
  return "Given a set of search results files, generate a pin file for input "
         "to crux percolator";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool MakePinApplication::needsOutputDirectory() {
  return true;
}


bool MakePinApplication:: hidden(){
  return false;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
