/**
 *\file MakePinApplication.cpp 
 *****************************************************************************/
#include "MakePinApplication.h"
#include "PinXMLWriter.h"
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
    "parameter-file"

  };

  int num_options = sizeof(option_list) / sizeof(char*);
  
  /* Define required command line arguments */
  const char* argument_list[] = {
    "target input", 
    "decoy input",
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

  //Get input: target
  string target_path = string(get_string_parameter_pointer("target input"));
  //Get input : decoy 
  string decoy_path = string(get_string_parameter_pointer("decoy input"));

  return main(target_path, decoy_path);
}

/**
 * \runs make-pin application
 */
int MakePinApplication::main(string target_path, string decoy_path) {
  //create MatchColletion 
  vector<MatchCollection*> decoys;
  MatchCollection* target_collection =
    MatchCollectionParser::create(target_path.c_str(), "__NULL_STR"); 
  MatchCollection* decoy_collection =
    MatchCollectionParser::create(decoy_path.c_str(), "__NULL_STR");

  // Mark decoy matches
  MatchIterator* decoy_iter = new MatchIterator(decoy_collection);
  while (decoy_iter->hasNext()) {
    Crux::Match* decoy_match = decoy_iter->next();
    decoy_match->setNullPeptide(true);
  }
  delete decoy_iter;

  PinXMLWriter* writer=new PinXMLWriter();
  decoys.push_back(decoy_collection);

  bool overwrite = false; 
  if(get_boolean_parameter("overwrite"))
    overwrite=true; 

  string output_dir=get_string_parameter("output-dir");
 
  //perpare output file 
  string output_filename = get_string_parameter_pointer("output-file");
  if(output_filename=="" || output_filename=="__NULL_STR"){
    string fileroot = get_string_parameter_pointer("fileroot");
    if (fileroot == "" || fileroot == "__NULL_STR")
      fileroot = "";
    else
      fileroot += ".";
    output_filename = fileroot + "make-pin.pin.xml";
  }
  writer->openFile(output_filename.c_str(),output_dir.c_str(),overwrite);

  //set process information 
  writer->setProcessInfo(target_path.c_str(), decoy_path.c_str());
 //write .pin.xml file 
  writer->write(target_collection, decoys, get_int_parameter("top-match"));
  writer->printFooter();

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

  return "Runs make-pin";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool MakePinApplication::needsOutputDirectory() {
  return true;
}


/**
 * hide sequest search 
*/

bool MakePinApplication:: hidden(){
  return true; 
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
