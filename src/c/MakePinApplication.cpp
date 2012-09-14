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
    "fileroot",
    "output-dir",
    "overwrite",
    "output-file",
    "verbosity",
    "parameter-file"

  };

  int num_options = sizeof(option_list) / sizeof(char*);
  
  /* Define required command line arguments */
  const char* argument_list[] = {"target input", 
    "decoy input",
    "protein database"
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
  char* target_path = get_string_parameter("target input");
  //Get input : decoy 
  char* decoy_path = get_string_parameter("decoy input");
  //Get input : protein 
  char* protein_dbase = get_string_parameter("protein database");
  //prepare output file

  //create MatchColletion 
  vector<MatchCollection*> decoys;
  MatchCollection* target_collection = MatchCollectionParser::create(target_path, protein_dbase); 
  MatchCollection* decoy_collection =  MatchCollectionParser::create(decoy_path, protein_dbase);
  PinXMLWriter* writer=new PinXMLWriter();
  decoys.push_back(decoy_collection);

  bool overwrite = false; 
  if(get_boolean_parameter("overwrite"))
    overwrite=true; 

  string output_dir=get_string_parameter("output-dir");
 
  //perpare output file 
  char* output_filename = get_string_parameter("output-file");
  if(output_filename==NULL){
    output_filename ="make-pin.pin.xml";
  }
  writer->openFile(output_filename,output_dir.c_str(),overwrite);

  //set process information 
  if(target_path!=NULL && decoy_path!=NULL )
    writer->setProcessInfo(target_path, decoy_path);
 //write .pin.xml file 
  writer->write(target_collection, decoys);
  writer->printFooter();

  //close file 
  writer->closeFile();

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
