#include <cstdio>

#include "io/carp.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "PSMConvertApplication.h"
#include "model/MatchCollection.h"
#include "model/ProteinMatchCollection.h"
#include "io/HTMLWriter.h"
#include "io/MatchFileReader.h"
#include "io/MzIdentMLReader.h"
#include "io/MzIdentMLWriter.h"
#include "io/PepXMLReader.h"
#include "io/PepXMLWriter.h"
#include "io/PinWriter.h"
#include "io/PMCDelimitedFileWriter.h"
#include "io/PMCPepXMLWriter.h"
#include "io/PMCSQTWriter.h"
#include "io/PSMReader.h"
#include "io/PSMWriter.h"
#include "io/SQTReader.h"

PSMConvertApplication::PSMConvertApplication() {
}

PSMConvertApplication::~PSMConvertApplication() {
}

void PSMConvertApplication::convertFile(string input_format, string output_format, string input_file, string output_file_base, string database_file, bool distinct_matches) {
  Database* data;
  if (database_file.empty()) {
    data = new Database();
    carp(CARP_INFO, "Database not provided; will use empty database.");
  } else {
    data = new Database(database_file.c_str(), false);
    carp(CARP_INFO, "Created database using %s.", database_file.c_str());
  }
  
  bool isTabDelimited = false;
  PSMReader* reader;
  
  if (input_format != "auto") {
    if (input_format == "tsv") {
      reader = new MatchFileReader(input_file.c_str(), data);
      isTabDelimited = true;
    } else if (input_format == "html") {
      carp(CARP_FATAL, "HTML format has not been implemented yet");
    } else if (input_format == "sqt") {
      reader = new SQTReader(input_file.c_str(), data);
    } else if (input_format == "pin") {
      carp(CARP_FATAL, "Pin format has not been implemented yet");
    } else if (input_format == "pepxml") {
      reader = new PepXMLReader(input_file.c_str(), data);
    } else if (input_format == "mzidentml") {
      reader = new MzIdentMLReader(input_file.c_str(), data);
    } else if (input_format == "barista-xml") {
      carp(CARP_FATAL, "Barista-XML format has not been implemented yet");
    } else {
      carp(CARP_FATAL, "Invalid input format. Valid formats are: tsv, html, "
           "sqt, pin, pepxml, mzidentml, barista-xml.");
    }
  } else {
    if (StringUtils::IEndsWith(input_file, ".txt")) {
      reader = new MatchFileReader(input_file.c_str(), data);
      isTabDelimited = true;
    } else if (StringUtils::IEndsWith(input_file, ".html")) {
      carp(CARP_FATAL, "HTML format has not been implemented yet");
    } else if (StringUtils::IEndsWith(input_file, ".sqt")) {
      reader = new SQTReader(input_file.c_str(), data);
    } else if (StringUtils::IEndsWith(input_file, ".pin")) {
      carp(CARP_FATAL, "Pin format has not been implemented yet");
    } else if (StringUtils::IEndsWith(input_file, ".xml")) {
      reader = new PepXMLReader(input_file.c_str(), data);
    } else if (StringUtils::IEndsWith(input_file, ".mzid")) {
      reader = new MzIdentMLReader(input_file.c_str(), data);
    } else if (StringUtils::IEndsWith(input_file, ".barista.xml")) {
      carp(CARP_FATAL, "Barista-XML format has not been implemented yet");
    } else {
      carp(CARP_FATAL, "Could not determine input format, "
           "Please name your files ending with .txt, .html, .sqt, .pin, "
           ".xml, .mzid, .barista.xml or use the --input-format option to "
           "specify file type");
    }
  }
  
  MatchCollection* collection = reader->parse();
  
  if (!isTabDelimited) {
    collection->setHasDistinctMatches(distinct_matches);
  } else if (collection->getHasDistinctMatches() != distinct_matches) {
    const char* matchType = collection->getHasDistinctMatches() ?
    "distinct" : "not distinct";
    carp(CARP_WARNING, "Parser has detected that matches are %s, but parameter "
         "distinct-matches is set to %s. We will assume that matches are %s",
         matchType, distinct_matches ? "distinct" : "not distinct",
         matchType);
  }
  
  carp(CARP_INFO, "Successfully read %d PSMs.", collection->getMatchTotal());
  
  // What will be used when PSMWriter is finished.
  
  PSMWriter* writer;
  stringstream output_file_name_builder;
  output_file_name_builder << output_file_base;
  
  if (output_format == "tsv") {
    output_file_name_builder << "txt";
    writer = new PMCDelimitedFileWriter();
  } else if (output_format == "html") {
    output_file_name_builder << "html";
    writer = new HTMLWriter();
  } else if (output_format == "sqt") {
    output_file_name_builder << "sqt";
    writer = new PMCSQTWriter();
  } else if (output_format == "pin") {
    output_file_name_builder << "pin";
    writer = new PinWriter();
  } else if (output_format == "pepxml") {
    output_file_name_builder << "pep.xml";
    writer = new PMCPepXMLWriter();
  } else if (output_format == "mzidentml") {
    output_file_name_builder << "mzid";
    writer = new MzIdentMLWriter();
  } else if (output_format == "barista-xml") {
    carp(CARP_FATAL, "Barista-XML format has not been implemented yet");
  } else {
    carp(CARP_FATAL, "Invalid output format.  Valid formats are: tsv, html, "
         "sqt, pin, pepxml, mzidentml, barista-xml.");
  }
  
  string output_file_name = make_file_path(output_file_name_builder.str());
  
  writer->openFile(this, output_file_name, PSMWriter::PSMS);
  writer->write(collection, database_file);
  writer->closeFile();
  
  // Clean Up
  delete collection;
  delete reader;
  delete writer;

}


int PSMConvertApplication::main(int argc, char** argv) {
  string database_file = Params::GetString("protein-database");
  string input_format = Params::GetString("input-format");
  string input_file = Params::GetString("input PSM file");
  string output_format = Params::GetString("output format");
  bool distinct_matches = Params::GetBool("distinct-matches");
  
  convertFile(input_format, output_format, input_file, "psm-convert.", database_file, distinct_matches);

  return 0;
}

string PSMConvertApplication::getName() const {
  return "psm-convert";
}

string PSMConvertApplication::getDescription() const {
  return
  "Reads in a file containing peptide-spectrum matches "
  "(PSMs) in one of the variety of supported formats and "
  "outputs the same PSMs in a different format";
}

vector<string> PSMConvertApplication::getArgs() const {
  string arr[] = {
    "input PSM file",
    "output format"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> PSMConvertApplication::getOptions() const {
  string arr[] = {
    "input-format",
    "distinct-matches",
    "protein-database",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity",
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > PSMConvertApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("psm-convert.log.txt",
    "a log file containing a copy of all messages that were printed to stderr."));
  outputs.push_back(make_pair("psm-convert.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  return outputs;
}

bool PSMConvertApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T PSMConvertApplication::getCommand() const {
  return PSM_CONVERT_COMMAND;
}

