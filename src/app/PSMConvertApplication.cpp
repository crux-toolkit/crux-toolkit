#include <cstdio>

#include "io/carp.h"
#include "parameter.h"
#include "PSMConvertApplication.h"
#include "io/PepXMLWriter.h"
#include "io/PSMReader.h"
#include "io/PSMWriter.h"
#include "model/MatchCollection.h"
#include "model/ProteinMatchCollection.h"
#include "io/PMCPepXMLWriter.h"

#include "io/MatchFileReader.h"
#include "io/MzIdentMLReader.h"
#include "io/PepXMLReader.h"
#include "io/SQTReader.h"

#include "io/HTMLWriter.h"
#include "io/MzIdentMLWriter.h"
#include "io/PinWriter.h"
#include "io/PMCDelimitedFileWriter.h"
#include "io/PMCPepXMLWriter.h"
#include "io/PMCSQTWriter.h"

PSMConvertApplication::PSMConvertApplication() {
}

PSMConvertApplication::~PSMConvertApplication() {
}

int PSMConvertApplication::main(int argc, char** argv) {

  bool overwrite = get_boolean_parameter("overwrite");

  carp(CARP_INFO, "Running psm-convert...");

  string cmd_line = "crux psm-convert";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  string database_file = get_string_parameter("protein-database");

  Database* data;

  if (database_file.empty() || database_file.compare("__NULL_STR") == 0) {
    data = new Database();
    carp(CARP_INFO, database_file);
    carp(CARP_INFO, "Database not provided, will use empty database");
  } else {
    data = new Database(database_file.c_str(), false);
    carp(CARP_INFO, "Created Database using Fasta File");
  }

  PSMReader* reader;
  string input_format = get_string_parameter("input-format");
  string input_file = get_string_parameter("input PSM file");

  bool isTabDelimited = false;

  if (input_format.compare("auto") == 0) {
    if (input_format.compare("tsv") == 0) {
      reader = new MatchFileReader(input_file.c_str(), data);
      isTabDelimited = true;
    } else if (input_format.compare("html") == 0) {
      carp(CARP_FATAL, "HTML format has not been implemented yet");
      reader = new MatchFileReader(input_file.c_str(), data);
    } else if (input_format.compare("sqt") == 0) {
      reader = new SQTReader(input_file.c_str(), data);
    } else if (input_format.compare("pin") == 0) {
      carp(CARP_FATAL, "Pin format has not been implemented yet");
      reader = new MatchFileReader(input_file.c_str(), data);
    } else if (input_format.compare("pepxml") == 0) {
      reader = new PepXMLReader(input_file.c_str(), data);
    } else if (input_format.compare("mzidentml") == 0) {
      reader = new MzIdentMLReader(input_file.c_str(), data);
    } else if (input_format.compare("barista-xml") == 0) {
      carp(CARP_FATAL, "Barista-XML format has not been implemented yet");
      reader = new MatchFileReader(input_file.c_str(), data);
    } else {
      carp(CARP_FATAL, "Invalid Input Format, valid formats are: tsv, html, "
        "sqt, pin, pepxml, mzidentml, barista-xml");
    }
  } else {
    if (endsWith(input_file, ".txt") == 0) {
      reader = new MatchFileReader(input_file.c_str(), data);
      isTabDelimited = true;
    } else if (endsWith(input_file, ".html") == 0) {
      carp(CARP_FATAL, "HTML format has not been implemented yet");
      reader = new MatchFileReader(input_file.c_str(), data);
    } else if (endsWith(input_file, ".sqt") == 0) {
      reader = new SQTReader(input_file.c_str(), data);
    } else if (endsWith(input_file, ".pin") == 0) {
      carp(CARP_FATAL, "Pin format has not been implemented yet");
      reader = new MatchFileReader(input_file.c_str(), data);
    } else if (endsWith(input_file, ".xml") == 0) {
      reader = new PepXMLReader(input_file.c_str(), data);
    } else if (endsWith(input_file, ".mzid") == 0) {
      reader = new MzIdentMLReader(input_file.c_str(), data);
    } else if (endsWith(input_file, ".barista.xml") == 0) {
      carp(CARP_FATAL, "Barista-XML format has not been implemented yet");
      reader = new MatchFileReader(input_file.c_str(), data);
    } else {
      carp(CARP_FATAL, "Could not determine input format, "
        "Please name your files ending with .txt, .html, .sqt, .pin, "
        ".xml, .mzid, .barista.xml or use the --input-format option to "
        "specify file type");
    }
  }

  MatchCollection* collection = reader->parse();

  if (!isTabDelimited) {
    collection->setHasDistinctMatches(get_boolean_parameter("distinct-matches"));
  } else if (collection->getHasDistinctMatches() != get_boolean_parameter("distinct-matches")) {
    carp(CARP_WARNING, "Parser has detected that matches are %s, but parameter "
      "distinct-matches is set to %s. We will assume that matches are %s",
      collection->getHasDistinctMatches() ? "distinct" : "not distinct",
      get_boolean_parameter("distinct-matches") ? "distinct" : "not distinct",
      collection->getHasDistinctMatches() ? "distinct" : "not distinct");
  }

  carp(CARP_INFO, "Reader has been succesfully parsed");

  string output_format = get_string_parameter("output format");

// What will be used when PSMWriter is finished.

  PSMWriter* writer;
  stringstream output_file_name_builder;
  output_file_name_builder << "psm-convert.";

  if (output_format.compare("tsv") == 0) {
    output_file_name_builder << "txt";
    writer = new PMCDelimitedFileWriter();
  } else if (output_format.compare("html") == 0) {
    output_file_name_builder << "html";
    writer = new HTMLWriter();
  } else if (output_format.compare("sqt") == 0) {
    output_file_name_builder << "sqt";
    writer = new PMCSQTWriter();
  } else if (output_format.compare("pin") == 0) {
    output_file_name_builder << "pin";
    writer = new PinWriter();
  } else if (output_format.compare("pepxml") == 0) {
    output_file_name_builder << "pep.xml";
    writer = new PMCPepXMLWriter();
  } else if (output_format.compare("mzidentml") == 0) {
    output_file_name_builder << "mzid";
    writer = new MzIdentMLWriter();
  } else if (output_format.compare("barista-xml") == 0) {
    carp(CARP_FATAL, "Barista-XML format has not been implemented yet");
    output_file_name_builder << "barista.xml";
    writer = new PMCDelimitedFileWriter();
  } else {
      carp(CARP_FATAL, "Invalid Input Format, valid formats are: tsv, html, "
        "sqt, pin, pepxml, mzidentml, barista-xml");
    writer = new PMCDelimitedFileWriter();
  }

  string output_file_name = make_file_path(output_file_name_builder.str());

  writer->openFile(this, output_file_name.c_str(), PSMWriter::PSMS);
  writer->write(collection, database_file);
  writer->closeFile();

// Clean Up
  delete collection;
  delete reader;
  delete writer;

  return 0;
}

// Accepts two strings, returns whether the first string ends with the
// second string
int PSMConvertApplication::endsWith(string s, string ending) {
  if (s.length() >= ending.length()) {
    return s.compare(s.length() - ending.length(), ending.length(), ending);
  }
  return 1;
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
    "protein-database",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity",
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

bool PSMConvertApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T PSMConvertApplication::getCommand() {
  return PSM_CONVERT_COMMAND;
}

