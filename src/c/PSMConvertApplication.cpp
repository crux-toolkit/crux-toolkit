#include <cstdio>

#include "carp.h"
#include "parameter.h"
#include "PSMConvertApplication.h"
#include "PepXMLWriter.h"
#include "PSMReader.h"
#include "PSMWriter.h"
#include "MatchCollection.h"
#include "ProteinMatchCollection.h"
#include "PMCPepXMLWriter.h"

#include "MatchFileReader.h"
#include "MzIdentMLReader.h"
#include "PepXMLReader.h"
#include "SQTReader.h"

#include "HTMLWriter.h"
#include "MzIdentMLWriter.h"
#include "PinWriter.h"
#include "PMCDelimitedFileWriter.h"
#include "PMCPepXMLWriter.h"
#include "PMCSQTWriter.h"

PSMConvertApplication::PSMConvertApplication() {
}

PSMConvertApplication::~PSMConvertApplication() {
}

int PSMConvertApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "input-format",
    "protein-database",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity",
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "input PSM file",
    "output format"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);
  bool overwrite = get_boolean_parameter("overwrite");

  carp(CARP_INFO, "Running psm-convert...");

  string cmd_line = "crux psm-convert";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  string database_file = get_string_parameter_pointer("protein-database");

  Database* data;

/*  Database data(database_file.c_str(), false);
  // MZID TESTING
  string input_file = get_string_parameter_pointer("input PSM file");
  MatchCollection* collection = MzIdentMLReader::parse(input_file.c_str(), &data, NULL);
*/

  if (database_file.empty() || database_file.compare("__NULL_STR") == 0) {
    data = new Database();
    carp(CARP_INFO, database_file);
    carp(CARP_INFO, "Database not provided, will use empty database");
  } else {
    data = new Database(database_file.c_str(), false);
    carp(CARP_INFO, "Created Database using Fasta File");
  }
  
// trying to use NULL database for mzid ?
  Database* decoydata = new Database();

  PSMReader* reader;
  string input_format = get_string_parameter_pointer("input-format");
  string input_file = get_string_parameter_pointer("input PSM file");
  // tsv, html, sqt, pin, pepxml, mzidentml, barista-xml

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
      reader = new MatchFileReader(input_file.c_str(), data, decoydata);
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
      reader = new MatchFileReader(input_file.c_str(), data);
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

  string output_format = get_string_parameter_pointer("output format");

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

//  string pinfilename = "psm-convert.pin";
//  string outputdir = string(get_string_parameter_pointer("output-dir"));

/*  vector<MatchCollection*> decoyvec;

  PinWriter pinwriter;
  pinwriter.openFile(pinfilename.c_str(), outputdir.c_str(), true);
  pinwriter.write(collection, decoyvec, 5);
  pinwriter.closeFile();*/

// Temporary, will remove when PSMWriter is finished, so it is ok to ignore this for now.
/*
  string output_file_name = make_file_path("psm-convert.txt");
  ProteinMatchCollection protein_collection(collection);

  if (output_format.compare("tsv") == 0) {
    PMCDelimitedFileWriter writer;
    writer.openFile(this, output_file_name, PMCDelimitedFileWriter::PSMS);
    writer.write(&protein_collection);
    writer.closeFile();
  } else if (output_format.compare("html") == 0) {
    carp(CARP_FATAL, "HTML format has not been implemented yet");
  } else if (output_format.compare("sqt") == 0) {
    PMCSQTWriter writer;
    writer.openFile(this, output_file_name, PMCDelimitedFileWriter::PSMS);
    writer.write(&protein_collection, database_file);
    writer.closeFile();
  } else if (output_format.compare("pin") == 0) {
    carp(CARP_FATAL, "Pin format has not been implemented yet");
  } else if (output_format.compare("pepxml") == 0) {
    PMCPepXMLWriter writer;
    writer.openFile(output_file_name.c_str(), overwrite);
    writer.write(&protein_collection);
    writer.closeFile();
  } else if (output_format.compare("mzidentml") == 0) {
    MzIdentMLWriter writer;
    writer.openFile(output_file_name.c_str(), overwrite);
    writer.addMatches(collection);
    writer.closeFile();
  } else if (output_format.compare("barista-xml") == 0) {
    carp(CARP_FATAL, "Barista-XML format has not been implemented yet");
  } else {
      carp(CARP_FATAL, "Invalid Input Format, valid formats are: tsv, html, "
        "sqt, pin, pepxml, mzidentml, barista-xml");
  }*/

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

string PSMConvertApplication::getName() {
  return "psm-convert";
}

string PSMConvertApplication::getDescription() {
  return
  "Reads in a file containing peptide-spectrum matches "
  "(PSMs) in one of the variety of supported formats and "
  "outputs the same PSMs in a different format";
}

bool PSMConvertApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T PSMConvertApplication::getCommand() {
  return PSM_CONVERT_COMMAND;
}

