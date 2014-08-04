/*************************************************************************//**
 * \file ExtractRows.cpp
 * \brief Given a tab delimited file and column name and value, print
 * out all rows that pass the relation operator (default equals).
 ****************************************************************************/

#include "ExtractRows.h"

#include <iostream>

#include "DelimitedFileReader.h"
#include "DelimitedFile.h"

#include "carp.h"
#include "parameter.h"

using namespace std;

/**
 * \returns a blank ExtractRows object
 */
ExtractRows::ExtractRows() {

}

/**
 * Destructor
 */
ExtractRows::~ExtractRows() {
}


/**
 * Determines whether the a op b is true/false
 */
template<typename TValue>
bool passesThreshold(
  TValue a,
  TValue& b,
  COMPARISON_T& comparison) {

  switch (comparison) {
    case COMPARISON_LT:
      return a < b;
    case COMPARISON_LTE:
      return a <= b;
    case COMPARISON_EQ:
      return a == b;
    case COMPARISON_GTE:
      return a >= b;
    case COMPARISON_GT:
      return a > b;
    case COMPARISON_NEQ:
      return a != b;
    case NUMBER_COMPARISONS:
    case COMPARISON_INVALID:
      carp(CARP_FATAL, "Invalid comparison type!");
      return false;
  }
  return false;

}

/**
 * main method for ExtractRows
 */
int ExtractRows::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "delimiter",
    "header",
    "comparison",
    "column-type",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {
    "tsv file", 
    "column name", 
    "column value"
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */
  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  /* Get parameters */
  const char* delimited_filename = 
    get_string_parameter_pointer("tsv file");

  string column_name = 
    string(get_string_parameter_pointer("column name"));
  string column_value = 
    string(get_string_parameter_pointer("column value"));

  COLTYPE_T column_type = get_column_type_parameter("column-type");
  COMPARISON_T comparison = get_comparison_parameter("comparison");
  char delimiter = get_delimiter_parameter("delimiter");


  DelimitedFileReader delimited_file(delimited_filename, true, delimiter);
  int column_idx = delimited_file.findColumn(column_name);

  if (column_idx == -1) {
    carp(CARP_FATAL,"column not found:%s\n\n:%s",
      column_name.c_str(),
      delimited_file.getAvailableColumnsString().c_str());
  }

  if (get_boolean_parameter("header")) {
    cout << delimited_file.getHeaderString() << endl;
  }

  string column_value_str =
    string(get_string_parameter_pointer("column value"));

  int column_value_int = 0;
  from_string(column_value_int, column_value_str);

  FLOAT_T column_value_real = 0;
  from_string(column_value_real, column_value_str);


  while (delimited_file.hasNext()) {

    bool passes = false;
    switch(column_type) {
      case COLTYPE_INT:
        passes = passesThreshold(
          delimited_file.getInteger(column_idx),
          column_value_int, 
          comparison);
        break;
      case COLTYPE_REAL:
        passes = passesThreshold(
          delimited_file.getFloat(column_idx),
          column_value_real,
          comparison);
        break;
      case COLTYPE_STRING:
        passes = passesThreshold(
          delimited_file.getString(column_idx),
          column_value_str,
          comparison);
        break;
      case NUMBER_COLTYPES:
      case COLTYPE_INVALID:
        carp(CARP_FATAL,"Unknown column type");
    }
    
    if (passes) {
      cout << delimited_file.getString() << endl;
    }
    delimited_file.next();
  }

  return 0;

}

/**
 * \returns the command name for ExtractRows
 */
string ExtractRows::getName() {
  return "extract-rows";
}

/**
 * \returns the description for ExtractRows
 */
string ExtractRows::getDescription() {
  return "Print specified rows from a tab-delimited file.";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
