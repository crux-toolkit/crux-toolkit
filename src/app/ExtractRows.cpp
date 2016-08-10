/*******************************************************************************
  * \file ExtractRows.cpp
  * \brief Given a tab delimited file and column name and value, print
  * out all rows that pass the relation operator (default equals).
  ****************************************************************************/

#include "ExtractRows.h"

#include <iostream>

#include "io/DelimitedFileReader.h"
#include "io/DelimitedFile.h"
#include "io/carp.h"
#include "parameter.h"
#include "util/Params.h"
#include "util/StringUtils.h"

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
  TValue b,
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
  /* Get parameters */
  string delimited_filename = Params::GetString("tsv file");

  string column_name = Params::GetString("column name");
  string column_value = Params::GetString("column value");

  COLTYPE_T column_type = get_column_type_parameter("column-type");
  COMPARISON_T comparison = get_comparison_parameter("comparison");
  char delimiter = get_delimiter_parameter("delimiter");


  DelimitedFileReader delimited_file(delimited_filename, true, delimiter);
  int column_idx = delimited_file.findColumn(column_name);

  if (column_idx == -1) {
    carp(CARP_FATAL, "column not found:%s\n\n:%s",
      column_name.c_str(),
      delimited_file.getAvailableColumnsString().c_str());
  }

  if (Params::GetBool("header")) {
    cout << delimited_file.getHeaderString() << endl;
  }

  string column_value_str = Params::GetString("column value");

  while (delimited_file.hasNext()) {
    bool passes = false;
    switch(column_type) {
      case COLTYPE_INT:
        passes = passesThreshold(
          delimited_file.getInteger(column_idx),
          StringUtils::FromString<int>(column_value_str), 
          comparison);
        break;
      case COLTYPE_REAL:
        passes = passesThreshold(
          delimited_file.getFloat(column_idx),
          StringUtils::FromString<FLOAT_T>(column_value_str), 
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
        carp(CARP_FATAL, "Unknown column type");
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
string ExtractRows::getName() const {
  return "extract-rows";
}

/**
 * \returns the description for ExtractRows
 */
string ExtractRows::getDescription() const {
  return
    "[[nohtml:Print specified rows from a tab-delimited file.]]"
    "[[html:Given a tab-delimited file, a column name and a column cell value, "
    "extract the rows that have the matching values for that column.]]";
}

/**
 * \returns the command arguments
 */
vector<string> ExtractRows::getArgs() const {
  string arr[] = {
    "tsv file", 
    "column name", 
    "column value"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> ExtractRows::getOptions() const {
  string arr[] = {
    "delimiter",
    "header",
    "comparison",
    "column-type",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > ExtractRows::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("stdout",
    "the rows for which the value observed in the specified column of the input "
    "file match the <column value> given on the command line."));
  return outputs;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
