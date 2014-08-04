/**
 * \file ExtractColumns.cpp 
 * \brief Give a tab delimited file and a comma-separated list of column names
 * print out a tab delimied file with only those columns
 *****************************************************************************/
#include "ExtractColumns.h"
#include "crux-utils.h"

using namespace std;


/**
 * \returns a blank ExtractRows object
 */
ExtractColumns::ExtractColumns() {

}

/**
 * Destructor
 */
ExtractColumns::~ExtractColumns() {
}

/**
 * main method for ExtractColumns
 */
int ExtractColumns::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "delimiter",
    "header",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"tsv file", "column names"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */
  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  const char* delimited_filename = get_string_parameter_pointer("tsv file");

  string column_names_string = string(get_string_parameter_pointer("column names"));

  char delimiter = get_delimiter_parameter("delimiter");

  DelimitedFileReader delimited_file(delimited_filename, true, delimiter);
  
  vector<string> column_name_list;
  tokenize(column_names_string, column_name_list, ',');

  vector<int> column_indices;
  for (unsigned int i=0;i<column_name_list.size();i++) {
    int col_idx = delimited_file.findColumn(column_name_list[i]);
    if (col_idx != -1) {
      column_indices.push_back(col_idx);
    } else {
      carp(CARP_ERROR,"column not found:%s\n\n%s", 
        column_name_list[i].c_str(),
        delimited_file.getAvailableColumnsString().c_str());
      return(-1);
    }
  }

  //print out the header if desired
  if (get_boolean_parameter("header")) {
    cout << column_name_list[0];
    for (unsigned int col_idx = 1;col_idx<column_name_list.size();col_idx++) {
      cout << delimiter << column_name_list[col_idx];
    }
    cout<<endl;
  }

  while(delimited_file.hasNext()) {

    int col_idx = column_indices[0];
    cout << delimited_file.getString(col_idx);
    for (unsigned int col_idx_idx = 1;col_idx_idx < column_indices.size();col_idx_idx++) {
      col_idx = column_indices[col_idx_idx];
      cout << delimiter << delimited_file.getString(col_idx);
    }
    cout <<endl;
    delimited_file.next();
  }

  return 0;

}

/**
 * \returns the command name for ExtractColumns
 */
string ExtractColumns::getName() {
  return "extract-columns";
}

/**
 * \returns the description for ExtractColumns
 */
string ExtractColumns::getDescription() {
  return "Print specified columns from a tab-delimited file.";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
