/**
 * \file ExtractColumns.cpp 
 * \brief Give a tab delimited file and a comma-separated list of column names
 * print out a tab delimied file with only those columns
 *****************************************************************************/
#include "ExtractColumns.h"
#include "util/crux-utils.h"
#include "util/Params.h"
#include "util/StringUtils.h"

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
  string delimited_filename = Params::GetString("tsv file");
  string column_names_string = Params::GetString("column names");

  char delimiter = get_delimiter_parameter("delimiter");

  DelimitedFileReader delimited_file(delimited_filename, true, delimiter);
  
  vector<string> column_name_list = StringUtils::Split(column_names_string, ',');

  vector<int> column_indices;
  for (unsigned int i = 0; i < column_name_list.size(); i++) {
    int col_idx = delimited_file.findColumn(column_name_list[i]);
    if (col_idx != -1) {
      column_indices.push_back(col_idx);
    } else {
      carp(CARP_ERROR, "column not found:%s\n\n%s", 
        column_name_list[i].c_str(),
        delimited_file.getAvailableColumnsString().c_str());
      return(-1);
    }
  }

  //print out the header if desired
  if (Params::GetBool("header")) {
    cout << column_name_list[0];
    for (unsigned int col_idx = 1; col_idx < column_name_list.size(); col_idx++) {
      cout << delimiter << column_name_list[col_idx];
    }
    cout << endl;
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
string ExtractColumns::getName() const {
  return "extract-columns";
}

/**
 * \returns the description for ExtractColumns
 */
string ExtractColumns::getDescription() const {
  return
    "[[nohtml:Print specified columns from a tab-delimited file.]]"
    "[[html:<p>Given a tab-delimited file and a comma-delimited list of column "
    "names, extract the requested columns.</p>]]";
}

/**
 * \returns the command arguments
 */
vector<string> ExtractColumns::getArgs() const {
  string arr[] = {
    "tsv file",
    "column names"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> ExtractColumns::getOptions() const {
  string arr[] = {
    "delimiter",
    "header",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > ExtractColumns::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("stdout",
    "the requested columns in tab-delimited format. The columns are printed in "
    "the same order that they appear in the input column list."));
  return outputs;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
