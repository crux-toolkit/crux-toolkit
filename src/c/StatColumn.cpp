/**
 * \file StatColumn.cpp 
 * \brief Given a delimited file and a column-name, print out statistics
 * for that column (n, min, max, sum, average, stddev, median).
 *****************************************************************************/
#include "StatColumn.h"

#include "DelimitedFile.h"

using namespace std;


/**
 * \returns a blank StatColumn object
 */
StatColumn::StatColumn() {

}

/**
 * Destructor
 */
StatColumn::~StatColumn() {
}

/**
 * main method for StatColumn
 */
int StatColumn::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "delimiter",
    "header",
    "precision",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"tsv file", "column name"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */
  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  delimited_filename_ = 
    string(get_string_parameter_pointer("tsv file"));

  column_name_string_ = 
    string(get_string_parameter_pointer("column name"));

  delimiter_ = get_delimiter_parameter("delimiter");

  DelimitedFileReader delimited_file(delimited_filename_, true, delimiter_);
  
  int col_idx = delimited_file.findColumn(column_name_string_);

  if (col_idx == -1) {
    carp(CARP_ERROR,"column not found:%s\n\n%s", 
      column_name_string_.c_str(),
      delimited_file.getAvailableColumnsString().c_str());
    return(-1);
  }

  vector<FLOAT_T> data;

  FLOAT_T sum = 0;

  while (delimited_file.hasNext()) {

    FLOAT_T current = delimited_file.getFloat(col_idx);
    data.push_back(current);
    sum += current;
    delimited_file.next();

  }
  
  sort(data.begin(), data.end(), less<FLOAT_T>());

  unsigned int num_points = data.size();

  FLOAT_T min = data.front();
  FLOAT_T max = data.back();
  
  FLOAT_T average = sum / (FLOAT_T)num_points;
 
  FLOAT_T std_dev = 0.0;

  if (num_points >= 2) {

    for (unsigned int idx = 0 ; idx < num_points ; idx++) {
      FLOAT_T temp = data.at(idx) - average;
      std_dev += temp * temp;
    }
    std_dev = std_dev / (1.0 / (double)(num_points - 1));
    std_dev = sqrt(std_dev);
  }
  


 
  FLOAT_T median = 0.0;

  if (data.size() == 0) {
    carp(CARP_WARNING, "Warning no data!");
  } else {
    int half = data.size() / 2;
    if (data.size() % 2 == 0) {
      median = (data[half] + data[half-1]) / 2.0;
    } else {
      median = data[half];
    }
  }

  //print out the header

  if (header_) {
    cout <<"N\tMin\tMax\tSum\tAverage\tStdDev\tMedian"<<endl;
  }

  cout << std::setprecision(get_int_parameter("precision"));
  cout << num_points << delimiter_;
  cout << min        << delimiter_;
  cout << max        << delimiter_;
  cout << sum        << delimiter_;
  cout << average    << delimiter_;
  if (num_points >= 2) {
    cout << std_dev << delimiter_;
  } else {
    cout << "N/A" << delimiter_;
  }
  cout << median << endl;

  return 0;
}

/**
 * \returns the command name for StatColumn
 */
string StatColumn::getName() {
  return "stat-column";
}

/**
 * \returns the description for StatColumn
 */
string StatColumn::getDescription() {
  return "Collect summary statistics from a column in a tab-delimited file.";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
