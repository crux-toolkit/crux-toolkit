/**
 * \file ExtractColumns.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Give a tab delimited file and a comma-separated list of column names
 * print out a tab delimied file with only those columns
 *****************************************************************************/
#ifndef EXTRACTCOLUMNS_H
#define EXTRACTCOLUMNS_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class ExtractColumns: public CruxApplication {

 public:

  /**
   * \returns a blank ExtractRows object
   */
  ExtractColumns();

  /**
   * Destructor
   */
  ~ExtractColumns();

  /**
   * main method for ExtractColumns
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for ExtractColumns
   */
  virtual std::string getName();

  /**
   * \returns the description for ExtractColumns
   */
  virtual std::string getDescription();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
