/**
 * \file PercolatorApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef PERCOLATORAPPLICATION_H
#define PERCOLATORAPPLICATION_H

#include "CruxApplication.h"

#include <string>
#include <fstream>

class PercolatorApplication: public CruxApplication {

 protected:

  //Calls the main method in HardklorApp
  static int hardklorMain(int argc, char* argv[]);



  /**
   * writes the ISOTOPE.DAT file for hardklor
   */
  static void writeIsotopeDat(
    std::string& filename ///<path for dat file
  );

  /**
   * write the ISOTOPE.DAT to an output stream
   */
  static void writeIsotopeDat(
    std::ostream& os ///< the output stream to use
  );

  /**
   * writes the Hardklor.dat to a path
   */
  static void writeHardklorDat(
    std::string& filename ///<path to write the Hardklor.dat to
    );  

  /**
   * writes the Hardklor.dat to a stream
   */
  static void writeHardklorDat(
    std::ostream& os ///< stream to write to.
    );  

 public:

  /**
   * \returns a blank PercolatorApplication object
   */
  PercolatorApplication();

  /**
   * Destructor
   */
  ~PercolatorApplication();

  /**
   * main method for PercolatorApplication
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for PercolatorApplication
   */
  virtual std::string getName();

  /**
   * \returns the description for PercolatorApplication
   */
  virtual std::string getDescription();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();

  /**
   * \brief runs hardklor on the input spectra
   * \returns whether hardklor was successful or not
   */
  static int main(
    const std::string& input_spectra ///< file path of spectra to process
  );
  
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
