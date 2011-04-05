/**
 * \file QRanker.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running qranker
 *****************************************************************************/
#ifndef QRANKER_H
#define QRANKER_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class QRanker: public CruxApplication {

 public:
  /**
   * \returns a blank QRanker object
   */
  QRanker();
  
  /**
   * Destructor
   */
  ~QRanker();

  /**
   * main method for QRanker
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for QRanker
   */
  virtual std::string getName();

  /**
   * \returns the description for QRanker
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
