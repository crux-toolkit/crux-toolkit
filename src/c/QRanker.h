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

  QRanker();
  ~QRanker();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
