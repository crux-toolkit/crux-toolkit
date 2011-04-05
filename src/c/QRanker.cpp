/**
 * \file QRanker.cpp
 * \brief Object for running qranker
 *****************************************************************************/

#include "QRanker.h"

#include "analyze_psms.h"


using namespace std;

/**
 * \returns a blank QRanker object
 */
QRanker::QRanker() {

}

/**
 * Destructor
 */
QRanker::~QRanker() {
}

/**
 * main method for QRanker
 */
int QRanker::main(int argc, char** argv) {
  analyze_matches_main(QRANKER_COMMAND, argc, argv);
  return 0;

}

/**
 * \returns the command name for QRanker
 */
string QRanker::getName() {
  return "q-ranker";
}

/**
 * \returns the description for QRanker
 */
string QRanker::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the q-ranker algorithm";

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
