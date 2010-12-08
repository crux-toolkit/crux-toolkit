#include "QRanker.h"

#include "analyze_psms.h"


using namespace std;

QRanker::QRanker() {

}

QRanker::~QRanker() {
}


int QRanker::main(int argc, char** argv) {
  analyze_matches_main(QRANKER_COMMAND, argc, argv);
  return 0;

}

string QRanker::getName() {
  return "q-ranker";
}

string QRanker::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the q-ranker algorithm";

}
