/**
 * \file CreateIndex.cpp 
 * \brief Object for running create-index
 *****************************************************************************/

#include "CreateIndex.h"

#include "create_index.h"

using namespace std;

/**
 * \returns a blank CreateIndex object
 */
CreateIndex::CreateIndex() {

}

/**
 * Destructor
 */
CreateIndex::~CreateIndex() {
}

/**
 * main method for CreateIndex
 */
int CreateIndex::main(int argc, char** argv) {

  return create_index_main(argc, argv);

}

/**
 * \returns the command name for CreateIndex
 */
string CreateIndex::getName() {
  return "create-index";
}

/**
 * \returns the description for CreateIndex
 */
string CreateIndex::getDescription() {
  return "Create an index for all peptides in a fasta file.";

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
