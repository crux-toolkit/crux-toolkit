#include "CreateIndex.h"

#include "create_index.h"

using namespace std;

CreateIndex::CreateIndex() {

}

CreateIndex::~CreateIndex() {
}

int CreateIndex::main(int argc, char** argv) {

  return create_index_main(argc, argv);

}

string CreateIndex::getName() {
  return "create-index";
}

string CreateIndex::getDescription() {
  return "Create an index for all peptides in a fasta file.";

}
