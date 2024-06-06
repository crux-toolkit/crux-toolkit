#include "CompactParser.h"

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  CompactParser *p = new CompactParser(argv[1]);
  delete p;

  return 0;

}
