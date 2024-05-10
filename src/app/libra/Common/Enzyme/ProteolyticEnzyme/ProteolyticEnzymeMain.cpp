#include <iostream>

#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#include "ProteolyticEnzyme.h"
#include "Common/Array.h"

int main(int argc, char** argv) {
  ProteolyticEnzymeFactory* factory = new ProteolyticEnzymeFactory();
  ProteolyticEnzyme* enz = NULL;

  if(! strcmp(argv[1], "digest")) {
    enz = factory->getProteolyticEnzyme(argv[2]);
    cout << "num tol term: " << enz->getNumTolTerm(argv[3][0], argv[4], argv[5][0]) << endl;
    return 0;
  }
  else if(! strcmp(argv[1], "nmc")) {
    enz = factory->getProteolyticEnzyme(argv[2]);
    if(argc > 4) {
      Array<int>* mod_pos = new Array<int>;
      int last = 0;
      char text[5];
      for(int k = 0; k < strlen(argv[4]); k++)
	if(argv[4][k] == ',') {
	  strncpy(text, argv[4]+last, k-last);
	  text[k-last] = 0;
	  mod_pos->insertAtEnd(atoi(text));
	  last = k + 1;
	}
      // last one
      strncpy(text, argv[4]+last, strlen(argv[4])-last);
      text[strlen(argv[4])-last] = 0;
      mod_pos->insertAtEnd(atoi(text));
      cout << "nmc: " << enz->getNumMissedCleavages(argv[3], mod_pos) << endl;
    }
    else
      cout << "nmc: " << enz->getNumMissedCleavages(argv[3], NULL) << endl;
    return 0;
  }

  return 0;
}
