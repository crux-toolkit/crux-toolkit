#include "CombineOut.h"
#include "Common/TPPVersion.h"

int main(int argc, char** argv) {
  hooks_tpp(argc,argv); // installdir issues etc
  if(argc < 3) {
    cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << endl;
    cerr << " usage: CombineOut <path to directory 1 with out files> <path to directory 2 with out files>  <output path>" << endl;
    //TODO implement feature support    cerr << "    or: Sequest2XML htmlfile (prob) (-Eenzyme) (-Pparamsfile) (-Sspectrometer) (-Vschemafile) (-pI)" << endl;
    exit(1);
  }
  
  CombineOut* conv = new CombineOut(argv[1], argv[2], argv[3]);
  conv->processData();
  delete conv;
}

