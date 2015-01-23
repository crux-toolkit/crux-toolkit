#include <fstream>
#include <string>

#include "CHardklor.h"
#include "CHardklor2.h"
#include "CModelLibrary.h"
#include "CHardklorParser.h"
#include "CHardklorSetting.h"
#include "CHardklorVariant.h"

#include "CruxHardklorApplication.h"

using namespace std;


#ifdef CRUX
int CruxHardklorApplication::hardklorMain(int argc, char* argv[]) {
#else
int main(int argc, char* argv[]) {
#endif
  int i;
	unsigned int j;
  bool bConf;
  char tstr[512]="\0";

	CAveragine *averagine;
	CMercury8 *mercury;
	CModelLibrary *models;

	if(argc==1){
		cout << "Hardklor v2.01, September 23, 2011\n";
		cout << "Usage:\t\thardklor -conf <config file>\n";
		cout << "See documentation for instructions to modify and use config files." << endl;
		exit(1);
	}
  
  bConf = false;
  for(i=0;i<argc;i++){
    if(strcmp(argv[i],"-conf")==0) {
      bConf = true;
      break;
    }
  }
  
  
  CHardklorParser hp;
  if(bConf) {
    hp.parseConfig(argv[i+1]);
  } else {
    for(i=1;i<argc;i++){
      strcat(tstr," ");
      strcat(tstr,argv[i]);
    }
    hp.parse(tstr);
  }
  

  //Create all the output files that will be used
  for(i=0;i<hp.size();i++){
    fstream fptr(&hp.queue(i).outFile[0],fstream::out);
    fptr.clear();
    fptr.close();
  }

	averagine = new CAveragine(hp.queue(0).MercuryFile,hp.queue(0).HardklorFile);
	mercury = new CMercury8(hp.queue(0).MercuryFile);
	models = new CModelLibrary(averagine,mercury);

  CHardklor h(averagine,mercury);
	CHardklor2 h2(averagine,mercury,models);
	vector<CHardklorVariant> pepVariants;
	CHardklorVariant hkv;

  for(i=0;i<hp.size();i++) {
		if(hp.queue(i).sna>3){
			pepVariants.clear();
			if(!hp.queue(i).noBase) pepVariants.push_back(hkv);
			for(j=0;j<hp.queue(i).variant->size();j++)  pepVariants.push_back(hp.queue(i).variant->at(j));

			models->eraseLibrary();
			models->buildLibrary(hp.queue(i).minCharge,hp.queue(i).maxCharge,pepVariants);
			h2.GoHardklor(hp.queue(i));
		} else {
			h.GoHardklor(hp.queue(i));
		}
  }

	delete models;
	delete averagine;
	delete mercury;
  
  return 0;
  
}
