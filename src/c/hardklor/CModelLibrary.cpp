#include "CModelLibrary.h"

CModelLibrary::CModelLibrary(CAveragine* avg, CMercury8* mer){
	averagine=avg;
	mercury=mer;

	chargeMin=0;
	chargeCount=0;
	varCount=0;
	merCount=0;
}

CModelLibrary::~CModelLibrary(){
	averagine=NULL;
	mercury=NULL;
	if(libModel!=NULL) {
		eraseLibrary();
		libModel=NULL;
	}
}

bool CModelLibrary::buildLibrary(int lowCharge, int highCharge, vector<CHardklorVariant>& pepVariants){

	int i,j,k;
	unsigned int n;

	vector<Peak_T> vMR;
	Peak_T p;
	float da;
	double mass;
	char av[64];

	if(libModel!=NULL) {
		cout << "library memory already in use." << endl;
		return false;
	}

	//Fill in boundaries
	chargeMin=lowCharge;
	chargeCount=highCharge+1;
	varCount=pepVariants.size();
	merCount=400;

	libModel = new mercuryModel**[chargeCount];
	for(i=chargeMin;i<chargeCount;i++){

		libModel[i] = new mercuryModel*[varCount];
		for(j=0;j<varCount;j++){

			libModel[i][j] = new mercuryModel[merCount];
			libModel[i][j][0].area=0.0f;
			libModel[i][j][0].size=0;
			libModel[i][j][0].zeroMass=0.0;
			libModel[i][j][0].peaks=NULL;
			for(k=1;k<merCount;k++){

				mass=k*5*i-(1.007276466*i);
				averagine->clear();
				averagine->calcAveragine(mass,pepVariants[j]);
				averagine->getAveragine(&av[0]);
				mercury->GoMercury(&av[0],i);

				vMR.clear();
				da=0.0f;
				for(n=0; n<mercury->FixedData.size(); n++) {
					if(mercury->FixedData[n].data<1.0) continue;
					p.intensity=(float)mercury->FixedData[n].data;
					p.mz=mercury->FixedData[n].mass;
					da+=p.intensity;
					vMR.push_back(p);
				}
				da/=100.0f;

				libModel[i][j][k].area = da;
				libModel[i][j][k].size = vMR.size();
				libModel[i][j][k].peaks = new Peak_T[vMR.size()];
				libModel[i][j][k].zeroMass = mercury->getZeroMass();

				for(n=0;n<vMR.size();n++) libModel[i][j][k].peaks[n]=vMR[n];
			}
		}
	}

	return true;

}

void CModelLibrary::eraseLibrary(){

	int i,j,k;

	if(libModel==NULL) return;

	for(i=chargeMin;i<chargeCount;i++){
		for(j=0;j<varCount;j++){
			for(k=0;k<merCount;k++){
				delete [] libModel[i][j][k].peaks;
			}
			delete [] libModel[i][j];
		}
		delete [] libModel[i];
	}
	delete [] libModel;

	libModel=NULL;
	
}

mercuryModel* CModelLibrary::getModel(int charge, int var, double mz){

	int intMZ=(int)(mz/5);
	return &libModel[charge][var][intMZ];

}