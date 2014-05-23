#ifndef _CMODELLIBRARY_H
#define _CMODELLIBRARY_H

#include "HardklorTypes.h"
#include "CAveragine.h"
#include "CMercury8.h"
#include "CHardklorVariant.h"
#include <vector>

using namespace std;

class CModelLibrary {
public:

	//Constructors & Destructors
	CModelLibrary(CAveragine* avg, CMercury8* mer);
	~CModelLibrary();

	//User functions
	bool buildLibrary(int lowCharge, int highCharge, vector<CHardklorVariant>& pepVariants);
	void eraseLibrary();
	mercuryModel* getModel(int charge, int var, double mz);

protected:

private:

	//Data Members
	int chargeMin;
	int chargeCount;
	int varCount;
	int merCount;

	CAveragine* averagine;
	CMercury8* mercury;
	mercuryModel*** libModel;

};

#endif