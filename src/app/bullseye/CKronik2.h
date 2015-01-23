#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>

using namespace std;

typedef struct sPep{
  int charge;
  float intensity;
  double monoMass;
	double basePeak;
  double xCorr;
	char mods[32];
} sPep;

typedef struct sScan {
  vector<sPep> *vPep;
	int scanNum;
	char file[256];
	float rTime;

	//Constructors & Destructor
	sScan(){vPep = new vector<sPep>;}
	sScan(const sScan& s){
		vPep = new vector<sPep>;
		for(unsigned int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
		scanNum = s.scanNum;
		strcpy(file,s.file);
		rTime=s.rTime;
	}
	~sScan(){delete vPep;}

  //Copy operator
	sScan& operator=(const sScan& s){
		if(this!=&s){
			delete vPep;
			vPep = new vector<sPep>;
			for(unsigned int i=0;i<s.vPep->size();i++) vPep->push_back(s.vPep->at(i));
			scanNum = s.scanNum;
			strcpy(file,s.file);
			rTime=s.rTime;
		}
		return *this;
	}

  //Clear
	void clear(){
		delete vPep;
		vPep = new vector<sPep>;
	}

  static int compareIntRev(const void *p1, const void *p2){
    const sPep d1 = *(sPep *)p1;
    const sPep d2 = *(sPep *)p2;
    if(d1.intensity>d2.intensity) return -1;
    else if(d1.intensity<d2.intensity) return 1;
    else return 0;
  }

  void sortIntRev(){
    if(vPep->size()==0) return;
    qsort(&vPep->at(0),vPep->size(),sizeof(sPep),compareIntRev);
  }

} sScan;

typedef struct sProfileData {
  bool interpolated;

  int scanNum;

  float intensity;
  float rTime;

  double monoMass;
  double xCorr;

  //Constructor
  sProfileData(){
    interpolated=false;
    scanNum=0;
    intensity=0.0f;
    rTime=0.0f;
    monoMass=0.0;
    xCorr=0.0;
  }
} sProfileData;

typedef struct sPepProfile {
  int lowScan;
  int highScan;
  int bestScan;
  int charge;
  int MS2Events;

  unsigned int datapoints;

  float rTime;
  float firstRTime;
  float lastRTime;
  float intensity;
  float sumIntensity;

  double monoMass;
  double basePeak;
  double xCorr;

  char mods[32];
  char gene[32];
  char sequence[64];

  sProfileData* profile;
  //sProfileData* gaussProfile;

  //Constructors
  sPepProfile() {
    profile = new sProfileData[1];
    //gaussProfile = new sProfileData[1];
    sProfileData d;
    profile[0]=d;
    //gaussProfile[0]=d;
    datapoints=0;
    lowScan=0;
    highScan=0;
    bestScan=0;
    MS2Events=0;
    charge=0;
    rTime=0.0f;
    firstRTime=0.0f;
    lastRTime=0.0f;
    intensity=0.0f;
    sumIntensity=0.0f;
    monoMass=0.0;
    basePeak=0.0;
    xCorr=0.0;
    strcpy(mods,"");
    strcpy(sequence,"");
    strcpy(gene,"");
  }
  sPepProfile(const unsigned int i) {
    profile = new sProfileData[i];
    //gaussProfile = new sProfileData[i];
    sProfileData d;
    unsigned int j;
    for(j=0;j<i;j++) profile[j]=d;
    datapoints=i;
    lowScan=0;
    highScan=0;
    bestScan=0;
    charge=0;
    MS2Events=0;
    rTime=0.0f;
    firstRTime=0.0f;
    lastRTime=0.0f;
    intensity=0.0f;
    sumIntensity=0.0f;
    monoMass=0.0;
    basePeak=0.0;
    xCorr=0.0;
    strcpy(mods,"");
    strcpy(sequence,"");
    strcpy(gene,"");
  }
  sPepProfile(const sPepProfile& p){
    unsigned int i;
    profile = new sProfileData[p.datapoints];
    for(i=0;i<p.datapoints;i++) profile[i]=p.profile[i];
    lowScan = p.lowScan;
    highScan = p.highScan;
    bestScan = p.bestScan;
    charge = p.charge;
    MS2Events = p.MS2Events;
    datapoints = p.datapoints;
    rTime = p.rTime;
    firstRTime = p.firstRTime;
    lastRTime = p.lastRTime;
    intensity = p.intensity;
    sumIntensity = p.sumIntensity;
    monoMass = p.monoMass;
    basePeak = p.basePeak;
    xCorr = p.xCorr;
    strcpy(mods,p.mods);
    strcpy(sequence,p.sequence);
    strcpy(gene,p.gene);
  }

  //Destructor
  ~sPepProfile(){
    delete [] profile;
  }

  //Copy operator
  sPepProfile& operator=(const sPepProfile& p){
    if(this!=&p){
      unsigned int i;
      delete [] profile;
      profile = new sProfileData[p.datapoints];
      for(i=0;i<p.datapoints;i++) profile[i]=p.profile[i];
      lowScan = p.lowScan;
      highScan = p.highScan;
      bestScan = p.bestScan;
      MS2Events = p.MS2Events;
      charge = p.charge;
      datapoints = p.datapoints;
      rTime = p.rTime;
      firstRTime = p.firstRTime;
      lastRTime = p.lastRTime;
      intensity = p.intensity;
      sumIntensity = p.sumIntensity;
      monoMass = p.monoMass;
      basePeak = p.basePeak;
      xCorr = p.xCorr;
      strcpy(mods,p.mods);
      strcpy(sequence,p.sequence);
      strcpy(gene,p.gene);
    }
    return *this;
  }

  void setPoints(const unsigned int& i){
    delete [] profile;
    profile = new sProfileData[i];
    datapoints=i;
    sProfileData d;
    for(unsigned int j=0;j<i;j++) profile[j]=d;
  }

  static int compareScanNum(const void *p1, const void *p2){
    const sProfileData d1 = *(sProfileData *)p1;
    const sProfileData d2 = *(sProfileData *)p2;
    if(d1.scanNum<d2.scanNum) return -1;
    else if(d1.scanNum>d2.scanNum) return 1;
    else return 0;
  }

  static int compareIntRev(const void *p1, const void *p2){
    const sProfileData d1 = *(sProfileData *)p1;
    const sProfileData d2 = *(sProfileData *)p2;
    if(d1.intensity>d2.intensity) return -1;
    else if(d1.intensity<d2.intensity) return 1;
    else return 0;
  }

  void sortIntRev(){
    if(datapoints==0) return;
    qsort(profile,datapoints,sizeof(sProfileData),compareIntRev);
  }

  void sortScanNum(){
    if(datapoints==0) return;
    qsort(profile,datapoints,sizeof(sProfileData),compareScanNum);
  }

} sPepProfile;

typedef struct iTwo{
  int scan;
  int pep;
} iTwo;

class CKronik2 {
public:

  //Constructors and destructors
  CKronik2();
  CKronik2(const CKronik2& c);
  ~CKronik2();

  //Operator overrides
  CKronik2& operator=(const CKronik2& c);
  sPepProfile& operator[ ](const unsigned int& i);

  //Statistics functions  
  bool pearson(int i1, int i2, bool byScan, bool interpolate, double& rval, double& pval, double& slope, double& intercept, float& sum1, float& sum2);
  
  //Parameter accessors and modifiers
  void setPPMTol(double d);
  void setMatchTol(int i);
  void setGapTol(int i);

  //Automation
  int getPercent();
  bool loadHK(char* in);
  bool processHK(char* in, char* out="\0");

  //Tools
  bool getRT(int scanNum, float& rt);

  //Filters
  void filterRT(float rt1, float rt2);
  void removeContaminants(float rt);
	void removeMass(double m1, double m2);

  //Data manipulation functions
  sPepProfile& at(unsigned int i);
  void add(sPepProfile& p);
  void clear();
  void clearHK();
  void erase(unsigned int i);
  unsigned int size();

  //Sorting functions
  void sortBasePeak();
  void sortMonoMass();
  void sortFirstRTime();
  void sortIntensityRev();

protected:
private:
  bool findMax(vector<sScan>& v, int& s, int& p);
  double interpolate(int x1, int x2, double y1, double y2, int x);
  
  //Statistics functions
  double betai(double a, double b, double x);
  double betacf(double a, double b, double x);
  double gammaln(double xx);

  //Data Members: data analysis
  vector<sPepProfile> vPeps;
  vector<sScan> hkData;

  //Data Members: parameters
  double dPPMTol;   //default 10.0
  int iGapTol;      //default 1
  int iMatchTol;    //Default 3
  int iPercent;

  //Sorting Functions
  void sortPeptide();
  static int compareBP(const void *p1, const void *p2);
  static int compareMM(const void *p1, const void *p2);
  static int compareFRT(const void *p1, const void *p2);
  static int compareIRev(const void *p1, const void *p2);

};
