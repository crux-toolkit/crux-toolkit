#include "CKronik2.h"

//-------------------------------------
//   Constructors and Destructors
//-------------------------------------
CKronik2::CKronik2(){
  dPPMTol   = 10.0; 
  iGapTol   = 1;   
  iMatchTol = 3;    
}
CKronik2::~CKronik2(){
}

CKronik2::CKronik2(const CKronik2& c){
  dPPMTol=c.dPPMTol;
  iGapTol=c.iGapTol;
  iMatchTol=c.iMatchTol;
  iPercent=c.iPercent;
  vPeps.clear();
  for(unsigned int i=0;i<c.vPeps.size();i++) vPeps.push_back(c.vPeps[i]);
  hkData.clear();
  for(unsigned int i=0;i<c.hkData.size();i++) hkData.push_back(c.hkData[i]);
}


//-------------------------------------
//    Operator Overloads
//-------------------------------------
CKronik2& CKronik2::operator=(const CKronik2& c){
  if(this!=&c){
    dPPMTol=c.dPPMTol;
    iGapTol=c.iGapTol;
    iMatchTol=c.iMatchTol;
    iPercent=c.iPercent;
    vPeps.clear();
    for(unsigned int i=0;i<c.vPeps.size();i++) vPeps.push_back(c.vPeps[i]);
    hkData.clear();
    for(unsigned int i=0;i<c.hkData.size();i++) hkData.push_back(c.hkData[i]);
  } 
  return *this;
}
sPepProfile& CKronik2::operator[ ](const unsigned int& i){
  return vPeps[i];
}

//-------------------------------------
//    Data Manipulation
//-------------------------------------
void CKronik2::add(sPepProfile &p){
  vPeps.push_back(p);
}

sPepProfile& CKronik2::at(unsigned int i){
  return vPeps[i];
}

void CKronik2::clear(){
  vPeps.clear();
  hkData.clear();
}

void CKronik2::clearHK(){
  hkData.clear();
}

void CKronik2::erase(unsigned int i){
  vPeps.erase(vPeps.begin()+i);
}

//------------------------------
//   Automation
//------------------------------
int CKronik2::getPercent(){
  return iPercent;
}
bool CKronik2::loadHK(char* in){
	FILE *hkr;
	sScan scan;
	sPep pep;
	double td;
	char tag;
	bool firstScan;
	char line[256];
	char* tok;

	int pepCount=0;

  sPepProfile s;
  sProfileData p;

  //Read in the Hardklor results
  firstScan=true;
	hkr = fopen(in,"rt");
	if(hkr==NULL) {
		cout << "Problem reading file." << endl;
		return false;
	}

  hkData.clear();
	while(!feof(hkr)){
    
		tag=fgetc(hkr);

    if(tag=='S') {
			if(firstScan) firstScan=false;
			else hkData.push_back(scan);	
			scan.clear();
			fgets(line,256,hkr);
			tok=strtok(line,"\t\n");
			scan.scanNum=atoi(tok);
			tok=strtok(NULL,"\t\n");
			scan.rTime=(float)atof(tok);
			tok=strtok(NULL,"\t\n");
			strcpy(scan.file,tok);
      //fscanf(hkr,"\t%d\t%f%s\n",&scan.scanNum,&scan.rTime,scan.file);
		} else {
			pepCount++;
			fscanf(hkr,"\t%lf\t%d\t%f\t%lf\t%lf-%lf\t%lf\t%s\t%lf\n", &pep.monoMass,&pep.charge,&pep.intensity,&pep.basePeak,&td,&td,&td,pep.mods,&pep.xCorr);
			scan.vPep->push_back(pep);
		}
	}
  hkData.push_back(scan);
	fclose(hkr);

  cout << pepCount << " peptides from " << hkData.size() << " scans." << endl;
}

bool CKronik2::processHK(char*  in, char* out) {
	FILE *hkr;
	sScan scan;
	sPep pep;
	double td;
	char tag;
	bool firstScan;
  int sIndex,pIndex;
  int i,j,k,k1,k2;

	char line[256];
	char* tok;

	int pepCount=0;
  vector<sScan> allScans;

  double mass;
  double ppm;
  int charge;
  int gap;
  int matchCount;
  bool bMatch;

  sPepProfile s;
  sProfileData p;

  //for tracking which peptides
  iTwo t;
  vector<iTwo> vLeft;
  vector<iTwo> vRight;

  //clear data
  vPeps.clear();

  //Read in the Hardklor results
  firstScan=true;
	hkr = fopen(in,"rt");
	if(hkr==NULL) {
		cout << "Problem reading file." << endl;
		return false;
	}

	while(!feof(hkr)){
    
		tag=fgetc(hkr);

    if(tag=='S') {
			if(firstScan) firstScan=false;
			else allScans.push_back(scan);	
			scan.clear();
			fgets(line,256,hkr);
			tok=strtok(line,"\t\n");
			scan.scanNum=atoi(tok);
			tok=strtok(NULL,"\t\n");
			scan.rTime=(float)atof(tok);
			tok=strtok(NULL,"\t\n");
			strcpy(scan.file,tok);
      //fscanf(hkr,"\t%d\t%f%s\n",&scan.scanNum,&scan.rTime,scan.file);
		} else {
			pepCount++;
			fscanf(hkr,"\t%lf\t%d\t%f\t%lf\t%lf-%lf\t%lf\t%s\t%lf\n", &pep.monoMass,&pep.charge,&pep.intensity,&pep.basePeak,&td,&td,&td,pep.mods,&pep.xCorr);
			scan.vPep->push_back(pep);
		}
	}
  allScans.push_back(scan);
	fclose(hkr);

  cout << pepCount << " peptides from " << allScans.size() << " scans." << endl;

  for(i=0;i<allScans.size();i++) allScans[i].sortIntRev();

  cout << "Finding persistent peptide signals:" << endl;

  int startCount=pepCount;
  int lastPercent=0;

  cout << lastPercent;

  //Perform the Kronik analysis
  while(pepCount>0){
    if(!findMax(allScans,sIndex,pIndex)) break;

    mass=allScans[sIndex].vPep->at(pIndex).monoMass;
    charge=allScans[sIndex].vPep->at(pIndex).charge;
    matchCount=1;

    //look left
    vLeft.clear();
    gap=0;
    i=sIndex-1;
    while(i>-1 && gap<=iGapTol){
      bMatch=false;
      t.scan=i;
      t.pep=-1;
      for(j=0;j<allScans[i].vPep->size();j++){
        //if(allScans[i].vPep->at(j).intensity<0.0) continue;
        ppm=(allScans[i].vPep->at(j).monoMass-mass)/mass*1000000;
        if(fabs(ppm)<dPPMTol && allScans[i].vPep->at(j).charge==charge){
          t.pep=j;
          gap=0;
          bMatch=true;
          matchCount++;
          break;
        }
      }
      if(!bMatch) gap++;
      vLeft.push_back(t);
      i--;
    }

    //look right
    vRight.clear();
    gap=0;
    i=sIndex+1;
    while(i<allScans.size() && gap<=iGapTol){    
      bMatch=false;
      t.scan=i;
      t.pep=-1;
      for(j=0;j<allScans[i].vPep->size();j++){
        //if(allScans[i].vPep->at(j).intensity<0.0) continue;
        ppm=(allScans[i].vPep->at(j).monoMass-mass)/mass*1000000;
        if(fabs(ppm)<dPPMTol && allScans[i].vPep->at(j).charge==charge){
          t.pep=j;
          gap=0;
          bMatch=true;
          matchCount++;
          break;
        }
      }
      if(!bMatch) gap++;
      vRight.push_back(t);
      i++;
    }

    //Only keep sufficient matches
    if(matchCount>=iMatchTol){
      //trim gaps
      while(vLeft.size()>0 && vLeft[vLeft.size()-1].pep<0) vLeft.pop_back();
      while(vRight.size()>0 && vRight[vRight.size()-1].pep<0) vRight.pop_back();

      //apply basic information
      s.rTime=allScans[sIndex].rTime;
      s.basePeak=allScans[sIndex].vPep->at(pIndex).basePeak;
      s.bestScan=allScans[sIndex].scanNum;
      s.charge=allScans[sIndex].vPep->at(pIndex).charge;
      s.intensity=allScans[sIndex].vPep->at(pIndex).intensity;
      s.monoMass=allScans[sIndex].vPep->at(pIndex).monoMass;
      strcpy(s.mods,allScans[sIndex].vPep->at(pIndex).mods);
      s.xCorr=allScans[sIndex].vPep->at(pIndex).xCorr;
      s.setPoints(vLeft.size()+vRight.size()+1);
      if(vLeft.size()==0) {
        s.lowScan=allScans[sIndex].scanNum;
        s.firstRTime=allScans[sIndex].rTime;
      } else {
        s.lowScan=allScans[vLeft[vLeft.size()-1].scan].scanNum;
        s.firstRTime=allScans[vLeft[vLeft.size()-1].scan].rTime;
      }
      if(vRight.size()==0) {
        s.highScan=allScans[sIndex].scanNum;
        s.lastRTime=allScans[sIndex].rTime;
      } else {
        s.highScan=allScans[vRight[vRight.size()-1].scan].scanNum;
        s.lastRTime=allScans[vRight[vRight.size()-1].scan].rTime;
      }

      //apply datapoints
      s.profile[0].intensity=allScans[sIndex].vPep->at(pIndex).intensity;
      s.profile[0].interpolated=false;
      s.profile[0].monoMass=allScans[sIndex].vPep->at(pIndex).monoMass;
      s.profile[0].rTime=allScans[sIndex].rTime;
      s.profile[0].scanNum=allScans[sIndex].scanNum;
      s.profile[0].xCorr=allScans[sIndex].vPep->at(pIndex).xCorr;

      i=1;
      j=0;
      while(j<vLeft.size()){
        //Handle interpolated data
        if(vLeft[j].pep<0) {
          k=0;
          while(vLeft[j+k].pep<0) k++;
          k2=j+k;   
          if(j==0){
            s.profile[i].intensity=interpolate(s.profile[0].scanNum,allScans[vLeft[k2].scan].scanNum,(double)s.profile[0].intensity,(double)allScans[vLeft[k2].scan].vPep->at(vLeft[k2].pep).intensity,allScans[vLeft[j].scan].scanNum);
            s.profile[i].monoMass =interpolate(s.profile[0].scanNum,allScans[vLeft[k2].scan].scanNum,s.profile[0].monoMass, allScans[vLeft[k2].scan].vPep->at(vLeft[k2].pep).monoMass, allScans[vLeft[j].scan].scanNum);
          } else {
            k=0;
            while(j+k>=0 && vLeft[j+k].pep<0) k--;
            k1=j+k;
            if(k1<0){
              s.profile[i].intensity=interpolate(s.profile[0].scanNum,allScans[vLeft[k2].scan].scanNum,(double)s.profile[0].intensity,(double)allScans[vLeft[k2].scan].vPep->at(vLeft[k2].pep).intensity,allScans[vLeft[j].scan].scanNum);
              s.profile[i].monoMass =interpolate(s.profile[0].scanNum,allScans[vLeft[k2].scan].scanNum,s.profile[0].monoMass, allScans[vLeft[k2].scan].vPep->at(vLeft[k2].pep).monoMass, allScans[vLeft[j].scan].scanNum);
            } else {
              s.profile[i].intensity=interpolate(allScans[vLeft[k1].scan].scanNum,allScans[vLeft[k2].scan].scanNum,(double)allScans[vLeft[k1].scan].vPep->at(vLeft[k1].pep).intensity,(double)allScans[vLeft[k2].scan].vPep->at(vLeft[k2].pep).intensity,allScans[vLeft[j].scan].scanNum);
              s.profile[i].monoMass =interpolate(allScans[vLeft[k1].scan].scanNum,allScans[vLeft[k2].scan].scanNum,allScans[vLeft[k1].scan].vPep->at(vLeft[k1].pep).monoMass, allScans[vLeft[k2].scan].vPep->at(vLeft[k2].pep).monoMass, allScans[vLeft[j].scan].scanNum);
            }
          }    
          s.profile[i].interpolated=true;
          s.profile[i].rTime=allScans[vLeft[j].scan].rTime;
          s.profile[i].scanNum=allScans[vLeft[j].scan].scanNum;
          s.profile[i].xCorr=0.0;
        } else {
          s.profile[i].intensity=allScans[vLeft[j].scan].vPep->at(vLeft[j].pep).intensity;
          s.profile[i].interpolated=false;
          s.profile[i].monoMass=allScans[vLeft[j].scan].vPep->at(vLeft[j].pep).monoMass;
          s.profile[i].rTime=allScans[vLeft[j].scan].rTime;
          s.profile[i].scanNum=allScans[vLeft[j].scan].scanNum;
          s.profile[i].xCorr=allScans[vLeft[j].scan].vPep->at(vLeft[j].pep).xCorr;
        }
        i++;
        j++;
      }

      j=0;
      while(j<vRight.size()){
        //Handle interpolated data
        if(vRight[j].pep<0) {
          k=0;
          while(vRight[j+k].pep<0)k++;
          k2=j+k;
          if(j==0){
            s.profile[i].intensity=interpolate(s.profile[0].scanNum,allScans[vRight[k2].scan].scanNum,s.profile[0].intensity,allScans[vRight[k2].scan].vPep->at(vRight[k2].pep).intensity,allScans[vRight[j].scan].scanNum);
            s.profile[i].monoMass =interpolate(s.profile[0].scanNum,allScans[vRight[k2].scan].scanNum,s.profile[0].monoMass, allScans[vRight[k2].scan].vPep->at(vRight[k2].pep).monoMass, allScans[vRight[j].scan].scanNum);
          } else {
            k=0;
            while((j+k>-1) && vRight[j+k].pep<0)k--;
            k1=j+k;
            if(k1<0){
              s.profile[i].intensity=interpolate(s.profile[0].scanNum,allScans[vRight[k2].scan].scanNum,s.profile[0].intensity,allScans[vRight[k2].scan].vPep->at(vRight[k2].pep).intensity,allScans[vRight[j].scan].scanNum);
              s.profile[i].monoMass =interpolate(s.profile[0].scanNum,allScans[vRight[k2].scan].scanNum,s.profile[0].monoMass, allScans[vRight[k2].scan].vPep->at(vRight[k2].pep).monoMass, allScans[vRight[j].scan].scanNum);
            } else {
              s.profile[i].intensity=interpolate(allScans[vRight[k1].scan].scanNum,allScans[vRight[k2].scan].scanNum,allScans[vRight[k1].scan].vPep->at(vRight[k1].pep).intensity,allScans[vRight[k2].scan].vPep->at(vRight[k2].pep).intensity,allScans[vRight[j].scan].scanNum);
              s.profile[i].monoMass =interpolate(allScans[vRight[k1].scan].scanNum,allScans[vRight[k2].scan].scanNum,allScans[vRight[k1].scan].vPep->at(vRight[k1].pep).monoMass, allScans[vRight[k2].scan].vPep->at(vRight[k2].pep).monoMass, allScans[vRight[j].scan].scanNum);
            }
          } 
          s.profile[i].interpolated=true;
          s.profile[i].rTime=allScans[vRight[j].scan].rTime;
          s.profile[i].scanNum=allScans[vRight[j].scan].scanNum;
          s.profile[i].xCorr=0.0;
        } else {
          s.profile[i].intensity=allScans[vRight[j].scan].vPep->at(vRight[j].pep).intensity;
          s.profile[i].interpolated=false;
          s.profile[i].monoMass=allScans[vRight[j].scan].vPep->at(vRight[j].pep).monoMass;
          s.profile[i].rTime=allScans[vRight[j].scan].rTime;
          s.profile[i].scanNum=allScans[vRight[j].scan].scanNum;
          s.profile[i].xCorr=allScans[vRight[j].scan].vPep->at(vRight[j].pep).xCorr;
        }
        i++;
        j++;
      }

      //Add this peptide
      s.sortScanNum();

      //summed intensity (including interpolation)
      s.sumIntensity=0.0f;
      for(i=0;i<s.datapoints;i++) s.sumIntensity+=s.profile[i].intensity;

      vPeps.push_back(s);

      //Erase datapoints already used
      for(i=0;i<vLeft.size();i++){
        if(vLeft[i].pep<0) continue;
        allScans[vLeft[i].scan].vPep->erase(allScans[vLeft[i].scan].vPep->begin()+vLeft[i].pep);
        pepCount--;
      }
      for(i=0;i<vRight.size();i++){
        if(vRight[i].pep<0) continue;
        allScans[vRight[i].scan].vPep->erase(allScans[vRight[i].scan].vPep->begin()+vRight[i].pep);
        pepCount--;
      }
    }

    //erase the one we're looking at
    allScans[sIndex].vPep->erase(allScans[sIndex].vPep->begin()+pIndex);
    pepCount--;

    //update percent
    iPercent=100-(int)((float)pepCount/(float)startCount*100.0);
    if(iPercent>lastPercent){
      cout << "\b\b\b" << iPercent;
      lastPercent=iPercent;
    }

  }
  cout << endl;

  if(out[0]!='\0'){
    FILE* f;
    f=fopen(out,"wt");

    //Heading line
	  fprintf(f,"File\tFirst Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\t");
	  fprintf(f,"Best Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications\n");

    for(i=0;i<vPeps.size();i++){
		  fprintf(f,"%s\t%d\t%d\t%d\t%d\t%lf\t%lf\t%f\t%f\t%f\t%f\t%f\t%lf\t%s\n","NULL",
																																		   vPeps[i].lowScan,
																																		   vPeps[i].highScan,
                                                                       vPeps[i].datapoints,
																																		   vPeps[i].charge,
																																		   vPeps[i].monoMass,
																																		   vPeps[i].basePeak,
																																		   vPeps[i].intensity,
																																		   vPeps[i].sumIntensity,
																																		   vPeps[i].firstRTime,
																																		   vPeps[i].lastRTime,
																																		   vPeps[i].rTime,
																																		   vPeps[i].xCorr,
																																		   vPeps[i].mods);
	  }
    
	  fclose(f);
  }

  return true;
}



bool CKronik2::findMax(vector<sScan>& v, int& s, int& p){
  bool found = false;
  float max=0;
  unsigned int i;
  p=0;
  for(i=0;i<v.size();i++){
    if(v[i].vPep->size()==0) continue;
    if(v[i].vPep->at(0).intensity>max){
      max=v[i].vPep->at(0).intensity;
      s=i;
      found = true;
    }
  }
  return found;
}


//-----------------------------------------  
//                 Tools
//-----------------------------------------
bool CKronik2::getRT(int scanNum, float& rt){

  unsigned int i;
  for(i=0;i<hkData.size();i++){
    if(hkData[i].scanNum==scanNum){
      rt=hkData[i].rTime;
      return true;
    }
  }

  rt=0.0f;
  return false;
}

//-----------------------------------------  
//                 Filters
//-----------------------------------------
void CKronik2::filterRT(float rt1, float rt2){
  unsigned int i;
  vector<sPepProfile> temp;

  for(i=0;i<vPeps.size();i++){
    if(vPeps[i].rTime>=rt1 && vPeps[i].rTime<=rt2){
      temp.push_back(vPeps[i]);
    }
  }

  vPeps.clear();
  for(i=0;i<temp.size();i++) vPeps.push_back(temp[i]);
}

void CKronik2::removeContaminants(float rt){
  unsigned int i;
  vector<sPepProfile> temp;

  for(i=0;i<vPeps.size();i++){
    if((vPeps[i].lastRTime-vPeps[i].firstRTime)<rt){
      temp.push_back(vPeps[i]);
    }
  }

  vPeps.clear();
  for(i=0;i<temp.size();i++) vPeps.push_back(temp[i]);
}

void CKronik2::removeMass(double m1, double m2){
	unsigned int i;
  vector<sPepProfile> temp;

  for(i=0;i<vPeps.size();i++){
		if( vPeps[i].monoMass>=m1 && vPeps[i].monoMass<=m2 ){
			temp.push_back(vPeps[i]);
    }
  }

  vPeps.clear();
  for(i=0;i<temp.size();i++) vPeps.push_back(temp[i]);
}

//-----------------------------------------
//    Parameter accessors and modifiers
//-----------------------------------------
void CKronik2::setPPMTol(double d){
  dPPMTol=d;
}

void CKronik2::setMatchTol(int i){
  iMatchTol=i;
}

void CKronik2::setGapTol(int i){
  iGapTol=i;
}

unsigned int CKronik2::size(){
  return vPeps.size();
}


//------------------------------------------------------
// STATISTICS FUNCTIONS (adapted from numerical recipes)
//------------------------------------------------------
bool CKronik2::pearson(int i1, int i2, bool byScan, bool interpolate, double& rval, double& pval, double& slope, double& intercept, float& sum1, float& sum2){

  int i,j;
  double sxx=0.0,syy=0.0,sxy=0.0;
  double ax=0.0,ay=0.0;
  double sumx=0.0,sumy=0.0,sumxy=0.0,sumxx=0.0;
  double yt,xt,t,df;
  double TINY=1.0e-20;
  double lrval;

  pval=1.0;
  rval=0.0;

  vector<float> p1;
  vector<float> p2;

  if(byScan){
    i=0;
    j=0;
    while(true){
      if(i==vPeps[i1].datapoints) return false;
      if(j==vPeps[i2].datapoints) return false;
      if(vPeps[i1].profile[i].scanNum<vPeps[i2].profile[j].scanNum){
        i++;
        continue;
      } else if(vPeps[i1].profile[i].scanNum>vPeps[i2].profile[j].scanNum){
        j++;
        continue;
      } else {
        break;
      }
    }
    while(i<vPeps[i1].datapoints && j<vPeps[i2].datapoints){
      //Build correlation vectors
      p1.push_back(vPeps[i1].profile[i++].intensity);
      p2.push_back(vPeps[i2].profile[j++].intensity);
    }
  } else {

    //Find center for each peptide
    float max;
    int max1, max2;

    max=0.0f;
    for(i=0;i<vPeps[i1].datapoints;i++){
      if(vPeps[i1].profile[i].intensity>max){
        max=vPeps[i1].profile[i].intensity;
        max1=i;
      }
    }

    max=0.0f;
    for(i=0;i<vPeps[i2].datapoints;i++){
      if(vPeps[i2].profile[i].intensity>max){
        max=vPeps[i2].profile[i].intensity;
        max2=i;
      }
    }

    //Build correlation vectors
    p1.push_back(vPeps[i1].profile[max1].intensity);
    p2.push_back(vPeps[i2].profile[max2].intensity);

    j=0;
    while((max1-j)>0 && (max2-j)>0){
      j++;
      if(vPeps[i1].profile[max1-j].interpolated && !interpolate) p1.push_back(0.0f);
      else p1.push_back(vPeps[i1].profile[max1-j].intensity);
      if(vPeps[i2].profile[max2-j].interpolated && !interpolate) p2.push_back(0.0f);
      else p2.push_back(vPeps[i2].profile[max2-j].intensity);
    }

    j=1;
    while((max1+j)<vPeps[i1].datapoints && (max2+j)<vPeps[i2].datapoints ){
      if(vPeps[i1].profile[max1+j].interpolated && !interpolate) p1.push_back(0.0f);
      else p1.push_back(vPeps[i1].profile[max1+j].intensity);
      if(vPeps[i2].profile[max2+j].interpolated && !interpolate) p2.push_back(0.0f);
      else p2.push_back(vPeps[i2].profile[max2+j].intensity);
      j++;
    }
  }

  //Get totals
  for(i=0;i<p1.size();i++){
    sumx+=p1[i];
    sumy+=p2[i];
    sumxy+=p1[i]*p2[i];
    sumxx+=p1[i]*p1[i];
  }

  //Pearsons
  ax=sumx/p1.size();
  ay=sumy/p2.size();
  for(i=0;i<p1.size();i++){
    xt=p1[i]-ax;
    yt=p2[i]-ay;
    sxx+=xt*xt;
    syy+=yt*yt;
    sxy+=xt*yt;
  }
  rval=sxy/(sqrt(sxx*syy)+TINY);
  if(rval<0.0) lrval=0.0;  //negative correlations are false hits in this analysis
  else lrval=rval; 

  df=p1.size()-2;
  if(df<1) {
    pval=1.0;
  } else {
    t=lrval*sqrt(df/((1-lrval+TINY)*(1+lrval+TINY)));
    pval=betai(0.5*df,0.5,df/(df+t*t));
  }

  //Line of best fit
  slope = (sumx*sumy - p1.size()*sumxy) / (sumx*sumx - p1.size()*sumxx);
  intercept = (sumx*sumxy - sumy*sumxx) / (sumx*sumx - p1.size()*sumxx);

  sum1=(float)sumx;
  sum2=(float)sumy;

  return true;

}

double CKronik2::betai(double a, double b, double x){
  double bt;
  if(x<0.0 || x>1.0) return 0.0;
  if(x==0.0 || x==1.0) {
    bt=0.0;
  } else {
    bt=exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(1.0-x));
  }
  if(x<(a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
  else return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double CKronik2::betacf(double a, double b, double x){
  int MAXIT = 100;
  double EPS = 3.0e-7;
  double FPMIN = 1.0e-30;

  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;

  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if(fabs(d)<FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for(m=1;m<=MAXIT;m++){
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if(fabs(d)<FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c)<FPMIN) c=FPMIN;
    d=1.0/d;
    h*=d*c;
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if(fabs(d)<FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c)<FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h*=del;
    if(fabs(del-1.0)<EPS) break;
  }
  if(m>MAXIT) return 0.0;
  return h;
}

double CKronik2::gammaln(double xx){
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp-=(x+0.5)*log(tmp);
  ser=1.000000000190015;
  for(j=0;j<=5;j++) ser+=cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

double CKronik2::interpolate(int x1, int x2, double y1, double y2, int x){
  double m;
  double b;

  m=(y2-y1)/(x2-x1);
  b=y2-m*x2;
  return m*x+b;

}

//--------------------------------
//  Sorting Functions
//--------------------------------
void CKronik2::sortBasePeak(){
	qsort(&vPeps[0],vPeps.size(),sizeof(sPepProfile),compareBP);
}

void CKronik2::sortMonoMass(){
	qsort(&vPeps[0],vPeps.size(),sizeof(sPepProfile),compareMM);
}

void CKronik2::sortFirstRTime(){
	qsort(&vPeps[0],vPeps.size(),sizeof(sPepProfile),compareFRT);
}

void CKronik2::sortIntensityRev(){
	qsort(&vPeps[0],vPeps.size(),sizeof(sPepProfile),compareIRev);
}

int CKronik2::compareBP(const void *p1, const void *p2){
  const sPepProfile d1 = *(sPepProfile *)p1;
  const sPepProfile d2 = *(sPepProfile *)p2;
  if(d1.basePeak<d2.basePeak) return -1;
  else if(d1.basePeak>d2.basePeak) return 1;
  else return 0;
}

int CKronik2::compareMM(const void *p1, const void *p2){
  const sPepProfile d1 = *(sPepProfile *)p1;
  const sPepProfile d2 = *(sPepProfile *)p2;
  if(d1.monoMass<d2.monoMass) return -1;
  else if(d1.monoMass>d2.monoMass) return 1;
  else return 0;
}

int CKronik2::compareFRT(const void *p1, const void *p2){
  const sPepProfile d1 = *(sPepProfile *)p1;
  const sPepProfile d2 = *(sPepProfile *)p2;
  if(d1.firstRTime<d2.firstRTime) return -1;
  else if(d1.firstRTime>d2.firstRTime) return 1;
  else return 0;
}

int CKronik2::compareIRev(const void *p1, const void *p2){
  const sPepProfile d1 = *(sPepProfile *)p1;
  const sPepProfile d2 = *(sPepProfile *)p2;
  //cout << d1.intensity << " " << d2.intensity << endl;
  if(d1.intensity>d2.intensity) return -1;
  else if(d1.intensity<d2.intensity) return 1;
  else return 0;
}
