#include "CHardklorParser.h"
#include <stdlib.h>
#include <string.h>

using namespace std;

CHardklorParser::CHardklorParser(){
  vQueue = new vector<CHardklorSetting>;
}

CHardklorParser::~CHardklorParser(){
  delete vQueue;
}

//Takes a command line and breaks it into tokens
//Tokens are then read and used to set global and local parameters
void CHardklorParser::parse(char* cmd) {
  bool isGlobal=true;
  sMolecule m;
  char *tok;
  char tmpstr[256];
	char upStr[64];
  unsigned int i;
  string tstr;
  vector<string> vs;

	//For modifications
	CHardklorVariant v;
	CPeriodicTable* PT;
	int j,k;
	string atom;
	string isotope;
	string percent;
	bool badMod = false;
	int atomNum;
	bool bNew;

  CHardklorSetting hs;

	//Replace first # with a terminator
	tok=strstr(cmd,"#");
	if(tok!=NULL) strncpy(tok,"\0",1);

  //Replace the space+hyphen with a delimiter
  while(true){
    tok = strstr(cmd," -");
    if(tok==NULL) break;
    strncpy(tok,"\n",1);
  }

  //Break parameter string into tokens
  tok = strtok(cmd,"\n");
  while(tok!=NULL){
    vs.push_back(tok);
    tok = strtok(NULL,"\n");
  }

  //Analyze each token
  for(i=0;i<vs.size();i++){

    //Grab first subtoken, which should be the parameter to change
    tok = strtok(&vs.at(i)[0]," \t\n");
    if(tok==NULL) continue;

    //-nb : No Base Molecule
    //followed by true/false. If not followed by anything, assumed true.
    //invalid characters assumed false
    if(strcmp(tok,"-nb")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        if(isGlobal) global.noBase=true;
			  else hs.noBase=true;
        continue;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.noBase=true;
			  else hs.noBase=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.noBase=false;
			  else hs.noBase=false;
			} else {
				if(isGlobal) global.noBase=true;
			  else hs.noBase=true;
      }

    //-ns : Do not split the spectrum
    } else if(strcmp(tok,"-ns")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        if(isGlobal) global.noSplit=true;
			  else hs.noSplit=true;
        continue;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.noSplit=true;
			  else hs.noSplit=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.noSplit=false;
			  else hs.noSplit=false;
			} else {
				if(isGlobal) global.noSplit=true;
			  else hs.noSplit=true;
      }

    //-i : Intersection analysis
    } else if(strcmp(tok,"-i")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        if(isGlobal) global.iAnalysis=true;
			  else hs.iAnalysis=true;
        continue;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.iAnalysis=true;
			  else hs.iAnalysis=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.iAnalysis=false;
			  else hs.iAnalysis=false;
			} else {
				if(isGlobal) global.iAnalysis=true;
			  else hs.iAnalysis=true;
      }

    //-da : Output distribution area instead of base peak intensity
    } else if(strcmp(tok,"-da")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        if(isGlobal) global.distArea=true;
			  else hs.distArea=true;
        continue;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.distArea=true;
			  else hs.distArea=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.distArea=false;
			  else hs.distArea=false;
			} else {
				if(isGlobal) global.distArea=true;
			  else hs.distArea=true;
      }

    //-da : Output distribution area instead of base peak intensity
    } else if(strcmp(tok,"-xml")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        if(isGlobal) global.xml=true;
			  else hs.xml=true;
        continue;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.xml=true;
			  else hs.xml=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.xml=false;
			  else hs.xml=false;
			} else {
				if(isGlobal) global.xml=true;
			  else hs.xml=true;
      }

    //-zero : Zero intensity analysis - allow zero intensity values
    } else if(strcmp(tok,"-zero")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        if(isGlobal) global.skipZero=false;
			  else hs.skipZero=false;
        continue;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.skipZero=false;
			  else hs.skipZero=false;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.skipZero=true;
			  else hs.skipZero=true;
			} else {
				if(isGlobal) global.skipZero=false;
			  else hs.skipZero=false;
      }
    
    //-d : Depth of recursion
    } else if(strcmp(tok,"-d")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atoi(tok)>0) {
				if(isGlobal) global.depth=atoi(tok);
				else hs.depth=atoi(tok);
      }

    //-corr : Correlation Threshold
    } else if(strcmp(tok,"-corr")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atof(tok)>0) {
				if(isGlobal) global.corr=atof(tok);
				else hs.corr=atof(tok);
      };

    //-c : Assume Centroid Data
    } else if(strcmp(tok,"-c")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        if(isGlobal) global.centroid=true;
				else hs.centroid=true;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.centroid=true;
				else hs.centroid=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.centroid=false;
				else hs.centroid=false;
			} else {
				if(isGlobal) global.centroid=true;
				else hs.centroid=true;
      }

	  //-mol :  Use specific molecular formula
    } else if(strcmp(tok,"-mol")==0) {
      tok = strtok(NULL,"\n");
      if(tok==NULL) continue;
			if(strlen(tok)>64) {
				cout << "User specified -mol formula is too long! Max = 64 characters." << endl;
			} else {
				if(isGlobal) strcpy(global.formula,tok);
				else strcpy(hs.formula,tok);
			}

    //-mdat : Path and name of mercury data file
    } else if(strcmp(tok,"-mdat")==0) {
      tok = strtok(NULL,"\n");
      if(tok==NULL) continue;
      if(isGlobal) strcpy(global.MercuryFile,tok);
			else strcpy(hs.MercuryFile,tok);

    //-hdat :  Path and name of hardklor data file
    } else if(strcmp(tok,"-hdat")==0) {
      tok = strtok(NULL,"\n");
      if(tok==NULL) continue;
			if(isGlobal) strcpy(global.HardklorFile,tok);
			else strcpy(hs.HardklorFile,tok);

    //-m : Parse averagine molecule variants
    } else if(strcmp(tok,"-m")==0) {
			v.clear();
			PT = new CPeriodicTable(global.HardklorFile);

			while(true){
				tok = strtok(NULL," \t\n");
				if(tok==NULL) break;

				badMod=false;
				if(isdigit(tok[0]) || tok[0]=='.'){
					//we have enrichment
					percent="";
					atom="";
					isotope="";

					//Get the APE
					for(j=0;j<(int)strlen(tok);j++){
						if(isdigit(tok[j]) || tok[j]=='.') {
							if(percent.size()==15){
								cout << "Malformed modification flag: too many digits." << endl;
								badMod=true;
								break;
							}
							percent+=tok[j];
						} else {
							break;
						}
					}
					if(badMod) break;

					//Get the atom
					for(j=j;j<(int)strlen(tok);j++){
						if(isalpha(tok[j])) {
							if(atom.size()==2){
								cout << "Malformed modification flag: invalid atom" << endl;
								badMod=true;
								break;
							}
							atom+=tok[j];
						} else {
							break;
						}
					}
					if(badMod) break;

					//Get the isotope
					for(j=j;j<(int)strlen(tok);j++){
						if(isotope.size()==2){
							cout << "Malformed modification flag: bad isotope" << endl;
							badMod=true;
							break;
						}
						isotope+=tok[j];
					}
					if(badMod) break;

					//format the atom properly
					atom.at(0)=toupper(atom.at(0));
				  if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

					//Get the array number for the atom
					atomNum=-1;
					for(j=0;j<PT->size();j++){
						if(strcmp(PT->at(j).symbol,&atom[0])==0){
							atomNum=j;
							break;
						}
					}

					if(atomNum==-1){
						cout << "Malformed modification flag: Atom not in periodic table: " << atom << endl;
						break;
					}

					v.addEnrich(atomNum,atoi(&isotope[0]),atof(&percent[0]));

				} else {
					//we have molecule
					percent="";
					atom="";
					bNew=true;

					//Get the atom
					for(j=0;j<(int)strlen(tok);j++){
      
						//Check for an atom symbol
						if(isalpha(tok[j])) {
	
							//Check if we were working on the count of the previous atom
							if(!bNew) {
								bNew=true;
	  
								//Make sure the atom has uppercase-lowercase letter format;
								atom.at(0)=toupper(atom.at(0));
								if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

								//Look up the new atom
								for(k=0;k<PT->size();k++){
									if(strcmp(PT->at(k).symbol,&atom[0])==0){
										//Add the new atom to the variant
										v.addAtom(k,atoi(&percent[0]));
										break;
									}
								}
	  
								//Clear the fields
								percent="";
								atom="";
							}
	
							//Add this letter to the atom symbol
							if(atom.size()==2){
								cout << "Malformed modification flag: invalid atom" << endl;
								badMod=true;
								break;
							}
							atom+=tok[j];

							if(badMod) break;
	
						} else if(isdigit(tok[j])){
	
							//Whenever we find a digit, we have already found an atom symbol
							bNew=false;
	
							//Add this letter to the atom count
							if(percent.size()==12){
								cout << "Malformed modification flag: unreasonable atom count" << endl;
								badMod=true;
								break;
							}
							percent+=tok[j];

							if(badMod) break;
	
						}      
					}

					//process the last atom
					//Make sure the atom has uppercase-lowercase letter format;
					atom.at(0)=toupper(atom.at(0));
					if(atom.size()==2) atom.at(1)=tolower(atom.at(1));

			    //Look up the new atom
					for(k=0;k<PT->size();k++){
						if(strcmp(PT->at(k).symbol,&atom[0])==0){

							//Add the new atom to the variant
							v.addAtom(k,atoi(&percent[0]));
							break;
	
						}
					}
				}

			}

			if(badMod) continue;
			if(isGlobal) global.variant->push_back(v);
			else hs.variant->push_back(v);

			delete PT;
      
    //-s : Savitsky-Golay Smoothing
    } else if(strcmp(tok,"-s")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atoi(tok)>0) {
				if(isGlobal) global.smooth=atoi(tok);
				else hs.smooth=atoi(tok);
      }
	  
    //-sc : Read specific scan number(s)
		} else if(strcmp(tok,"-sc")==0) {
			tok = strtok(NULL," \t\n");
			if(tok==NULL) continue;
			if(atoi(tok)>0) {
				if(isGlobal) global.scan.iLower=atoi(tok);
				else hs.scan.iLower=atoi(tok);
			}
			tok = strtok(NULL," \t\n");
			if(tok==NULL) continue;
			if(atoi(tok)>0) {
				if(isGlobal) global.scan.iUpper=atoi(tok);
				else hs.scan.iUpper=atoi(tok);
			}
	 
    //-w : Spectrum window boundaries
    } else if(strcmp(tok,"-w")==0) {
			tok = strtok(NULL," \t\n");
			if(tok==NULL) continue;
			if(atof(tok)>=1.0) {
				if(isGlobal) global.window.dLower=atof(tok);
				else hs.window.dLower=atof(tok);
			}
			tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atof(tok)>=1.0) {
				if(isGlobal) global.window.dUpper=atof(tok);
				else hs.window.dUpper=atof(tok);
      }

    //-sl : sensitivity level
    } else if(strcmp(tok,"-sl")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atof(tok)>=0) {
				if(isGlobal) global.sl=atoi(tok);
				else hs.sl=atoi(tok);
      }

    //-sna : Signal to Noise Algorithm
    } else if(strcmp(tok,"-sna")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "Unspecified signal to noise algorithm. Defaulting to STD method." << endl;
        continue;
      }
      if(strcmp(tok,"STD")==0 || strcmp(tok,"std")==0) {
				if(isGlobal) global.sna=0;
				else hs.sna=0;
			} else if (strcmp(tok,"ASN")==0 || strcmp(tok,"asn")==0) {
				if(isGlobal) global.sna=1;
				else hs.sna=1;
      } else if (strcmp(tok,"AVG")==0 || strcmp(tok,"avg")==0) {
				if(isGlobal) global.sna=2;
				else hs.sna=2;
      } else if (strcmp(tok,"MIX")==0 || strcmp(tok,"mix")==0) {
				if(isGlobal) global.sna=3;
				else hs.sna=3;
			} else if (strcmp(tok,"V2")==0 || strcmp(tok,"v2")==0) {
				if(isGlobal) global.sna=4;
				else hs.sna=4;
			} else if (strcmp(tok,"V2ASN")==0 || strcmp(tok,"v2asn")==0) {
				if(isGlobal) global.sna=5;
				else hs.sna=5;
			} else if (strcmp(tok,"V2MIX")==0 || strcmp(tok,"v2mix")==0) {
				if(isGlobal) global.sna=6;
				else hs.sna=6;
			} else {
				cout << "Unkown signal to noise algorithm: " << tok << ". Defaulting to STD method." << endl;
				if(isGlobal) global.sna=0;
				else hs.sna=0;
			}

    //-ppMatch : Match threshold for preprocessing
    } else if(strcmp(tok,"-ppMatch")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "Unspecified ppMatch parameter. Using default value: 1" << endl;
        continue;
      }
      if(atoi(tok)>=0) {
        if(isGlobal) global.ppMatch=atoi(tok);
        else hs.ppMatch=atoi(tok);
      } else {
        cout << "Invalid ppMatch parameter. Using default value: 1" << endl;
      }

    //-ppm : PPM threshold for preprocessing
    } else if(strcmp(tok,"-ppm")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "Unspecified -ppm parameter. Using default value: 10.0" << endl;
        continue;
      }
      if(atof(tok)>=0) {
        if(isGlobal) global.ppm=atof(tok);
        else hs.ppm=atof(tok);
      } else {
        cout << "Invalid -ppm parameter. Using default value: 10.0" << endl;
      }

    //-ppWin : Window size for preprocessing
    } else if(strcmp(tok,"-ppWin")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "Unspecified pWin parameter. Using default value: 1" << endl;
        continue;
      }
      if(atoi(tok)>=1) {
        if(isGlobal) global.ppWin=atoi(tok);
        else hs.ppWin=atoi(tok);
      } else {
        cout << "Invalid ppWin parameter. Using default value: 1" << endl;
      }

    } else if(strcmp(tok,"-rf")==0) {
      tok = strtok(NULL,"\n");
      if(tok==NULL) {
        cout << "WARNING: Unspecified rf parameter. No filter applied." << endl;
        continue;
      }
      if(isGlobal) strcpy(global.rawFilter,tok);
      else strcpy(hs.rawFilter,tok);

    //-avg : Average RAW scan data
      /*
    } else if(strcmp(tok,"-avg")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "WARNING: Unknown value for -avg. Setting to default value: false" << endl;
        if(isGlobal) global.rawAvg=false;
				else hs.rawAvg=false;
      }
			if(strcmp(tok,"true")==0){
				if(isGlobal) global.rawAvg=true;
				else hs.rawAvg=true;
			} else if(strcmp(tok,"false")==0 || atoi(tok)==0) {
				if(isGlobal) global.rawAvg=false;
				else hs.rawAvg=false;
			} else {
        cout << "WARNING: Unknown value for -avg. Setting to default value: false" << endl;
				if(isGlobal) global.rawAvg=false;
				else hs.rawAvg=false;
      }
      */

    //-avgWin : Window size for RAW scan averaging (+/-)
      /*
    } else if(strcmp(tok,"-avgWin")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "WARNING: Unknown value for -avgWin. Using default value: 1" << endl;
        if(isGlobal) global.rawAvgWidth=1;
        else hs.rawAvgWidth=1;
        continue;
      }
      if(atoi(tok)>=1) {
        if(isGlobal) global.rawAvgWidth=atoi(tok);
        else hs.rawAvgWidth=atoi(tok);
      } else {
        cout << "WARNING: Unknown value for -avgWin. Using default value: 1" << endl;
        if(isGlobal) global.rawAvgWidth=1;
        else hs.rawAvgWidth=1;
      }
      */

    //-avgCut : Intensity cutoff for RAW scan averaging
      /*
    } else if(strcmp(tok,"-avgCut")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) {
        cout << "WARNING: Unknown value for -avgCut. Using default value: 1000" << endl;
        if(isGlobal) global.rawAvgCutoff=1000;
        else hs.rawAvgCutoff=1000;
        continue;
      }
      if(atoi(tok)>=1) {
        if(isGlobal) global.rawAvgCutoff=atoi(tok);
        else hs.rawAvgCutoff=atoi(tok);
      } else {
        cout << "WARNING: Unknown value for -avgCut. Using default value: 1000" << endl;
        if(isGlobal) global.rawAvgCutoff=1000;
        else hs.rawAvgCutoff=1000;
      }
      */
      
    //-sn : Signal to Noise Ratio
    } else if(strcmp(tok,"-sn")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atof(tok)>=0) {
				if(isGlobal) global.sn=atof(tok);
				else hs.sn=atof(tok);
      }

    //-snWin : Signal to Noise Window Size
    } else if(strcmp(tok,"-snWin")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atof(tok)>0) {
				if(isGlobal) global.snWindow=atof(tok);
				else hs.snWindow=atof(tok);
      }
			tok = strtok(NULL," \t\n");
			if(tok==NULL) {
				if(isGlobal) global.staticSN=false;
				else hs.staticSN=false;
				continue;
			}
      if(strcmp(tok,"true")==0 || atoi(tok)>0) {
				if(isGlobal) global.staticSN=true;
				else hs.staticSN=true;
      }

    //-win : Window size when dividing spectrum
    } else if(strcmp(tok,"-win")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atof(tok)>0) {
				if(isGlobal) global.winSize=atof(tok);
				else hs.winSize=atof(tok);
      }

    //-chMax : Maximum charge
    } else if(strcmp(tok,"-chMax")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atoi(tok)>0) {
				if(isGlobal) global.maxCharge=atoi(tok);
				else hs.maxCharge=atoi(tok);
      }

    //-chMin : Minimum charge
    } else if(strcmp(tok,"-chMin")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atoi(tok)>0) {
				if(isGlobal) global.minCharge=atoi(tok);
				else hs.minCharge=atoi(tok);
      }

    //-p : Maximum number of peptide models
    } else if(strcmp(tok,"-p")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atoi(tok)>0) {
				if(isGlobal) global.peptide=atoi(tok);
				else hs.peptide=atoi(tok);
      }

    //-res : Resolution
    } else if(strcmp(tok,"-res")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(atoi(tok)>0) {
				if(isGlobal) global.res400=atof(tok);
				else hs.res400=atof(tok);
      }
			tok = strtok(NULL," \t\n");
			for(j=0;j<(int)strlen(tok);j++) upStr[j]=toupper(tok[j]);
			upStr[j]='\0';
			if(strcmp(upStr,"FTICR")==0) {
				if(isGlobal) global.msType=FTICR;
				else hs.msType=FTICR;
			} else if (strcmp(upStr,"ORBITRAP")==0) {
				if(isGlobal) global.msType=OrbiTrap;
				else hs.msType=OrbiTrap;
			} else if (strcmp(upStr,"TOF")==0) {
				if(isGlobal) global.msType=TOF;
				else hs.msType=TOF;
			} else if (strcmp(upStr,"QIT")==0) {
				if(isGlobal) global.msType=QIT;
				else hs.msType=QIT;
			} else {
				if(isGlobal) global.msType=FTICR;
				else hs.msType=FTICR;
			}

    //-mF : mzXML Filter
    } else if(strcmp(tok,"-mF")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(strcmp(tok,"MS1")==0) {
				if(isGlobal) global.mzXMLFilter=MS1;
				else hs.mzXMLFilter=MS1;
      } else if(strcmp(tok,"MS2")==0){
				if(isGlobal) global.mzXMLFilter=MS2;
				else hs.mzXMLFilter=MS2;
			} else if(strcmp(tok,"MS3")==0){
				if(isGlobal) global.mzXMLFilter=MS3;
				else hs.mzXMLFilter=MS3;
			} else if(strcmp(tok,"ZS")==0){
				if(isGlobal) global.mzXMLFilter=ZS;
				else hs.mzXMLFilter=ZS;
			} else if(strcmp(tok,"UZS")==0){
				if(isGlobal) global.mzXMLFilter=UZS;
				else hs.mzXMLFilter=UZS;
      }

		//-a : Algorithm
    } else if(strcmp(tok,"-a")==0) {
			tok = strtok(NULL," \t\n");
			if(strcmp(tok,"Basic")==0) {
				if(isGlobal) global.algorithm=Basic;
				else hs.algorithm=Basic;
			} else if (strcmp(tok,"SemiComplete")==0) {
				if(isGlobal) global.algorithm=SemiComplete;
				else hs.algorithm=SemiComplete;
			} else if (strcmp(tok,"SemiCompleteFast")==0) {
				if(isGlobal) global.algorithm=SemiCompleteFast;
				else hs.algorithm=SemiCompleteFast;
			} else if (strcmp(tok,"Dynamic")==0) {
				if(isGlobal) global.algorithm=Dynamic;
				else hs.algorithm=Dynamic;
			} else if (strcmp(tok,"DynamicSemiComplete")==0) {
				if(isGlobal) global.algorithm=DynamicSemiComplete;
				else hs.algorithm=DynamicSemiComplete;
			} else if (strcmp(tok,"SemiSubtractive")==0) {
				if(isGlobal) global.algorithm=SemiSubtractive;
				else hs.algorithm=SemiSubtractive;
			} else if (strcmp(tok,"FewestPeptides")==0) {
				if(isGlobal) global.algorithm=FewestPeptides;
				else hs.algorithm=FewestPeptides;
			} else if (strcmp(tok,"FewestPeptidesChoice")==0) {
				if(isGlobal) global.algorithm=FewestPeptidesChoice;
				else hs.algorithm=FewestPeptidesChoice;
			} else if (strcmp(tok,"FastFewestPeptides")==0) {
				if(isGlobal) global.algorithm=FastFewestPeptides;
				else hs.algorithm=FastFewestPeptides;
			} else if (strcmp(tok,"FastFewestPeptidesChoice")==0) {
				if(isGlobal) global.algorithm=FastFewestPeptidesChoice;
				else hs.algorithm=FastFewestPeptidesChoice;
			} else {
				if(isGlobal) global.algorithm=Basic;
				else hs.algorithm=Basic;
			}

    //-cdm : Charge Determination Mode
    } else if(strcmp(tok,"-cdm")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
			switch(tok[0]){
				case 'C':
					if(isGlobal) global.chargeMode='C';
					else hs.chargeMode='C';
					break;
				case 'P':
					if(isGlobal) global.chargeMode='P';
					else hs.chargeMode='P';
					break;
				case 'F':
					if(isGlobal) global.chargeMode='F';
					else hs.chargeMode='F';
					break;
				case 'S':
					if(isGlobal) global.chargeMode='S';
					else hs.chargeMode='S';
					break;
				case 'Q':
				default:
					if(isGlobal) global.chargeMode='Q';
					else hs.chargeMode='Q';
					break;
			}

    /*
    //-o : Output File (deprecated)
    } else if(strcmp(tok,"-o")==0) {
      tok = strtok(NULL," \t\n");
      if(tok==NULL) continue;
      if(isGlobal) strcpy(global.outFile,tok);
      else strcpy(hs.outFile,tok);
    */
      
    //Reading file or invalid parameter
    } else {
      if(tok[0]=='-') {
	      //we have invalid parameter
        cout << "WARNING: Unknown parameter " << tok << endl;
      } else {
        isGlobal=false;
				hs = global;
        
        //on systems that allow a space in the path, require quotes (") to capture
        //complete file name
        strcpy(tmpstr,tok);

        //Check for quote
        if(tok[0]=='\"'){

          if(tok[strlen(tok)-1]=='\"') {
            strcpy(tmpstr,tok);
          } else {

            //continue reading tokens until another quote is found
            while(true){
              tok=strtok(NULL," \t\n");
              if(tok==NULL) {
                cout << "Invalid input file." << endl;
                exit(-1);
              }
              strcat(tmpstr," ");
              strcat(tmpstr,tok);
              if(tok[strlen(tok)-1]=='\"') break;
            }
          }
          //For some reason, the null string terminator is not places, so place it manually
          strncpy(hs.inFile,&tmpstr[1],strlen(tmpstr)-2);
          hs.inFile[strlen(tmpstr)-2]='\0';

        //If no quote, assume file name has no spaces
        } else {
					strcpy(hs.inFile,tok);
        }

				tok = strtok(NULL," \t\n");
				if(tok==NULL) {
					cout << "Invalid output file." << endl;
					exit(-1);
				}

        //Repeat the entire process for output file
        strcpy(tmpstr,tok);
        if(tok[0]=='\"'){

          if(tok[strlen(tok)-1]=='\"') {
            strcpy(tmpstr,tok);
          } else {

            //continue reading tokens until another quote is found
            while(true){
              tok=strtok(NULL," \t\n");
              if(tok==NULL) {
                cout << "Invalid output file." << endl;
                exit(-1);
              }
              strcat(tmpstr," ");
              strcat(tmpstr,tok);
              if(tok[strlen(tok)-1]=='\"') break;
            }
          }
          //For some reason, the null string terminator is not places, so place it manually
          strncpy(hs.outFile,&tmpstr[1],strlen(tmpstr)-2);
          hs.outFile[strlen(tmpstr)-2]='\0';
          
        } else {
					strcpy(hs.outFile,tok);
        }

				hs.fileFormat = getFileFormat(hs.inFile);
      }
    }
    
  }
  
  //After all tokens are parsed, push back settings if file was parsed
  if(!isGlobal){
    isGlobal=true;
    vQueue->push_back(hs);
  }

}

//Reads in a config file and passes it to the parser
void CHardklorParser::parseConfig(char* c){
  fstream fptr;
  char tstr[512];

  fptr.open(c,ios::in);
  if(!fptr.good()){
    cout << "Cannot open config file!" << endl;
    return;
  }

  while(!fptr.eof()) {
    fptr.getline(tstr,512);
    if(tstr[0]==0) continue;
    if(tstr[0]=='#') continue;
    parse(tstr);
  }

  fptr.close();
}

CHardklorSetting& CHardklorParser::queue(int i){
  return vQueue->at(i);
}

int CHardklorParser::size(){
  return vQueue->size();
}

//Identifies file format from extension - Must conform to these conventions
MSFileFormat CHardklorParser::getFileFormat(char* c){

	char file[256];
	char ext[256];
	char *tok;

	strcpy(file,c);
	tok=strtok(file,".\n");
	while(tok!=NULL){
		strcpy(ext,tok);
		tok=strtok(NULL,".\n");
	}

	if(strcmp(ext,"ms1")==0 || strcmp(ext,"MS1")==0) return ms1;
	if(strcmp(ext,"ms2")==0 || strcmp(ext,"MS2")==0) return ms2;
	if(strcmp(ext,"bms1")==0 || strcmp(ext,"BMS1")==0) return bms1;
	if(strcmp(ext,"bms2")==0 || strcmp(ext,"BMS2")==0) return bms2;
	if(strcmp(ext,"cms1")==0 || strcmp(ext,"CMS1")==0) return cms1;
	if(strcmp(ext,"cms2")==0 || strcmp(ext,"CMS2")==0) return cms2;
	if(strcmp(ext,"zs")==0 || strcmp(ext,"ZS")==0) return zs;
	if(strcmp(ext,"uzs")==0 || strcmp(ext,"UZS")==0) return uzs;
	if(strcmp(ext,"mzML")==0 || strcmp(ext,"MZML")==0) return mzML;
	if(strcmp(ext,"mzXML")==0 || strcmp(ext,"MZXML")==0) return mzXML;
  if(strcmp(ext,"raw")==0 || strcmp(ext,"RAW")==0) return raw;
	return dunno;

}
