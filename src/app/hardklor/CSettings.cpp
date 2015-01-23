#include "CSettings.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string.h>

using namespace std;

CSettings::CSettings(){
  chlor=0;
  selen=0;
  maxPep=5;
  maxCharge=0;
  maxDepth=100;
  smooth=0;
  signal=0;
  scan=0;
  S2N=0;
  windowLower=0;
  windowUpper=0;
  enrichAtom[0]='\0';
  enrich=0;
  enrichTimes=0;
  QAR = false;
  ROC = false;
  readFile((char*)"settings.conf");
};

CSettings::CSettings(char *c){
  chlor=0;
  selen=0;
  maxPep=5;
  maxCharge=0;
  maxDepth=100;
  smooth=0;
  signal=0;
  scan=0;
  S2N=0;
  windowLower=0;
  windowUpper=0;
  enrichAtom[0]='\0';
  enrich=0;
  enrichTimes=0;
  QAR = false;
  ROC = false;
  readFile(c);
};

void CSettings::readFile(char* c) {
  fstream fptr;
  FileName* f;
  char tstr[256];
  char *tok;
  int i;
  varType v;

  fptr.open(c,ios::in);
  if(!fptr.good()){
    cout << "Cannot read settings file: " << c << endl;
    fptr.close();
    return;
  };

  while(!fptr.eof()){
    fptr.getline(tstr,100);
    if(tstr[0]=='#' || tstr[0]==0) continue;

    tok = strtok(tstr," \n\t\r");
    if(tok==0) continue;
    if(tok[0]=='#') continue;

    if(strcmp(tok,"[CHLORINE]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      chlor = atoi(tok);
      for(i=0;i<chlor;i++){
	v.code=1;
	v.value=i+1;
	variants.push_back(v);
      };
    } else if(strcmp(tok,"[SELENIUM]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      selen = atoi(tok);
      for(i=0;i<selen;i++){
	v.code=2;
	v.value=i+1;
	variants.push_back(v);
      };
    } else if(strcmp(tok,"[ENRICH]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      strcpy(enrichAtom,tok);
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      enrich = atof(tok);
      if(strcmp(enrichAtom,"N")==0){
	v.code=3;
	v.value=enrich;
	variants.push_back(v);
      };
      enrichTimes=1;
    } else if(strcmp(tok,"[SMOOTH]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      smooth = atoi(tok);
    } else if(strcmp(tok,"[CORR]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      corrThresh = atof(tok);
    } else if(strcmp(tok,"[S2N]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      S2N = atof(tok);
    } else if(strcmp(tok,"[WINDOW]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      windowLower = atof(tok);
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      windowUpper = atof(tok);
    } else if(strcmp(tok,"[MAXPEP]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      maxPep = atoi(tok);
    } else if(strcmp(tok,"[MAXCHARGE]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      maxCharge = atoi(tok);
    } else if(strcmp(tok,"[MAXDEPTH]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      maxDepth = atoi(tok);
    } else if(strcmp(tok,"[QAR]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      if(atoi(tok)!=0) QAR = true;
    } else if(strcmp(tok,"[ROC]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      if(atoi(tok)!=0) ROC = true;
    } else if(strcmp(tok,"[SCAN]")==0) {
      tok = strtok(NULL, " \n\t\r");
      if(tok==NULL) continue;
      scan = atoi(tok);
    } else if(strcmp(tok,"[FILE]")==0){
      f = new FileName;
      tok = strtok(NULL,"\t\n\r");
      if(tok==NULL) continue;
      if(strcmp(tok,"Zoom")==0) f->st=Zoom;
      else if(strcmp(tok,"UltraZoom")==0) f->st=UltraZoom;
      else if(strcmp(tok,"IonSpec")==0) f->st=IonSpec2;
      else if(strcmp(tok,"Other")==0) f->st=Other;
      tok = strtok(NULL,"\t\n\r");
      if(tok==NULL) continue;
      else strcpy(f->infile,tok);
      tok = strtok(NULL,"\t\n\r");
      if(tok==NULL) continue;
      else strcpy(f->outfile,tok);
      files.push_back(f);
    } else if(strcmp(tok,"[SIGNAL]")==0){
      tok = strtok(NULL," \t\n\r");
      if(tok==NULL) continue;
      signal = atof(tok);
    };
  };

};

  
