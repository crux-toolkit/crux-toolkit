#ifndef _STPETER_H
#define _STPETER_H

#include <algorithm>
#include <ctime>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

//Custom Data Structures
typedef struct stpDB {
  string name;
  string sequence;
} stpDB;

typedef struct stpMS2 {
  int     fileID;
  string  xLabel;
  int     scanNum;
  double  sumIntensity;
} stpMS2;

typedef struct stpPSMIndex {
  string  id;
  size_t  protID;
  size_t  pepID;
} stpPSMIndex;

/*
typedef struct stpPeptide {
  int             charge;
  double          mass;
  double          spectralIndex;
  string          sequence;
  vector<stpMS2>* scans;
  stpPeptide(){
    charge=0;
    mass=0;
    spectralIndex=0;
    sequence="";
    scans=new vector<stpMS2>;
  }
  stpPeptide(const stpPeptide& s){
    charge=s.charge;
    mass=s.mass;
    spectralIndex=s.spectralIndex;
    sequence=s.sequence;
    scans=new vector<stpMS2>;
    for(size_t i=0;i<s.scans->size();i++) scans->push_back(s.scans->at(i));
  }
  ~stpPeptide(){
    delete scans;
  }
  stpPeptide& operator=(const stpPeptide& s){
    if(this!=&s){
      charge=s.charge;
      mass=s.mass;
      spectralIndex = s.spectralIndex;
      sequence=s.sequence;
      delete scans;
      scans=new vector<stpMS2>;
      for(size_t i=0;i<s.scans->size();i++) scans->push_back(s.scans->at(i));
    }
    return *this;
  }
} stpPeptide;

typedef struct stpProtein{
  int                 length;
  double              logSIN;
  double              ng;
  double              SIN;
  string              name;
  string              description;
  vector<stpPeptide>* peptides;
  stpProtein(){
    length=0;
	  logSIN=0;
	  ng=0;
    SIN=0;
    name="";
    description="";
    peptides = new vector<stpPeptide>;
  }
  stpProtein(const stpProtein& p){
    length=p.length;
	  logSIN=p.logSIN;
	  ng=p.ng;
    SIN=p.SIN;
    name=p.name;
    description=p.description;
    peptides = new vector<stpPeptide>;
    for(size_t i=0;i<p.peptides->size();i++) peptides->push_back(p.peptides->at(i));
  }
  ~stpProtein(){
    delete peptides;
  }
  stpProtein& operator=(const stpProtein& p){
    if(this!=&p){
      length=p.length;
	    logSIN=p.logSIN;
	    ng=p.ng;
      SIN=p.SIN;
      name=p.name;
      description=p.description;
      delete peptides;
      peptides = new vector<stpPeptide>;
      for(size_t i=0;i<p.peptides->size();i++) peptides->push_back(p.peptides->at(i));
    }
    return *this;
  }
} stpProtein;
*/

typedef struct stpPeptide{
  int             charge;
  double          mass;
  double          spectralCount;
  double          spectralIndex;
  string          sequence;
  string          preciseSequence;
  vector<size_t>* protIndex;
  vector<stpMS2>* scans;
  stpPeptide(){
    charge = 0;
    mass = 0;
    spectralCount = 0;
    spectralIndex = 0;
    sequence = "";
    preciseSequence="";
    scans = new vector<stpMS2>;
    protIndex = new vector<size_t>;
  }
  stpPeptide(const stpPeptide& s){
    size_t i;
    charge = s.charge;
    mass = s.mass;
    spectralCount = s.spectralCount;
    spectralIndex = s.spectralIndex;
    sequence = s.sequence;
    preciseSequence = s.preciseSequence;
    scans = new vector<stpMS2>;
    protIndex = new vector<size_t>;
    for (i = 0; i<s.scans->size(); i++) scans->push_back(s.scans->at(i));
    for (i = 0; i<s.protIndex->size(); i++) protIndex->push_back(s.protIndex->at(i));
  }
  ~stpPeptide(){
    delete protIndex;
    delete scans;
  }
  stpPeptide& operator=(const stpPeptide& s){
    if (this != &s){
      size_t i;
      charge = s.charge;
      mass = s.mass;
      spectralCount = s.spectralCount;
      spectralIndex = s.spectralIndex;
      sequence = s.sequence;
      preciseSequence = s.preciseSequence;
      delete scans;
      scans = new vector<stpMS2>;
      for (i = 0; i<s.scans->size(); i++) scans->push_back(s.scans->at(i));
      delete protIndex;
      protIndex = new vector<size_t>;
      for (i = 0; i<s.protIndex->size(); i++) protIndex->push_back(s.protIndex->at(i));
    }
    return *this;
  }
  bool operator==(const stpPeptide& s){
    if(s.charge!=charge) return false;
    if(s.sequence.compare(sequence)!=0) return false;
    return true;
  }
}stpPeptide;

typedef struct stpdPep{
  size_t index;
  double dSC;
  double dSI;
  double dSIn;
  stpdPep(){
    clear();
  }
  void clear(){
    index=0;
    dSC=0;
    dSI=0;
    dSIn=0;
  }
}stpdPep;

typedef struct stpPepRef{
  vector<stpdPep>* v;
  stpPepRef(){
    v=new vector<stpdPep>;
  }
  stpPepRef(const stpPepRef& s){
    v = new vector<stpdPep>;
    for(size_t i=0;i<s.v->size();i++)v->push_back(s.v->at(i));
  }
  ~stpPepRef(){
    delete v;
  }
  stpPepRef& operator=(const stpPepRef& s){
    if(this!=&s){
      delete v;
      v = new vector<stpdPep>;
      for (size_t i = 0; i<s.v->size(); i++)v->push_back(s.v->at(i));
    }
    return *this;
  }
  stpdPep& operator[](const size_t& index){
    return v->at(index);
  }
  void clear(){
    v->clear();
  }
  void push_back(stpdPep& p){
    v->push_back(p);
  }
  size_t size(){
    return v->size();
  }
}stpPepRef;

typedef struct stpProtein{
  int    length;
  double logNSAF;
  double logSIN;
  double ng;
  double ngC;
  double NSAF;
  double SC;
  double SAF;
  double SI;
  double SIN;
  string name;
  string description;
  stpPepRef peptides;
  stpProtein(){
    length=0;
    logNSAF=0;
    logSIN=0;
    NSAF=0;
    ng=0;
    ngC=0;
    SAF=0;
    SC=0;
    SI=0;
    SIN=0;
    name.clear();
    description.clear();
    peptides.clear();
  }
}stpProtein;

typedef struct stpGroup{
  int groupNumber;
  vector<stpProtein>* proteins;
  vector<stpPeptide>* peptides;
  stpGroup(){
    groupNumber=0;
    proteins=new vector<stpProtein>;
    peptides=new vector<stpPeptide>;
  }
  stpGroup(const stpGroup& s){
    size_t i;
    groupNumber = s.groupNumber;
    proteins = new vector<stpProtein>;
    peptides = new vector<stpPeptide>;
    for(i=0;i<s.proteins->size();i++) proteins->push_back(s.proteins->at(i));
    for(i=0;i<s.peptides->size();i++) peptides->push_back(s.peptides->at(i));
  }
  ~stpGroup(){
    delete proteins;
    delete peptides;
  }
  stpGroup& operator=(const stpGroup& s){
    if(this!=&s){
      size_t i;
      groupNumber = s.groupNumber;
      delete proteins;
      delete peptides;
      proteins = new vector<stpProtein>;
      peptides = new vector<stpPeptide>;
      for (i = 0; i<s.proteins->size(); i++) proteins->push_back(s.proteins->at(i));
      for (i = 0; i<s.peptides->size(); i++) peptides->push_back(s.peptides->at(i));
    }
    return *this;
  }
}stpGroup;

#endif
