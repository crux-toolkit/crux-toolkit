#include "CPepXMLSpectrum.h"

using namespace std;

CPepXMLSpectrum::CPepXMLSpectrum(){
  searchID=0;
  charge=0;
  fileID=0;
  scanEnd=0;
  scanNumber=0;
  scanStart=0;
  rTimeSec=0;
  precursorNeutMass=0;
  ID.clear();
  nativeID.clear();
  ID.shrink_to_fit();
  nativeID.shrink_to_fit();
  psms = new vector<CPepXMLPSM>;
}

CPepXMLSpectrum::CPepXMLSpectrum(const CPepXMLSpectrum& p){
  searchID=p.searchID;
  charge=p.charge;
  fileID=p.fileID;
  scanEnd=p.scanEnd;
  scanNumber=p.scanNumber;
  scanStart=p.scanStart;
  rTimeSec=p.rTimeSec;
  precursorNeutMass=p.precursorNeutMass;
  ID=p.ID;
  nativeID=p.nativeID;
  ID.shrink_to_fit();
  nativeID.shrink_to_fit();
  psms = new vector<CPepXMLPSM>(*p.psms);
  //for(size_t i=0;i<p.psms->size();i++) psms->push_back(p.psms->at(i));
}

CPepXMLSpectrum::~CPepXMLSpectrum(){
  delete psms;
}

CPepXMLSpectrum& CPepXMLSpectrum::operator=(const CPepXMLSpectrum& p){
  if(this!=&p){
    searchID = p.searchID;
    charge=p.charge;
    fileID=p.fileID;
    scanEnd=p.scanEnd;
    scanNumber=p.scanNumber;
    scanStart=p.scanStart;
    rTimeSec=p.rTimeSec;
    precursorNeutMass=p.precursorNeutMass;
    ID=p.ID;
    nativeID=p.nativeID;
    ID.shrink_to_fit();
    nativeID.shrink_to_fit();
    delete psms;
    psms = new vector<CPepXMLPSM>(*p.psms);
    //for(size_t i=0;i<p.psms->size();i++) psms->push_back(p.psms->at(i));
  }
  return *this;
}

CPepXMLPSM& CPepXMLSpectrum::operator [](const size_t &index){
  return psms->at(index);
}

void CPepXMLSpectrum::clear(){
  searchID=0;
  charge=0;
  fileID=0;
  scanEnd=0;
  scanNumber=0;
  scanStart=0;
  rTimeSec=0;
  precursorNeutMass=0;
  ID.clear();
  nativeID.clear();
  ID.shrink_to_fit();
  nativeID.shrink_to_fit();
  delete psms;
  psms = new vector<CPepXMLPSM>;
}

size_t CPepXMLSpectrum::size(){
  return psms->size();
}

size_t CPepXMLSpectrum::sizeOf(){
  size_t i;
  i=sizeof(*this);
  i+=ID.capacity();
  i+=nativeID.capacity();
  for(size_t j=0;j<psms->size();j++) i+=psms->at(j).sizeOf();
  return i;
}
