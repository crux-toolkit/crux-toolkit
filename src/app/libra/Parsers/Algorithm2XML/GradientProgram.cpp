#include "Parsers/Algorithm2XML/GradientProgram.h"

GradientProgram::GradientProgram() { 
  secs_= new Array<double>(0);
  acns_= new Array<double>(0);
}

GradientProgram:: GradientProgram(const char* file) { 
  double secs, acn;
  ifstream fin(file);
  secs_= new Array<double>();
  acns_= new Array<double>();
  while (1) {
    fin >> secs >> acn;

    if (fin.fail()) break;
    
    secs_->insertAtEnd(secs);
    acns_->insertAtEnd(acn);

  }
}


double GradientProgram::getAcn(double rtsec) {
  int i; 

  if (secs_->size() < 1) return rtsec;

  for (i=0; (*secs_)[i] < rtsec && i < secs_->size() ; i++);
  

  if (rtsec >= (*secs_)[ secs_->size()-1]) {
    //cerr << "WARNING: RT " << rtsec << " outside of GradientProgram range." << endl;
    return (*acns_)[ acns_->size()-1];
  }
  
  if (i==0) {
    return (*acns_)[i];
  }

  double leftY = (*acns_)[i-1];
  double rightY  = (*acns_)[i];

  double leftX = (*secs_)[i-1];
  double rightX  = (*secs_)[i];



  double A = (rightY - leftY) / (rightX - leftX);

  double B = rightY - A * rightX;

  return A * rtsec + B;
  
   

}

GradientProgram::~GradientProgram() { 
  delete secs_;
  delete acns_;
}
