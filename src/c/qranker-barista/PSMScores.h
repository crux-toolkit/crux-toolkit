#ifndef PSMSCORES_H_
#define PSMSCORES_H_
#include <vector>
#include <algorithm>
using namespace std;
#include "DataSet.h"


class PSMScoreHolder{
public:
  double score; 
  int psmind;
  double q;
  double PEP;
  int label;
 PSMScoreHolder():score(0.0),psmind(0),q(0.0),PEP(0.0),label(0){;}
  virtual ~PSMScoreHolder() {;}
};

class PSMScores
{
public:
    PSMScores();
    ~PSMScores();
    void clear();
    vector<PSMScoreHolder>::iterator begin() {return scores.begin();}
    vector<PSMScoreHolder>::iterator end() {return scores.end();}    
    static double pi0;
    double factor;

    int calcOverFDR(double fdr);
    void calcMultiOverFDR(vector<double> &fdr, vector<int> &overFDR);
    inline PSMScoreHolder& operator[](int ix){return scores[ix];}    
    void static fillFeaturesSplit(PSMScores& train,PSMScores& test,Dataset &d, double ratio);
    void static fillFeaturesFull(PSMScores& full,Dataset &d);
    inline int size(){return scores.size();}
protected:
    int neg,pos,posNow;
    vector<PSMScoreHolder> scores;
};



#endif /*PSMSCORES_H_*/
