#ifndef PEPSCORES_H_
#define PEPSCORES_H_
#include <vector>
#include <algorithm>
using namespace std;
#include "DataSet.h"


class PepScoreHolder{
public:
  double score; 
  int pepind;
  double q;
  int label;
  double nsaf;
  double PEP;
 PepScoreHolder():score(0.0),pepind(0),q(0.0),label(0),nsaf(0.0),PEP(0.0){;}
  virtual ~PepScoreHolder() {;}
};

class PepScores
{
public:
	PepScores();
	~PepScores();
    vector<PepScoreHolder>::iterator begin() {return scores.begin();}
    vector<PepScoreHolder>::iterator end() {return scores.end();}    
    static double pi0;
    double factor;
    void clear();

    int calcOverFDR(double fdr);
    void calcMultiOverFDR(vector<double> &fdr, vector<int> &overFDR);
    inline PepScoreHolder& operator[](int ix){return scores[ix];}    
    void static fillFeaturesSplit(PepScores& train,PepScores& test,Dataset &d, double ratio);
    void static fillFeaturesFull(PepScores& full,Dataset &d);
    inline int size(){return scores.size();}
protected:
    int neg,pos,posNow;
    vector<PepScoreHolder> scores;
};



#endif /*PEPSCORES_H_*/
