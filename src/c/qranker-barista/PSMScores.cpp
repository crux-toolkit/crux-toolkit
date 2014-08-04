#include <assert.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <set>
#include <vector>
#include <string>
#include <math.h>
using namespace std;
#include "PSMScores.h"
#include "utils.h"

inline bool operator>(const PSMScoreHolder &one, const PSMScoreHolder &other) 
    {return (one.score>other.score);}

inline bool operator<(const PSMScoreHolder &one, const PSMScoreHolder &other) 
    {return (one.score<other.score);}



PSMScores::PSMScores()
{
    factor=1;
    neg=0;
    pos=0;
    posNow=0;
}

PSMScores::~PSMScores()
{

}

void PSMScores::clear()
{
  scores.clear();
  pos=neg=posNow = 0;
}


/**
 * Percentage of target scores that are drawn according to the null.
 */
double PSMScores::pi0 = 0.9;


/**
 * Calculate the number of targets that score above a specified FDR.
 */
int PSMScores::calcOverFDR(double fdr) {
  
  vector<PSMScoreHolder>::iterator it = scores.begin();

  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
  
  int positives=0,nulls=0;
  double efp=0.0,q;
  posNow = 0;
  register unsigned int ix=0;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      positives++;
    if (it->label==-1) {
      nulls++;
      efp=pi0*nulls*factor;
    }
    if (positives)
      q=efp/(double)positives;
    else
      q=pi0;
    if (q>pi0)
      q=pi0;
    it->q=q;
    if (fdr>=q)
      posNow = positives;
  }
  for (ix=scores.size();--ix;) {
    if (scores[ix-1].q > scores[ix].q)
      scores[ix-1].q = scores[ix].q;  
  }
  return posNow;
}


void PSMScores::calcMultiOverFDR(vector<double> &fdr, vector<int> &overFDR) {
  
  vector<PSMScoreHolder>::iterator it = scores.begin();

  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
  
  int positives=0,nulls=0;
  double efp=0.0,q;
  register unsigned int ix=0;
  
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label==1)
      positives++;
    if (it->label==-1) {
      nulls++;
      efp=pi0*nulls*factor;
    }
    if (positives)
      q=efp/(double)positives;
    else
      q=pi0;
    it->q=q;
    if (q>pi0)
      q=pi0;
       
    for(unsigned int ct = 0; ct < fdr.size(); ct++)
      if (fdr[ct]>=q)
	overFDR[ct] = positives;
  }
  for (ix=scores.size();--ix;) {
    if (scores[ix-1].q > scores[ix].q)
      scores[ix-1].q = scores[ix].q;  
  }
}




void PSMScores::fillFeaturesSplit(PSMScores& train,PSMScores& test, Dataset& d, double ratio) {
 
  int num_psms = d.get_num_psms();
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  int n = num_psms;
  int k = (int)(n*ratio);
  int l = n-k;
 
  PSMScoreHolder s;
  
  //collect everybody 
  vector<PSMScoreHolder> all_examples;
  all_examples.resize(n,s);
  for (int i = 0; i < num_psms; i++)
    {
      all_examples[i].psmind = i;
      all_examples[i].label = d.psmind2label(i);
    }
  
  //mix up the examples
  for(int i = 0; i < n; i++)
    {
      int p1 = (int)((double)myrandom()/UNIFORM_INT_DISTRIBUTION_MAX*(n-1)); 
      int p2 = (int)((double)myrandom()/UNIFORM_INT_DISTRIBUTION_MAX*(n-1)); 
      s = all_examples[p1];
      all_examples[p1] = all_examples[p2];
      all_examples[p2] = s;
    }
  
  train.scores.resize(k,s);
  test.scores.resize(l,s);
  

  int num_pos_train = 0;
  int num_neg_train = 0;
  int num_pos_test = 0;
  int num_neg_test = 0;
  //distribute the normal set between train and test
  for (int i = 0; i < k; i++)
    {
      train.scores[i] = all_examples[i];
      if (train.scores[i].label == 1)
	num_pos_train++;
      else
	num_neg_train++;
    }
  for(int i=0; i < l; i++)
    {
      test.scores[i] = all_examples[i+k];
      if (test.scores[i].label == 1)
	num_pos_test++;
      else
	num_neg_test++;
    }
  //cout << num_pos_train << " " << num_neg_train << " " << num_pos_test << " " << num_neg_test << "\n";
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=num_neg_train;
  test.neg=num_neg_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}


void PSMScores::fillFeaturesFull(PSMScores& full, Dataset& d) {
 
  int n = d.get_num_psms();
 
  PSMScoreHolder s;
  full.scores.resize(n,s);
  int num_pos = 0;
  int num_neg = 0;

  for (int i = 0; i < n; i++)
    {
      full.scores[i].psmind = i;
      full.scores[i].label = d.psmind2label(i);
      if (full.scores[i].label == 1)
	num_pos++;
      else
	num_neg++;
  }
  
  full.pos=num_pos;
  full.neg=num_neg;
  full.factor = full.pos/(double)full.neg;
}









