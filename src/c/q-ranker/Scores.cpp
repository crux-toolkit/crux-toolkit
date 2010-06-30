/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Scores.cpp,v 1.69.2.5 2009/08/01 20:29:50 jasonw Exp $
 *******************************************************************************/
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
#include "DataSet.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Scores.h"
#include "Globals.h"
#include "PosteriorEstimator.h"
#include "ssl.h"

namespace qranker {

inline bool operator>(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score>other.score);}

inline bool operator<(const ScoreHolder &one, const ScoreHolder &other) 
    {return (one.score<other.score);}

Scores::Scores()
{
    factor=1;
    neg=0;
    pos=0;
    posNow=0;
}

Scores::~Scores()
{
}

/**
 * Percentage of target scores that are drawn according to the null.
 */
double Scores::pi0 = 0.9;

void Scores::merge(vector<Scores>& sv) {
  scores.clear();
  for(vector<Scores>::iterator a = sv.begin();a!=sv.end();a++) {
  	a->normalizeScores();
    copy(a->begin(),a->end(),back_inserter(scores));
  }
  sort(scores.begin(),scores.end(), greater<ScoreHolder>() );
}


void Scores::printRoc(string & fn){
 ofstream rocStream(fn.data(),ios::out);
 vector<ScoreHolder>::iterator it;
 for(it=scores.begin();it!=scores.end();it++) {
   rocStream << (it->label==-1?-1:1) << endl;
 }
 rocStream.close();
}	

double Scores::calcScore(const double * feat) const{
  register int ix=FeatureNames::getNumFeatures();
  register double score = w_vec[ix];
  for(;ix--;) {
  	score += feat[ix]*w_vec[ix];
  }
  return score;
}

ScoreHolder* Scores::getScoreHolder(const double *d){
  if (scoreMap.size()==0) {
    vector<ScoreHolder>::iterator it;
    for(it=scores.begin();it!=scores.end();it++) {
      scoreMap[it->pPSM->features] = &(*it);
    }
  }
  return scoreMap[d];  
}


void Scores::fillFeatures(SetHandler& norm,SetHandler& shuff) {
  scores.clear();
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
  while((pPSM=normIter.getNext())!=NULL) {
    scores.push_back(ScoreHolder(.0,1,pPSM));
  }
  while((pPSM=shuffIter.getNext())!=NULL) {
    scores.push_back(ScoreHolder(.0,-1,pPSM));
  }
  pos = norm.getSize(); neg = shuff.getSize();
  factor=norm.getSize()/(double)shuff.getSize();
}


void Scores::createXvalSets(vector<Scores>& train,vector<Scores>& test, const unsigned int xval_fold) {
  train.resize(xval_fold);
  test.resize(xval_fold);
  vector<size_t> remain(xval_fold);
  size_t fold = xval_fold, ix = scores.size();
  while (fold--) {
    remain[fold] = ix / (fold + 1);
    ix -= remain[fold];
  }
  
  for(unsigned int j=0;j<scores.size();j++) {
    ix = rand()%(scores.size()-j);
    fold = 0;
    while(ix>remain[fold])
      ix-= remain[fold++];
    for(unsigned int i=0;i<xval_fold;i++) {
      if(i==fold) {
        test[i].scores.push_back(scores[j]);
      } else {
        train[i].scores.push_back(scores[j]);
      }
    }
    --remain[fold];  
  }
  vector<ScoreHolder>::const_iterator it;
  for(unsigned int i=0;i<xval_fold;i++) {
  	train[i].pos=0;train[i].neg=0;
  	for(it=train[i].begin();it!=train[i].end();it++) {
      if (it->label==1) train[i].pos++;
      else train[i].neg++;
    }
    train[i].factor=train[i].pos/(double)train[i].neg;
    test[i].pos=0;test[i].neg=0;
  	for(it=test[i].begin();it!=test[i].end();it++) {
      if (it->label==1) test[i].pos++;
      else test[i].neg++;
  	}
    test[i].factor=test[i].pos/(double)test[i].neg;
  }
}

void Scores::normalizeScores() {
  // Normalize scores so that distance between 1st and 3rd quantile of the null scores are 1
  unsigned int q1index = neg/4,q3index = neg*3/4,decoys=0;
  vector<ScoreHolder>::iterator it = scores.begin();
  double q1 = it->score;
  double q3 = q1 + 1.0;
  while(it!=scores.end()) {
  	if (it->label == -1) {
  	  if(++decoys==q1index)
  	    q1=it->score;
  	  else if (decoys==q3index) {
  	    q3=it->score;
  	    break;
  	  }
  	}
    ++it;
  }
  double diff = q1-q3;
  assert(diff>0);
  for(it=scores.begin();it!=scores.end();++it) { 
    it->score -= q1;
    it->score /= diff;
  }
}

int Scores::calcScores(vector<double>& w,double fdr) {
  w_vec=w;
  const double * features;
  unsigned int ix;
  vector<ScoreHolder>::iterator it = scores.begin();
  while(it!=scores.end()) {
    features = it->pPSM->features;
  	it->score = calcScore(features);
    it++;
  }
  sort(scores.begin(),scores.end(),greater<ScoreHolder>());
  if (VERB>3) {
    cerr << "10 best scores and labels" << endl;
    for (ix=0;ix < 10;ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
    cerr << "10 worst scores and labels" << endl;
    for (ix=scores.size()-10;ix < scores.size();ix++) {
  	  cerr << scores[ix].score << " " << scores[ix].label << endl;
    }
  }
  int positives=0,nulls=0;
  double efp=0.0,q;
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
    it->pPSM->q=q;
    if (fdr>=q)
      posNow = positives;
  }
  for (ix=scores.size();--ix;) {
    if (scores[ix-1].pPSM->q > scores[ix].pPSM->q)
      scores[ix-1].pPSM->q = scores[ix].pPSM->q;  
  }
  return posNow;
}

void Scores::generateNegativeTrainingSet(AlgIn& data,const double cneg) {
  unsigned int ix1=0,ix2=0;
  for(ix1=0;ix1<size();ix1++) {
    if (scores[ix1].label==-1) {
      data.vals[ix2]=scores[ix1].pPSM->features;
      data.Y[ix2]=-1;
      data.C[ix2++]=cneg;
    }
  }
  data.negatives=ix2;
}


void Scores::generatePositiveTrainingSet(AlgIn& data,const double fdr,const double cpos) {
  unsigned int ix1=0,ix2=data.negatives,p=0;
  for(ix1=0;ix1<size();ix1++) {
    if (scores[ix1].label==1) {
      if (fdr<scores[ix1].pPSM->q) {
        posNow=p;
        break;
      }
      data.vals[ix2]=scores[ix1].pPSM->features;
      data.Y[ix2]=1;
      data.C[ix2++]=cpos;
      ++p;
    }
  }
  data.positives=p;
  data.m=ix2;
}

void Scores::recalculateDescriptionOfGood(const double fdr) {
  doc.clear();
  unsigned int ix1=0;
  for(ix1=0;ix1<size();ix1++) {
    if (scores[ix1].label==1) {
      if (fdr>scores[ix1].pPSM->q) {
        doc.registerCorrect(scores[ix1].pPSM);
      }
    }
  }
  doc.trainCorrect();
  for(ix1=0;ix1<size();ix1++) {
    doc.setFeatures(scores[ix1].pPSM);
  }
}

int Scores::getInitDirection(const double fdr, vector<double>& direction, bool findDirection) {
  int bestPositives = -1;
  int bestFeature =-1;
  bool lowBest = false;
  
  if (findDirection) { 
    for (unsigned int featNo=0;featNo<FeatureNames::getNumFeatures();featNo++) {
      vector<ScoreHolder>::iterator it = scores.begin();
      while(it!=scores.end()) {
        it->score = it->pPSM->features[featNo];
        it++;
      }
      sort(scores.begin(),scores.end());
      for (int i=0;i<2;i++) {
        int positives=0,nulls=0;
        double efp=0.0,q;
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
          if (fdr<=q) {
            if (positives>bestPositives && scores.begin()->score!=it->score) {
              bestPositives=positives;
              bestFeature = featNo;
              lowBest = (i==0);
            }
            if (i==0) {
              reverse(scores.begin(),scores.end());
            }
            break;
          }
        }
      }
    }
    for (int ix=FeatureNames::getNumFeatures();ix--;) {
      direction[ix]=0;
    }
    direction[bestFeature]=(lowBest?-1:1);
    if (VERB>1) {
      cerr << "Selected feature number " << bestFeature +1 << " as initial search direction, could separate " << 
              bestPositives << " positives in that direction" << endl;
    }
  } else {
    bestPositives = calcScores(direction,fdr);
    if (VERB>1) {
      cerr << "Found " << 
              bestPositives << " positives in the initial search direction" << endl;
    }    
  }
  return bestPositives;
}

double Scores::estimatePi0() {
  vector<pair<double,bool> > combined;
  transform(scores.begin(),scores.end(),back_inserter(combined),  mem_fun_ref(&ScoreHolder::toPair));

  // Estimate pi0
  pi0 = PosteriorEstimator::estimatePi0(combined);
  return pi0;

}

void Scores::calcPep() {

  vector<pair<double,bool> > combined;
  transform(scores.begin(),scores.end(),back_inserter(combined),  mem_fun_ref(&ScoreHolder::toPair));
  vector<double> peps;                                                                                                                  
  
  // Logistic regression on the data
  PosteriorEstimator::estimatePEP(combined,pi0,peps);

  size_t pix=0;
  for(size_t ix=0; ix<scores.size(); ix++) {
    (scores[ix]).pPSM->pep = peps[pix];
    if(scores[ix].label==1)
      ++pix;
  }
}



/**
 * Calculate the number of targets that score above a specified FDR.
 */
int Scores::calcOverFDR(double fdr) {
  
  vector<ScoreHolder>::iterator it = scores.begin();

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
    it->pPSM->q=q;
    if (fdr>=q)
      posNow = positives;
  }
  //for (ix=scores.size();--ix;) {
  //if (scores[ix-1].pPSM->q > scores[ix].pPSM->q)
  //  scores[ix-1].pPSM->q = scores[ix].pPSM->q;  
  //}
  return posNow;
}


void Scores::calcMultiOverFDR(vector<double> &fdr, vector<int> &overFDR) {
  
  vector<ScoreHolder>::iterator it = scores.begin();

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
    it->pPSM->q=q;
    if (q>pi0)
      q=pi0;
       
    for(unsigned int ct = 0; ct < fdr.size(); ct++)
      if (fdr[ct]>=q)
	overFDR[ct] = positives;
  }
  for (ix=scores.size();--ix;) {
    if (scores[ix-1].pPSM->q > scores[ix].pPSM->q)
      scores[ix-1].pPSM->q = scores[ix].pPSM->q;  
  }
}




void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
 
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  //shuffled
  int n = shuff.getSize();
  int k = (int)(shuff.getSize()*ratio);
  int l = shuff.getSize() - k;
  //normal
  int nz  = norm.getSize();
  
 
  ScoreHolder s;
  
  int i=0;
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);
  
  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }
      

  //1 for trainset, 2 for testset
  vector<int> pos_assignments(nz,0);
  int num_pos_train = 0;
  int num_pos_test = 0;
  for (int count = 0; count < nz; count++)
    {
      if (num_pos_train < num_pos_test)
	{
	  //want to add to the trainset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      pos_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_pos_train++;
	    }
	  else
	    {
	      pos_assignments[count] = 2;
	      num_pos_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      pos_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_pos_test++;
	    }
	  else
	    {
	      pos_assignments[count] = 1;
	      num_pos_train++;
	    }
	}
    }
  
  
  //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n,s);
  while((pPSM=shuffIter.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  assert(i == n);
 
  //mix up the examples
   for(i = 0; i < n; i++)
     {
       int p1 = (int)((double)rand()/RAND_MAX*(n-1)); 
       int p2 = (int)((double)rand()/RAND_MAX*(n-1)); 
       s = all_neg_examples[p1];
       all_neg_examples[p1] = all_neg_examples[p2];
       all_neg_examples[p2] = s;
     }
  
  train.scores.resize(num_pos_train+k,s);
  test.scores.resize(num_pos_test+l,s);
  
    
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  for (int count = 0; count < nz; count++)
    {
      if (pos_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_pos_examples[count];
	  ++ix1;
      }
      else if (pos_assignments[count] == 2)
	{
	  test.scores[ix2] = all_pos_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid pos_assignments_value " << pos_assignments[count] << "\n";
    }
  

  //distribute the shuffled set between train and test
  int count = 0;
  while(count < k)
    {
      train.scores[ix1] = all_neg_examples[count];
      ++ix1;
      ++count;
    }
  while (count < n)
    {
      test.scores[ix2] = all_neg_examples[count];
      ++ix2;    
      count++;
    }
  
        
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=ix1-num_pos_train;
  test.neg=ix2-num_pos_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}



void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff,SetHandler& shuff1, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff),shuffIter1(&shuff1), normIter(&norm);
 
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  //shuffled
  int n = shuff.getSize();
  int k = (int)(shuff.getSize()*ratio);
  int l = shuff.getSize() - k;
  //shuffled1
  int nt = shuff1.getSize();
  int kkk = (int)(shuff1.getSize()*ratio);
  int lll = shuff1.getSize() - kkk;
  //normal
  int nz  = norm.getSize();
  
  ScoreHolder s;
  
  int i=0;
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);
  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }
  //1 for trainset, 2 for testset
  vector<int> pos_assignments(nz,0);
  int num_pos_train = 0;
  int num_pos_test = 0;
  for (int count = 0; count < nz; count++)
    {
      if (num_pos_train < num_pos_test)
	{
	  //want to add to the trainset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      pos_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_pos_train++;
	    }
	  else
	    {
	      pos_assignments[count] = 2;
	      num_pos_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      pos_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_pos_test++;
	    }
	  else
	    {
	      pos_assignments[count] = 1;
	      num_pos_train++;
	    }
	}
    }
  
  
  //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n+nt,s);
  while((pPSM=shuffIter.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  while((pPSM=shuffIter1.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
   assert(i == n+nt);
  //mix up the examples
   for(i = 0; i < n+nt; i++)
     {
       int p1 = (int)((double)rand()/RAND_MAX*(n+nt-1)); 
       int p2 = (int)((double)rand()/RAND_MAX*(n+nt-1)); 
       s = all_neg_examples[p1];
       all_neg_examples[p1] = all_neg_examples[p2];
       all_neg_examples[p2] = s;
     }
     

  
 
  train.scores.resize(num_pos_train+k+kkk,s);
  test.scores.resize(num_pos_test+l+lll,s);
  
  
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  for (int count = 0; count < nz; count++)
    {
      if (pos_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_pos_examples[count];
	  ++ix1;
      }
      else if (pos_assignments[count] == 2)
	{
	  test.scores[ix2] = all_pos_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid pos_assignments_value " << pos_assignments[count] << "\n";
    }
  

  
  //distribute the shuffled set between train and test
  int count = 0;
  while(count < k+kkk)
    {
      train.scores[ix1] = all_neg_examples[count];
      ++ix1;
      ++count;
    }
  while (count < n+nt)
    {
      test.scores[ix2] = all_neg_examples[count];
      ++ix2;    
      count++;
    }
  

        
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=ix1-num_pos_train;
  test.neg=ix2-num_pos_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}




void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff,SetHandler& shuff1,SetHandler& shuff2, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff),shuffIter1(&shuff1),shuffIter2(&shuff2), normIter(&norm);
 
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  //shuffled
  int n = shuff.getSize();
  int k = (int)(shuff.getSize()*ratio);
  int l = shuff.getSize() - k;
  //shuffled1
  int nt = shuff1.getSize();
  int kkk = (int)(shuff1.getSize()*ratio);
  int lll = shuff1.getSize() - kkk;
  //shuffled2
  int nt1 = shuff2.getSize();
  int kkk1 = (int)(shuff2.getSize()*ratio);
  int lll1 = shuff2.getSize() - kkk1;
  //normal
  int nz  = norm.getSize();
  
  ScoreHolder s;
  
  int i=0;
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);
  
  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }
      

  //1 for trainset, 2 for testset
  vector<int> pos_assignments(nz,0);
  int num_pos_train = 0;
  int num_pos_test = 0;
  for (int count = 0; count < nz; count++)
    {
      if (num_pos_train < num_pos_test)
	{
	  //want to add to the trainset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      pos_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_pos_train++;
	    }
	  else
	    {
	      pos_assignments[count] = 2;
	      num_pos_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      pos_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_pos_test++;
	    }
	  else
	    {
	      pos_assignments[count] = 1;
	      num_pos_train++;
	    }
	}
    }
  
   //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n+nt+nt1,s);
  while((pPSM=shuffIter.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  while((pPSM=shuffIter1.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  while((pPSM=shuffIter2.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
   assert(i == n+nt+nt1);
 
  //mix up the examples
   for(i = 0; i < n+nt+nt1; i++)
     {
       int p1 = (int)((double)rand()/RAND_MAX*(n+nt+nt1-1)); 
       int p2 = (int)((double)rand()/RAND_MAX*(n+nt+nt1-1)); 
       s = all_neg_examples[p1];
       all_neg_examples[p1] = all_neg_examples[p2];
       all_neg_examples[p2] = s;
     }
     

  
  train.scores.resize(num_pos_train+k+kkk+kkk1,s);
  test.scores.resize(num_pos_test+l+lll+lll1,s);
  
  
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  for (int count = 0; count < nz; count++)
    {
      if (pos_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_pos_examples[count];
	  ++ix1;
      }
      else if (pos_assignments[count] == 2)
	{
	  test.scores[ix2] = all_pos_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid pos_assignments_value " << pos_assignments[count] << "\n";
    }
  

  //distribute the shuffled set between train and test
  int count = 0;
  while(count < k+kkk+kkk1)
    {
      train.scores[ix1] = all_neg_examples[count];
      ++ix1;
      ++count;
    }
  while (count < n+nt+nt1)
    {
      test.scores[ix2] = all_neg_examples[count];
      ++ix2;    
      count++;
    }
  
  

      
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=ix1-num_pos_train;
  test.neg=ix2-num_pos_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}







} // qranker namspace


/*********************************************************
 *
 * Old split functions taking into account peptide sequences
 *
 */


/*
void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
 
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  //shuffled
  int n = shuff.getSize();
  //normal
  int nz  = norm.getSize();
  
 
  ScoreHolder s;
  
  int i=0;
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);
  
  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }
      

  //1 for trainset, 2 for testset
  vector<int> pos_assignments(nz,0);
  int num_pos_train = 0;
  int num_pos_test = 0;
  for (int count = 0; count < nz; count++)
    {
      if (num_pos_train < num_pos_test)
	{
	  //want to add to the trainset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      pos_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_pos_train++;
	    }
	  else
	    {
	      pos_assignments[count] = 2;
	      num_pos_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      pos_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_pos_test++;
	    }
	  else
	    {
	      pos_assignments[count] = 1;
	      num_pos_train++;
	    }
	}
    }
  
  
  //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n,s);
  while((pPSM=shuffIter.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  assert(i == n);
 
  //mix up the examples
   for(i = 0; i < n; i++)
     {
       int p1 = (int)((double)rand()/RAND_MAX*(n-1)); 
       int p2 = (int)((double)rand()/RAND_MAX*(n-1)); 
       s = all_neg_examples[p1];
       all_neg_examples[p1] = all_neg_examples[p2];
       all_neg_examples[p2] = s;
     }
     
  //1 for trainset, 2 for testset
  vector<int> neg_assignments(n,0);
  
  int num_neg_train = 0;
  int num_neg_test = 0;
  for (int count = 0; count < n; count++)
    {
      if (num_neg_train < num_neg_test)
	{
	  //want to add to the trainset
	  string pep = all_neg_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      neg_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_neg_train++;
	    }
	  else
	    {
	      neg_assignments[count] = 2;
	      num_neg_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_neg_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      neg_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_neg_test++;
	    }
	  else
	    {
	      neg_assignments[count] = 1;
	      num_neg_train++;
	    }
	}
    }
  
  train.scores.resize(num_pos_train+num_neg_train,s);
  test.scores.resize(num_pos_test+num_neg_test,s);
  
    
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  for (int count = 0; count < nz; count++)
    {
      if (pos_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_pos_examples[count];
	  ++ix1;
      }
      else if (pos_assignments[count] == 2)
	{
	  test.scores[ix2] = all_pos_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid pos_assignments_value " << pos_assignments[count] << "\n";
    }
  

  for (int count = 0; count < n; count++)
    {
      if (neg_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_neg_examples[count];
	  ++ix1;
      }
      else if (neg_assignments[count] == 2)
	{
	  test.scores[ix2] = all_neg_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid neg_assignments_value " << neg_assignments[count] << "\n";
    }
        
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=num_neg_train;
  test.neg=num_neg_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}


void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff,SetHandler& shuff1, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff),shuffIter1(&shuff1), normIter(&norm);
 
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  //shuffled
  int n = shuff.getSize();
  //shuffled1
  int nt = shuff1.getSize();
  //normal
  int nz  = norm.getSize();
  
  ScoreHolder s;
  
  int i=0;
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);
  
  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }
      

  //1 for trainset, 2 for testset
  vector<int> pos_assignments(nz,0);
  int num_pos_train = 0;
  int num_pos_test = 0;
  for (int count = 0; count < nz; count++)
    {
      if (num_pos_train < num_pos_test)
	{
	  //want to add to the trainset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      pos_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_pos_train++;
	    }
	  else
	    {
	      pos_assignments[count] = 2;
	      num_pos_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      pos_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_pos_test++;
	    }
	  else
	    {
	      pos_assignments[count] = 1;
	      num_pos_train++;
	    }
	}
    }
  
  
  //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n+nt,s);
  while((pPSM=shuffIter.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  while((pPSM=shuffIter1.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
   assert(i == n+nt);
 
  //mix up the examples
   for(i = 0; i < n+nt; i++)
     {
       int p1 = (int)((double)rand()/RAND_MAX*(n+nt-1)); 
       int p2 = (int)((double)rand()/RAND_MAX*(n+nt-1)); 
       s = all_neg_examples[p1];
       all_neg_examples[p1] = all_neg_examples[p2];
       all_neg_examples[p2] = s;
     }
     
  //1 for trainset, 2 for testset
  vector<int> neg_assignments(n+nt,0);
  
  int num_neg_train = 0;
  int num_neg_test = 0;
  for (int count = 0; count < n+nt; count++)
    {
      if (num_neg_train < num_neg_test)
	{
	  //want to add to the trainset
	  string pep = all_neg_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      neg_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_neg_train++;
	    }
	  else
	    {
	      neg_assignments[count] = 2;
	      num_neg_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_neg_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      neg_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_neg_test++;
	    }
	  else
	    {
	      neg_assignments[count] = 1;
	      num_neg_train++;
	    }
	}
    }
  
 
  train.scores.resize(num_pos_train+num_neg_train,s);
  test.scores.resize(num_pos_test+num_neg_test,s);
  
  
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  for (int count = 0; count < nz; count++)
    {
      if (pos_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_pos_examples[count];
	  ++ix1;
      }
      else if (pos_assignments[count] == 2)
	{
	  test.scores[ix2] = all_pos_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid pos_assignments_value " << pos_assignments[count] << "\n";
    }
  

  for (int count = 0; count < n+nt; count++)
    {
      if (neg_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_neg_examples[count];
	  ++ix1;
      }
      else if (neg_assignments[count] == 2)
	{
	  test.scores[ix2] = all_neg_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid neg_assignments_value " << neg_assignments[count] << "\n";
    }
        
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=num_neg_train;
  test.neg=num_neg_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}




void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff,SetHandler& shuff1,SetHandler& shuff2, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff),shuffIter1(&shuff1),shuffIter2(&shuff2), normIter(&norm);
 
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  //shuffled
  int n = shuff.getSize();
  //shuffled1
  int nt = shuff1.getSize();
  //shuffled2
  int nt1 = shuff2.getSize();
  //normal
  int nz  = norm.getSize();
  
  ScoreHolder s;
  
  int i=0;
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);
  
  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }
      

  //1 for trainset, 2 for testset
  vector<int> pos_assignments(nz,0);
  int num_pos_train = 0;
  int num_pos_test = 0;
  for (int count = 0; count < nz; count++)
    {
      if (num_pos_train < num_pos_test)
	{
	  //want to add to the trainset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      pos_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_pos_train++;
	    }
	  else
	    {
	      pos_assignments[count] = 2;
	      num_pos_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_pos_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      pos_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_pos_test++;
	    }
	  else
	    {
	      pos_assignments[count] = 1;
	      num_pos_train++;
	    }
	}
    }
  
   //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n+nt+nt1,s);
  while((pPSM=shuffIter.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  while((pPSM=shuffIter1.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
  while((pPSM=shuffIter2.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
   }
   assert(i == n+nt+nt1);
 
  //mix up the examples
   for(i = 0; i < n+nt+nt1; i++)
     {
       int p1 = (int)((double)rand()/RAND_MAX*(n+nt+nt1-1)); 
       int p2 = (int)((double)rand()/RAND_MAX*(n+nt+nt1-1)); 
       s = all_neg_examples[p1];
       all_neg_examples[p1] = all_neg_examples[p2];
       all_neg_examples[p2] = s;
     }
     
  //1 for trainset, 2 for testset
  vector<int> neg_assignments(n+nt+nt1,0);
  int num_neg_train = 0;
  int num_neg_test = 0;
  for (int count = 0; count < n+nt+nt1; count++)
    {
      if (num_neg_train < num_neg_test)
	{
	  //want to add to the trainset
	  string pep = all_neg_examples[count].pPSM->peptide;
	  //but make sure the pep is not already in the testset
	  if (peptides_test.count(pep) == 0)
	    {
	      neg_assignments[count] = 1;
	      peptides_train.insert(pep);
	      num_neg_train++;
	    }
	  else
	    {
	      neg_assignments[count] = 2;
	      num_neg_test++;
	    }
	}
      else
	{
	  //want to add to the testset
	  string pep = all_neg_examples[count].pPSM->peptide;
	  //but make sure it's not already in the trainset
	  if(peptides_train.count(pep) == 0)
	    {
	      neg_assignments[count] = 2;
	      peptides_test.insert(pep);
	      num_neg_test++;
	    }
	  else
	    {
	      neg_assignments[count] = 1;
	      num_neg_train++;
	    }
	}
    }
  
  train.scores.resize(num_pos_train+num_neg_train,s);
  test.scores.resize(num_pos_test+num_neg_test,s);
  
  
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  for (int count = 0; count < nz; count++)
    {
      if (pos_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_pos_examples[count];
	  ++ix1;
      }
      else if (pos_assignments[count] == 2)
	{
	  test.scores[ix2] = all_pos_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid pos_assignments_value " << pos_assignments[count] << "\n";
    }
  

  for (int count = 0; count < n+nt+nt1; count++)
    {
      if (neg_assignments[count] == 1 )
	{
	  train.scores[ix1] = all_neg_examples[count];
	  ++ix1;
      }
      else if (neg_assignments[count] == 2)
	{
	  test.scores[ix2] = all_neg_examples[count];
	  ++ix2;
	}
      else
	cerr << "invalid neg_assignments_value " << neg_assignments[count] << "\n";
    }
        
  train.pos=num_pos_train;
  test.pos=num_pos_test;
  train.neg=num_neg_train;
  test.neg=num_neg_test;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;
  
}

*/


/*************************************************************
 *
 *
 *  Old split functions NOT taking into account the peptide sequences
 *
 */

/*



void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff), normIter(&norm);
 
  assert(ratio>0 && ratio < 1);
  //cerr << "fill features small\n";
  //shuffled
  int n = shuff.getSize();
  int k = (int)(shuff.getSize()*ratio);
  int l = shuff.getSize() - k;
  //normal
  int nz  = norm.getSize();
  int kk = (int)(norm.getSize()*ratio);
  int ll = norm.getSize() - kk;

  ScoreHolder s;
  train.scores.resize(kk+k,s);
  test.scores.resize(ll+l,s);
   
  int i=0;
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);

  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }

  
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  int count = 0;
  while (count < kk)
    {
      train.scores[ix1] = all_pos_examples[count];
      ++ix1;
      ++count;
    }
  while (count < nz)
    {
      test.scores[ix2] = all_pos_examples[count];
      ++ix2;    
      count++;
    }
  assert(ix1==kk);
  assert(ix2==ll);
  
  //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n,s);
   while((pPSM=shuffIter.getNext())!=NULL) {
     all_neg_examples[i].label=-1;
     all_neg_examples[i].pPSM=pPSM;
      ++i;
   }
   assert(i == n);
 
  //mix up the examples
  for(i = 0; i < n; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(n-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(n-1)); 
      s = all_neg_examples[p1];
      all_neg_examples[p1] = all_neg_examples[p2];
      all_neg_examples[p2] = s;
    }

  
  //distribute the shuffled set between train and test
  count = 0;
  while(count < k){
    train.scores[ix1] = all_neg_examples[count];
      ++ix1;
      ++count;
  }
  while (count < n)
    {
      test.scores[ix2] = all_neg_examples[count];
      ++ix2;    
      count++;
    }
    
  train.pos=kk;
  test.pos=ll;
  train.neg=ix1-kk;
  test.neg=ix2-ll;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;

}



void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff,SetHandler& shuff1, const double ratio) {
  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff),shuffIter1(&shuff1), normIter(&norm);
 
  assert(ratio>0 && ratio < 1);
  cerr << "fill features medium\n";
  //shuffled
  int n = shuff.getSize();
  int k = (int)(shuff.getSize()*ratio);
  int l = shuff.getSize() - k;
  //shuffled1
  int nt = shuff1.getSize();
  int kkk = (int)(shuff1.getSize()*ratio);
  int lll = shuff1.getSize() - kkk;
  //normal
  int nz  = norm.getSize();
  int kk = (int)(norm.getSize()*ratio);
  int ll = norm.getSize() - kk;

  ScoreHolder s;
  train.scores.resize(kk+k+kkk,s);
  test.scores.resize(ll+l+lll,s);
   
  
  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);

  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }

  
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  int count = 0;
  while (count < kk)
    {
      train.scores[ix1] = all_pos_examples[count];
      ++ix1;
      ++count;
    }
  while (count < nz)
    {
      test.scores[ix2] = all_pos_examples[count];
      ++ix2;    
      count++;
    }
  assert(ix1==kk);
  assert(ix2==ll);

  //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n,s);
   while((pPSM=shuffIter.getNext())!=NULL) {
     all_neg_examples[i].label=-1;
     all_neg_examples[i].pPSM=pPSM;
      ++i;
   }
   assert(i == n);
 
  //mix up the examples
  for(i = 0; i < n; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(n-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(n-1)); 
      s = all_neg_examples[p1];
      all_neg_examples[p1] = all_neg_examples[p2];
      all_neg_examples[p2] = s;
    }

  
  //distribute the shuffled set between train and test
  count = 0;
  while(count < k){
    train.scores[ix1] = all_neg_examples[count];
      ++ix1;
      ++count;
  }
  while (count < n)
    {
      test.scores[ix2] = all_neg_examples[count];
      ++ix2;    
      count++;
    }

  //collect everybody negative again
  i = 0;
  vector<ScoreHolder> all_neg_examples1;
  all_neg_examples1.resize(nt,s);
   while((pPSM=shuffIter1.getNext())!=NULL) {
     all_neg_examples1[i].label=-1;
     all_neg_examples1[i].pPSM=pPSM;
      ++i;
   }
   assert(i == nt);
 
  //mix up the examples
  for(i = 0; i < nt; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nt-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nt-1)); 
      s = all_neg_examples1[p1];
      all_neg_examples1[p1] = all_neg_examples1[p2];
      all_neg_examples1[p2] = s;
    }

  
  //distribute the shuffled set between train and test
  count = 0;
  while(count < kkk){
    train.scores[ix1] = all_neg_examples1[count];
      ++ix1;
      ++count;
  }
  while (count < nt)
    {
      test.scores[ix2] = all_neg_examples1[count];
      ++ix2;    
      count++;
    }
  
  
  train.pos=kk;
  test.pos=ll;
  train.neg=ix1-kk;
  test.neg=ix2-ll;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;

}


void Scores::fillFeaturesSplit(Scores& train,Scores& test,SetHandler& norm,SetHandler& shuff,SetHandler& shuff1,SetHandler& shuff2, const double ratio) {

  PSMDescription * pPSM;
  SetHandler::Iterator shuffIter(&shuff),shuffIter1(&shuff1),shuffIter2(&shuff2), normIter(&norm);
 
  set<string> peptides_train;
  set<string> peptides_test;

  assert(ratio>0 && ratio < 1);
  cerr << "fill features large\n";
  //shuffled
  int n = shuff.getSize();
  int k = (int)(shuff.getSize()*ratio);
  int l = shuff.getSize() - k;
  //shuffled1
  int nt = shuff1.getSize();
  int kkk = (int)(shuff1.getSize()*ratio);
  int lll = shuff1.getSize() - kkk;
  //shuffled2
  int ntt = shuff2.getSize();
  int kkkk = (int)(shuff2.getSize()*ratio);
  int llll = shuff2.getSize() - kkkk;
  //normal
  int nz  = norm.getSize();
  int kk = (int)(norm.getSize()*ratio);
  int ll = norm.getSize() - kk;

  
  ScoreHolder s;
  train.scores.resize(kk+k+kkk+kkkk,s);
  test.scores.resize(ll+l+lll+llll,s);

  //collect everybody positive
  vector<ScoreHolder> all_pos_examples;
  all_pos_examples.resize(nz,s);
  while((pPSM=normIter.getNext())!=NULL) {
    all_pos_examples[i].label=1;
    all_pos_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nz);
  //mix up the examples
  for(i = 0; i < nz; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nz-1)); 
      s = all_pos_examples[p1];
      all_pos_examples[p1] = all_pos_examples[p2];
      all_pos_examples[p2] = s;
    }
  //distribute the normal set between train and test
  int ix1=0;
  int ix2=0;
  int count = 0;
  while (count < kk)
    {
      train.scores[ix1] = all_pos_examples[count];
      ++ix1;
      ++count;
    }
  while (count < nz)
    {
      test.scores[ix2] = all_pos_examples[count];
      ++ix2;    
      count++;
    }
  assert(ix1==kk);
  assert(ix2==ll);


  //collect everybody negative
  i = 0;
  vector<ScoreHolder> all_neg_examples;
  all_neg_examples.resize(n,s);
  while((pPSM=shuffIter.getNext())!=NULL) {
    all_neg_examples[i].label=-1;
    all_neg_examples[i].pPSM=pPSM;
    ++i;
  }
  assert(i == n);
  //mix up the examples
  for(i = 0; i < n; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(n-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(n-1)); 
      s = all_neg_examples[p1];
      all_neg_examples[p1] = all_neg_examples[p2];
      all_neg_examples[p2] = s;
    }
  //distribute the shuffled set between train and test
  count = 0;
  while(count < k){
    train.scores[ix1] = all_neg_examples[count];
    ++ix1;
    ++count;
  }
  while (count < n)
    {
      test.scores[ix2] = all_neg_examples[count];
      ++ix2;    
      count++;
    }

  //collect everybody negative again
  i = 0;
  vector<ScoreHolder> all_neg_examples1;
  all_neg_examples1.resize(nt,s);
  while((pPSM=shuffIter1.getNext())!=NULL) {
    all_neg_examples1[i].label=-1;
    all_neg_examples1[i].pPSM=pPSM;
    ++i;
  }
  assert(i == nt);
  //mix up the examples
  for(i = 0; i < nt; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(nt-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(nt-1)); 
      s = all_neg_examples1[p1];
      all_neg_examples1[p1] = all_neg_examples1[p2];
      all_neg_examples1[p2] = s;
    }
  //distribute the shuffled set between train and test
  count = 0;
  while(count < kkk){
    train.scores[ix1] = all_neg_examples1[count];
    ++ix1;
    ++count;
  }
  while (count < nt)
    {
      test.scores[ix2] = all_neg_examples1[count];
      ++ix2;    
      count++;
    }
  
  //collect everybody negative again
  i = 0;
  vector<ScoreHolder> all_neg_examples2;
  all_neg_examples2.resize(ntt,s);
  while((pPSM=shuffIter2.getNext())!=NULL) {
     all_neg_examples2[i].label=-1;
     all_neg_examples2[i].pPSM=pPSM;
      ++i;
   }
   assert(i == ntt);
   //mix up the examples
   for(i = 0; i < ntt; i++)
    {
      int p1 = (int)((double)rand()/RAND_MAX*(ntt-1)); 
      int p2 = (int)((double)rand()/RAND_MAX*(ntt-1)); 
      s = all_neg_examples2[p1];
      all_neg_examples2[p1] = all_neg_examples2[p2];
      all_neg_examples2[p2] = s;
    }
  //distribute the shuffled set between train and test
  count = 0;
  while(count < kkkk){
    train.scores[ix1] = all_neg_examples2[count];
      ++ix1;
      ++count;
  }
  while (count < ntt)
    {
      test.scores[ix2] = all_neg_examples2[count];
      ++ix2;    
      count++;
    }


  train.pos=kk;
  test.pos=ll;
  train.neg=ix1-kk;
  test.neg=ix2-ll;
  train.factor = train.pos/(double)train.neg;
  test.factor = train.pos/(double)train.neg;


}






*/

