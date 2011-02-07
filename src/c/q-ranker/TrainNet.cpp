#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <math.h>
using namespace std;
#include "Option.h"
#include "SanityCheck.h"
#include "SqtSanityCheck.h"
#include "DataSet.h"
#include "IntraSetRelation.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SetHandler.h"
#include "Net.h"
#include "ssl.h"
#include "Caller.h"
#include "Globals.h"
#include "MSReader.h"
#include "Spectrum.h"



namespace qranker {


string res_prefix="yeast_trypsin";





/********************* writing out results functions***************************************/
int Caller :: getOverFDR(Scores &set, NeuralNet &n, double fdr)
{
  double r = 0.0;
  const double* featVec;
  int label = 0;
  loss = 0;
  for(unsigned int i = 0; i < set.size(); i++)
    {
      featVec = set[i].pPSM->features;
      label = set[i].label;
      r = n.classify(featVec);
      set[i].score = r;
      set[i].pPSM->sc=r;
      loss += n.get_err(label);
    }
  return set.calcOverFDR(fdr);
}
 

void Caller :: getMultiFDR(Scores &set, NeuralNet &n, vector<double> &qvalues)
{
  double r = 0.0;
  const double* featVec;
  int label = 0;
  loss = 0;
 
  for(unsigned int i = 0; i < set.size(); i++)
    {
      featVec = set[i].pPSM->features;
      label = set[i].label;
      r = n.classify(featVec);
      set[i].score = r;
      loss += n.get_err(label);
     }
 

  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti);
}


void Caller :: printNetResults(vector<int> &scores)
{
  double qv; int fdr;
  cerr << "QVALS SCORES:: ";
  for(unsigned int count = 0; count < qvals.size();count++)
    {  
      qv=qvals[count]; 
      fdr = scores[count];
      cerr << qv << ":" << fdr << " ";
    }
  //cerr << "loss " << loss << " ";
  cerr << endl;
}

void Caller :: write_max_nets(string filename, NeuralNet* max_net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  
  for(int count = 0; count < num_qvals; count++)
    {
      net = max_net[count];
      int r = getOverFDR(testset,net, qvals[count]);
      int r1 = getOverFDR(trainset,net, qvals[count]);
      f1 << qvals[count] << " " << r1 << " " << r << "\n";
      
      double qn;
      if(qvals[count] < 0.01)
	qn = 0.0012;
      else
	qn = 0.005;
      r = getOverFDR(testset,net, qvals[count]+qn);
      r1 = getOverFDR(trainset,net, qvals[count]+qn);
      f1 << qvals[count]+qn << " " << r1 << " " << r << "\n";
      
      //write out the net for a certain qvals
      ostringstream s1_net;
      s1_net << filename <<  "_net_" << qvals[count] << ".txt";
      ofstream f1_net(s1_net.str().c_str());
      net.write_to_file(f1_net);
      f1_net.close();
    }
  f1.close();
}


/********** xvalidation functions********************************************/

void Caller :: train_general_net(Scores &train, Scores &thresh, double qv)
{
  int overFDR = 0;
  int best_overFDR = 0;
  
  overFDR = getOverFDR(thresh,net, qv);
  ind_low = -1;

  NeuralNet best_net;
  best_net = net;
  best_overFDR = overFDR;
  
  for(int i=0;i<switch_iter;i++) {
    train_net_two(train);
    
    //record the best result
    overFDR = getOverFDR(thresh,net, qv);
    if(overFDR > best_overFDR)
      {
	best_overFDR = overFDR;
	best_net = net;
      }
      
    /*
    if((i % 100) == 0)
      {
	if(VERB>1) cerr << "Iteration " << i << " : \n";
	cerr << "train: ";
	cerr << net.get_mu() << " " << qv << " " << getOverFDR(train,net,qv) << "\n";
	cerr << "thresh: ";
	cerr << net.get_mu() << " " << qv << " " << getOverFDR(thresh,net,qv) << "\n";
	cerr << "testset: ";
	cerr << net.get_mu() << " " << qv << " " << getOverFDR(testset,net,qv) << "\n";
      }
    */
    
  }
  net = best_net;
  best_net.clear();
}


void Caller :: train_target_net(Scores &train, Scores &thresh, double qv)
{
  
  int overFDR = 0;
  int best_overFDR = 0;
  
  overFDR = getOverFDR(train,net, qv);
  ind_low = overFDR;
  
  NeuralNet best_net;
  best_net = net;
  best_overFDR = overFDR;
  //cerr << "target_q " << qv << " starting " << overFDR << "\n";
 
  for(int i = switch_iter; i < niter; i++)
    {
      overFDR = getOverFDR(train,net,qv);
      train_net_two(train);
      
      overFDR = getOverFDR(thresh,net,qv);
      if(overFDR > best_overFDR)
	{
	  best_overFDR = overFDR;
	  best_net = net;
	}
      /*
      if(i%100 == 0)
	{
	  cerr << i << " " << interval << "\n";
	  //cerr << "trainset: ";
	  //cerr << net.get_mu() << " " << qv << " " << getOverFDR(train,net,qv) << "\n";
	  cerr << "thresh: ";
	  cerr << net.get_mu() << " " << qv << " " << getOverFDR(thresh,net,qv) << "\n";
	  cerr << "testset: ";
	  cerr << net.get_mu() << " " << qv << " " << getOverFDR(testset,net,qv) << "\n";
	}
      */
    }
  //cerr << "end on interval " << interval << "\n";
  net = best_net;
  best_net.clear();
}

 
/**
 * Do cross-validation to select two hyperparameters: the learning
 * rate and the weight decay.  The optimization criterion is the
 * number of target scores below a specified q-value threshold.
 */
void Caller :: xvalidate_net(
  double qv ////< The q-value threshold -in
)
{
  
  vector <double> xv_wt;
  xv_wt.push_back(0.00001);xv_wt.push_back(0.00005);
  
  vector <double> xv_mu;
  //xv_mu.push_back(0.005);xv_mu.push_back(0.007);xv_mu.push_back(0.01);
  xv_mu.push_back(0.005);xv_mu.push_back(0.01);

  

  //split the trainset
  //for(unsigned int i = 0; i < xval_fold; i++)
   //cerr << "size " << trainset.size() << " " << xv_train[i].size() << " " << xv_test[i].size() << "\n";
  
  //set the linear flag: 1 if linear, 0 otherwise
  int lf = 0;
  //set whether there is bias in the linear units: 1 if yes, 0 otherwise
  int bs = 1;
  //cost linear flag indicating whether to use the sigmoid(0) or linear loss(1)
  int clf = 0;
 
  double best_sum = 0;
  int best_wt = 0;
  int best_mu = 0;
  cerr << "xvalidating \n";
  for(unsigned int h = 0; h < xv_wt.size();h++)
    {
      for(unsigned int m = 0; m < xv_mu.size(); m++)
	{
	  mu = xv_mu[m];
	  int overFDR = 0;
	  double sm=0;
	  for(unsigned int i = 0; i < xval_fold; i++)
	    {
	      net.initialize(FeatureNames::getNumFeatures(),num_hu,mu,clf,lf,bs);
	      net.set_weightDecay(xv_wt[h]);
	      cerr << "learning rate " << mu << " weight decay "<<xv_wt[h] << " xval set " << i << " size xval set " << xv_train[i].size() << "\n";
	      train_general_net(xv_train[i], xv_train[i],qv);

	      net.set_cost_flag(1);;
	      net.remove_bias();
	     
	      train_target_net(xv_train[i],xv_train[i],qv);
	      overFDR = getOverFDR(xv_test[i],net, qv);
	      sm +=overFDR;
	      net.clear();
	    }
	  sm /=xval_fold;
	  //cerr << "best_sum " << best_sum << " sum " << sm << "\n";
	  if(sm > best_sum)
	    {
	      best_sum = sm;
	      best_mu = m;
	      best_wt = h;
	    }
	}
    }
  mu = xv_mu[best_mu];
  weightDecay = xv_wt[best_wt];    
  
  cerr << "training target: mu " << best_mu << " " << mu << " wt_decay " << weightDecay << "\n";
   
}




/*********************** training net functions*******************************************/

void Caller :: train_net_two(Scores &set)
{
  double r1 = 0;
  double r2 = 0;
  double diff = 0;
  int label = 1;
  const double *featVec;
  
  for(unsigned int i = 0; i < set.size(); i++)
    { 
      if(ind_low<0)
	{
	  unsigned int ind;
	  ind = random()%set.size();
	  featVec = set[ind].pPSM->features;
	  label = set[ind].label;
	  net.train(featVec,label);
	}
      else
	{
	  interval = ind_low*2;
	  unsigned int ind1, ind2;
	  int label_flag = 1;
	  //get the first index
	  if(interval == 0)
	    ind1 = 0;
	  else
	    ind1 = random()%(interval);
	  if(ind1<0 || ind1>set.size()-1) continue;
	  if(set[ind1].label == 1)
	    label_flag = -1;
	  else
	    label_flag = 1;
	  
	  int cn = 0;
	  while(1)
	    {
	      if(interval == 0)
		ind2 = random()%set.size();
	      else
		ind2 = random()%(interval);
	      if(ind2<0 || ind2>set.size()-1) continue;
	      if(set[ind2].label == label_flag) break;
	      if(cn > 1000)
		{
		  ind2 = random()%set.size();
		  break;
	      	}
	      cn++;
	    }
	  
	  
	  //pass both through the net
	  r1 = net.classify(set[ind1].pPSM->features);
	  r2 = net.classify(set[ind2].pPSM->features);
	  diff = r1-r2;


	  label=0;

	  if(  set[ind1].label==1 && set[ind2].label==-1)
	    label=1;
	  
	  if( set[ind1].label==-1 && set[ind2].label==1)
	    label=-1;

	  
	  if(label != 0)
	    {
	      if(label*diff<1)
		{
		  net.train1(set[ind1].pPSM->features, set[ind1].label);
		  net.train1(set[ind2].pPSM->features, set[ind2].label);
		}
	    }
	}
    }
}




void Caller :: train_many_general_nets()
{
  ind_low = -1;
  for(int i=0;i<switch_iter;i++) {
    train_net_two(trainset);
    
    //record the best result
    getMultiFDR(thresholdset,net,qvals);
    for(int count = 0; count < num_qvals;count++)
      {
	if(overFDRmulti[count] > max_overFDR[count])
	  {
	    max_overFDR[count] = overFDRmulti[count];
	    max_net_gen[count] = net;
	  }
      }
    if((i % 30) == 0)
      {
	
	cerr << "Iteration " << i << " : \n";
	cerr << "trainset: ";
	getMultiFDR(trainset,net,qvals);
	printNetResults(overFDRmulti);
	cerr << "\n";
	cerr << "testset: ";
	getMultiFDR(testset,net,qvals);
	printNetResults(overFDRmulti);
	cerr << "\n";
      }
  }

}

void Caller :: train_many_target_nets_ave()
{

  int  thr_count = num_qvals-1;
  while (thr_count > 0)
    {
      net = max_net_gen[thr_count];
      net.set_cost_flag(1);;
      net.remove_bias();
           
      //cerr << "training thresh " << thr_count << " mu " << net.get_mu() << "\n";
      cerr << "training thresh " << thr_count  << "\n";
      
      ind_low = getOverFDR(trainset, net, qvals[thr_count]);
      for(int i=switch_iter;i<niter;i++) {
		
	//sorts the examples in the training set according to the current net scores
	getMultiFDR(trainset,net,qvals);
	train_net_two(trainset);
	
	//record best result according to the average around q-value of interest
	for(int count = 0; count < num_qvals; count++)
	  ave_overFDR[count] = 0;
	//record the best result
	getMultiFDR(thresholdset,net, qvals);
	for(int count = 0; count < num_qvals; count++)
	  ave_overFDR[count] += overFDRmulti[count];
		
	getMultiFDR(thresholdset,net, qvals1);
	for(int count = 0; count < num_qvals; count++)
	  ave_overFDR[count] += overFDRmulti[count];
	
	getMultiFDR(thresholdset,net, qvals2);
	for(int count = 0; count < num_qvals; count++)
	  ave_overFDR[count] += overFDRmulti[count];
	
	for(int count = 0; count < num_qvals; count++)
	  ave_overFDR[count] /=3;
	  
	for(int count = 0; count < num_qvals;count++)
	  {
	    if(ave_overFDR[count] > max_overFDR[count])
	      {
		max_overFDR[count] = ave_overFDR[count];
		max_net_targ[count] = net;
	      }
	  }

	if((i % 50) == 0)
	  {
	    cerr << "Iteration " << i << " : \n";
	    getMultiFDR(testset,net,qvals);
	    cerr << "testset: ";
	    printNetResults(overFDRmulti);
	    cerr << "\n";
	    //cerr << "ave ";
	    //printNetResults(ave_overFDR);
	    //cerr << "\n";
	  }

      }
      thr_count -= 2;
    }
}


void Caller::train_many_nets(
  bool do_xval ////< Select hyperparameters via cross-validation. -in
  ) 
{
  
  switch_iter =100;
  niter = 200;
   
  num_qvals = 14;
  qvals.resize(num_qvals,0.0);
  qvals1.resize(num_qvals,0.0);
  qvals2.resize(num_qvals,0.0);
 
  overFDRmulti.resize(num_qvals,0);
  ave_overFDR.resize(num_qvals,0);
  max_overFDR.resize(num_qvals,0); 
  max_net_gen = new NeuralNet[num_qvals];
  max_net_targ = new NeuralNet[num_qvals];
  
  double q = 0.0;
  for(int count = 0; count < num_qvals; count++)
    {
      qvals[count] = q;
      if (count < 2)
	qvals1[count] = q;
      else
	qvals1[count] = q-0.005;
      qvals2[count] = q+0.005;

      if(q < 0.01)
	q+=0.0025;
      else
	q+=0.01;
    }
  
  //set the linear flag: 1 if linear, 0 otherwise
  int lf = 0;
  if(num_hu == 1)
    lf = 1;
  //set whether there is bias in the linear units: 1 if yes, 0 otherwise
  int bs = 1;
  //cost linear flag indicating whether to use the sigmoid(0) or linear loss(1)
  int clf = 0;

  if (do_xval) {
    xvalidate_net(selectionfdr);
  }

  net.initialize(FeatureNames::getNumFeatures(),num_hu,mu,clf,lf,bs);
  net.set_weightDecay(weightDecay);
  for(int count = 0; count < num_qvals; count++){
    max_net_gen[count] = net;
  }
  
  cerr << "Before iterating\n";
  cerr << "trainset: ";
  getMultiFDR(trainset,net,qvals);
  printNetResults(overFDRmulti);
  cerr << "\n";
  cerr << "testset: ";
  getMultiFDR(testset,net,qvals);
  printNetResults(overFDRmulti);
  cerr << "\n";
  

  train_many_general_nets();
  

  //write out the results of the general net
  if (0){
    cerr << "general net results: ";
    ostringstream filename;
    filename << res_prefix << "_hu" << num_hu;
    write_max_nets(filename.str(), max_net_gen);
  }
  //copy the general net into target nets;
  for(int count = 0; count < num_qvals; count++)
    max_net_targ[count] = max_net_gen[count];

  //calculate average of target q-values
  //and q-value +0.05 and -0.05
  getMultiFDR(thresholdset,net, qvals);
  for(int count = 0; count < num_qvals; count++)
    ave_overFDR[count] += overFDRmulti[count];
  getMultiFDR(thresholdset,net, qvals1);
  for(int count = 0; count < num_qvals; count++)
    ave_overFDR[count] += overFDRmulti[count];
  getMultiFDR(thresholdset,net, qvals2);
  for(int count = 0; count < num_qvals; count++)
    ave_overFDR[count] += overFDRmulti[count];
  for(int count = 0; count < num_qvals; count++)
    ave_overFDR[count] /=3;
  //printNetResults(ave_overFDR);
  for(int count = 0; count < num_qvals;count++)
    max_overFDR[count] = ave_overFDR[count];
  

  train_many_target_nets_ave();  
 
  //write out the results of the target net
  //cerr << "target net results: ";
  //ostringstream s2;
  //s2 << res_prefix << "_hu" << num_hu  << "_targ";
  //write_max_nets(s2.str(), max_net_targ);

  //choose the best net for the selectionfdr
  int max_fdr = 0;
  int fdr = 0;
  int ind = 0;
  for(unsigned int count = 0; count < qvals.size();count++)
    {
      fdr = getOverFDR(thresholdset, max_net_targ[count], selectionfdr);
      if(fdr > max_fdr)
	{
	  max_fdr = fdr;
	  ind = count;
	}
    }
  net = max_net_targ[ind];
  
  cerr << " Found " << getOverFDR(fullset, net, selectionfdr) << " over q<" << selectionfdr << "\n";
  
  delete [] max_net_gen;
  delete [] max_net_targ;

}

} // qranker namspace

