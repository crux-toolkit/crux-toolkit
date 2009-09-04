/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Caller.h,v 1.47.2.4 2009/08/01 20:29:50 jasonw Exp $
 *******************************************************************************/
#ifndef CALLER_H_
#define CALLER_H_

#include <time.h>
#include "SanityCheck.h"

namespace qranker {

class Caller
{
public:
    enum XvType {NO_XV=0, EACH_STEP, WHOLE};
    enum SetHandlerType {NORMAL=0,SHUFFLED,SHUFFLED1,SHUFFLED2};
public:
	Caller();
	virtual ~Caller();
    void readRetentionTime(string filename);
    void step(Scores& train, vector<double>& w, double Cpos, double Cneg, double fdr);
    void train();
    void trainEm(vector<vector<double> >& w);
    int xvalidate_step(vector<vector<double> >& w);
    int xv_step(vector<vector<double> >& w);
    void xvalidate(vector<vector<double> >& w);
	static string greeter();
	string extendedGreeter();
    bool parseOptions(int argc, char **argv);
    void printWeights(ostream & weightStream, vector<double>& w);
    void readWeights(istream & weightStream, vector<double>& w);
    void readFiles(bool &doSingleFile);
    void filelessSetup(unsigned int nsets,const unsigned int numFeatures, const unsigned int numSpectra, char ** fetureNames, double pi0);
    void fillFeatureSets();    
    int preIterationSetup();
    Scores* getFullSet() {return &fullset;}    
    int run();
    SetHandler * getSetHandler(SetHandlerType sh) {
        switch(sh) {
           case NORMAL: return &normal;
           case SHUFFLED: return &shuffled;
           case SHUFFLED1: return &shuffled1;
           case SHUFFLED2: return &shuffled2;
           default: return NULL;
        }
    }

    void train_net_two(Scores &set);
    void train_many_general_nets();
    void train_many_target_nets_ave();
    void train_many_nets();
    void three_plot(double qv);
    
    int getOverFDR(Scores &set, NeuralNet &n, double fdr);
    void getMultiFDR(Scores &set, NeuralNet &n, vector<double> &qval);
    void printNetResults(vector<int> &scores);
    void write_max_nets(string filename, NeuralNet *max_net);
    
    void xvalidate_net(double qv);
    void train_general_net(Scores &train, Scores &thresh, double qv);
    void train_target_net(Scores &train, Scores &thresh, double qv);

protected:
    Normalizer * pNorm;
    SanityCheck * pCheck;
    AlgIn *svmInput;
    string modifiedFN;
    string modifiedDecoyFN;
    string forwardFN;
    string decoyFN;
    string decoyWC;
    string rocFN;
    string gistFN;
    string tabFN;
    string weightFN;
    string call;
    string spectrumFile;
    bool gistInput;
    bool tabInput;
    bool dtaSelect;
    bool docFeatures;
    bool reportPerformanceEachIteration;
    double test_fdr;
    double selectionfdr;
    double selectedCpos;
    double selectedCneg;
    double threshTestRatio;    
    double trainRatio;    
    unsigned int niter;
    unsigned int seed;
    time_t startTime;
    clock_t startClock;
    const static unsigned int xval_fold;
    XvType xv_type; 
    vector<Scores> xv_train,xv_test;
    vector<double> xv_cposs,xv_cfracs;
    SetHandler normal,shuffled;
    Scores fullset; 
    map<int,double> scan2rt; 

    
    string decoyFN1,decoyFN2;
    bool shuffled1_present, shuffled2_present;
    SetHandler shuffled1,shuffled2;



    Scores trainset,testset,thresholdset;
    int trainNN;
    int splitPositives;

    NeuralNet net;
    int num_hu;
    double mu;
    double weightDecay;

    int ind_low;
    int interval;
    int switch_iter;
    double loss;

    int num_qvals;
    vector<double> qvals;
    vector<double> qvals1;
    vector<double> qvals2;
    
    vector<int> overFDRmulti;
    vector<int> max_overFDR;
    vector<int> ave_overFDR;
    NeuralNet* max_net_gen;
    NeuralNet* max_net_targ;


};

} // qranker namspace

#endif /*CALLER_H_*/
