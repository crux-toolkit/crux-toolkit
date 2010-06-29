/*******************************************************************************
 * Percolator unofficial version
 * Copyright (c) 2006-8 University of Washington. All rights reserved.
 * Written by Lukas Käll (lukall@u.washington.edu) in the 
 * Department of Genome Science at the University of Washington. 
 *
 * $Id: Caller.cpp,v 1.97.2.6 2009/08/01 20:29:50 jasonw Exp $
 *******************************************************************************/
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

const unsigned int Caller::xval_fold = 3;

Caller::Caller() : pNorm(NULL), pCheck(NULL), svmInput(NULL),
  modifiedFN(""), modifiedDecoyFN(""), forwardFN(""), decoyFN (""), //shuffledThresholdFN(""), shuffledTestFN(""),
  decoyWC(""), rocFN(""), gistFN(""), tabFN(""), weightFN(""),
  gistInput(false), tabInput(false), dtaSelect(false), docFeatures(false), reportPerformanceEachIteration(false),
  test_fdr(0.01), selectionfdr(0.01), selectedCpos(0), selectedCneg(0), threshTestRatio(0.3), trainRatio(0.5),
		   niter(10), seed(0), xv_type(EACH_STEP), trainNN(1), splitPositives(1), num_hu(5),mu(0.005),weightDecay(0.00005),
		   shuffled1_present(false),shuffled2_present(false)
{
}

Caller::~Caller()
{
    if (pNorm)
      delete pNorm;
    pNorm=NULL;
    if (pCheck)
      delete pCheck;
    pCheck=NULL;
    if (svmInput)
      delete svmInput;
    svmInput=NULL;
}

string Caller::extendedGreeter() {
  ostringstream oss;
  char * host = getenv("HOST");
  if (!host)
    host="unknown_host";
  oss << greeter();
  oss << "Issued command:" << endl << call;
  oss << "Started " << ctime(&startTime);
  oss.seekp(-1, ios_base::cur);
  oss << " on " << host << endl;
  oss << "Hyperparameters fdr=" << selectionfdr;
  oss << ", Cpos=" << selectedCpos << ", Cneg=" << selectedCneg << ", maxNiter=" << niter << endl;
  return oss.str();
}

string Caller::greeter() {
  ostringstream oss;
  oss << "Percolator unofficial version, ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2006-8 University of Washington. All rights reserved." << endl;
  oss << "Written by Lukas Käll (lukall@u.washington.edu) in the" << endl; 
  oss << "Department of Genome Science at the University of Washington." << endl; 
  return oss.str();
}

bool Caller::parseOptions(int argc, char **argv){
  ostringstream callStream;
  callStream << argv[0];
  for (int i=1;i<argc;i++) callStream << " " << argv[i];
  callStream << endl;
  call = callStream.str();
  ostringstream intro,endnote;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   percolator [options] normal shuffle [[shuffled_treshhold] shuffled_test]" << endl;
  intro << "or percolator [options] -P pattern normal_and_shuffled.sqt" << endl;
  intro << "or percolator [options] -g gist.data gist.label" << endl << endl;
  intro << "   where normal is the normal sqt-file," << endl;
  intro << "         shuffle the shuffled sqt-file used in the training," << endl;
  intro << "         shuffle_test is an otional second shuffled sqt-file for q-value calculation" << endl;
  intro << "         shuffle_treshhold is an otional shuffled sqt-file for determine q-value treshold" << endl << endl;
  intro << "To be able to merge small data set one can replace the sqt-files with meta" << endl;
  intro << "files. Meta files are text files containing the paths of sqt-files, one path" << endl;
  intro << "per line. For successful result, the different runs should be generated under" << endl;
  intro << "similair condition. Particulary, they need to be generated with the same protease." << endl;
  // init
  CommandLineParser cmd(intro.str());
  cmd.defineOption("o","sqt-out",
    "Create an SQT file with the specified name from the given target SQT file, \
replacing the XCorr value the learned score and Sp with the negated q-value.","filename");
  cmd.defineOption("s","shuffled",
    "Same as -o, but for the decoy SQT file",
    "filename");
  cmd.defineOption("P","pattern",
    "Option for single SQT file mode defining the name pattern used for shuffled data base. \
Typically set to random_seq","pattern");
  cmd.defineOption("p","Cpos",
    "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.",
    "value");
  cmd.defineOption("n","Cneg",
    "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified or -p not specified.",
    "value");
  cmd.defineOption("F","trainFDR",
    "False discovery rate threshold to define positive examples in training. Set by cross validation if 0. Default is 0.01.",
    "value");
  cmd.defineOption("t","testFDR",
    "False discovery rate threshold for evaluating best cross validation result and the reported end result. Default is 0.01.",
    "value");
  cmd.defineOption("i","maxiter",
    "Maximal number of iteratins",
    "number");
  cmd.defineOption("m","matches",
    "Maximal number of matches to take in consideration per spectrum when using sqt-files",
    "number");
  cmd.defineOption("f","train-ratio",
    "Fraction of the negative data set to be used as train set when only providing one negative set, remaining examples will be used as test set. Set to 0.6 by default.",
    "value");
  cmd.defineOption("G","gist-out",
    "Output the computed features to the given file in tab-delimited format. A file with the features, named <trunc name>.data, and a file with the labels named <trunc name>.label will be created",
    "trunc name");
  cmd.defineOption("g","gist-in",
    "Input files are given as gist files. In this case first argument should be a file name \
of the data file, the second the label file. Labels are interpreted as 1 -- positive train \
and test set, -1 -- negative train set, -2 -- negative in test set.","",TRUE_IF_SET);
  cmd.defineOption("J","tab-out",
    "Output the computed features to the given file in tab-delimited format. A file with the features with the given file name will be created",
    "file name");
  cmd.defineOption("j","tab-in",
    "Input files are given as a tab delimetered file. In this case the only argument should be a file name\
of the data file. The tab delimeterad fields should be id <tab> label <tab> feature1 \
<tab> ... <tab> featureN <tab> peptide <tab> proteinId1 <tab> .. <tab> proteinIdM \
Labels are interpreted as 1 -- positive train \
and test set, -1 -- negative train set, -2 -- negative in test set.","",TRUE_IF_SET);
  cmd.defineOption("w","weights",
    "Output final weights to the given file",
    "filename");
  cmd.defineOption("W","init-weights",
    "Read initial weights from the given file",
    "filename");
  cmd.defineOption("v","verbose",
    "Set verbosity of output: 0=no processing info, 5=all, default is 2",
    "level");
  cmd.defineOption("r","result",
    "Output result file (score ranked labels) to given filename",
    "filename");
  cmd.defineOption("u","unitnorm",
    "Use unit normalization [0-1] instead of standard deviation normalization","",TRUE_IF_SET);
  cmd.defineOption("a","aa-freq","Calculate amino acid frequency features","",TRUE_IF_SET);
  cmd.defineOption("b","PTM","Calculate feature for number of post-translational modifications","",TRUE_IF_SET);
  cmd.defineOption("d","DTASelect",
    "Add an extra hit to each spectra when writing sqt files","",TRUE_IF_SET);
  cmd.defineOption("R","test-each-itteration","Measure performance on test set each itteration","",TRUE_IF_SET);
  cmd.defineOption("Q","quadratic",
    "Calculate quadratic feature terms","",TRUE_IF_SET);
  cmd.defineOption("O","override",
    "Overide error check and do not fall back on default score vector in case of suspect score vector",
    "",TRUE_IF_SET);
  cmd.defineOption("I","intra-set",
    "Turn Off calculation of intra set features","",TRUE_IF_SET);
  cmd.defineOption("y","notryptic",
    "Turn off calculation of tryptic/chymo-tryptic features.","",TRUE_IF_SET);
  cmd.defineOption("c","chymo",
    "Replace tryptic features with chymo-tryptic features.","",TRUE_IF_SET);
  cmd.defineOption("e","elastase",
    "Replace tryptic features with elastase features.","",TRUE_IF_SET);
  cmd.defineOption("x","whole-xval",
    "Select hyper parameter cross validation to be performed on whole itterating procedure, rather than on each iteration step."
    ,"",TRUE_IF_SET);
  cmd.defineOption("q","pi-zero",
    "Estimated proportion of PSMs generated by a random model. Default value is 1-0.1/<value of -m option>."
    ,"value");
  cmd.defineOption("S","seed",
    "Setting seed of the random number generator. Default value is 0"
    ,"value");
  cmd.defineOption("2","ms2-file",
    "File containing spectra and retention time. The file could be in mzXML, MS2 or compressed MS2 file.",
    "filename");
  cmd.defineOption("M","isotope",
    "Mass difference calculated to closest isotope mass rather than to the average mass.","",TRUE_IF_SET);
  cmd.defineOption("D","doc",
    "Include description of correct features.","",TRUE_IF_SET);


  cmd.defineOption("N","trainNet",
		   "Use neural net for training","",TRUE_IF_SET);
  cmd.defineOption("hu","numHU",
		   "Number of hidden units in neural net",
		   "number");
  
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);

  
  // now query the parsing results
  if (cmd.optionSet("o"))
    modifiedFN = cmd.options["o"];
  if (cmd.optionSet("s"))
    modifiedDecoyFN = cmd.options["s"];
  if (cmd.optionSet("P"))
    decoyWC = cmd.options["P"];
  if (cmd.optionSet("p")) {
    selectedCpos = cmd.getDouble("p",0.0,1e127);
  }
  if (cmd.optionSet("n")) {
    selectedCneg = cmd.getDouble("n",0.0,1e127);
  }
  if (cmd.optionSet("G"))
    gistFN = cmd.options["G"];
  if (cmd.optionSet("g")) {
    gistInput=true;
    if (cmd.arguments.size()!=2) {
      cerr << "Provide exactly two arguments when using gist input" << endl;
      exit(-1);
    }   
  }
  if (cmd.optionSet("J"))
    tabFN = cmd.options["J"];
  if (cmd.optionSet("j")) {
    tabInput=true;
    if (cmd.arguments.size()!=1) {
      cerr << "Provide exactly one arguments when using tab delimetered input" << endl;
      exit(-1);
    }   
  }
  if (cmd.optionSet("w"))
    weightFN = cmd.options["w"];
  if (cmd.optionSet("W"))
    SanityCheck::setInitWeightFN(cmd.options["W"]);
  if (cmd.optionSet("f")) {
    double frac = cmd.getDouble("f", 0.0 ,1.0);
    trainRatio=frac;
  }
  if (cmd.optionSet("r"))
    rocFN = cmd.options["r"];
  if (cmd.optionSet("u"))
    Normalizer::setType(Normalizer::UNI);
  if (cmd.optionSet("d"))
    dtaSelect=true;
  if (cmd.optionSet("Q"))
    DataSet::setQuadraticFeatures(true);
  if (cmd.optionSet("O"))
    SanityCheck::setOverrule(true);
  if (cmd.optionSet("I"))
    DataSet::setCalcIntraSetFeatures(false);
  if (cmd.optionSet("y"))
    DataSet::setEnzyme(NO_ENZYME);
  if (cmd.optionSet("R"))
    reportPerformanceEachIteration=true;
  if (cmd.optionSet("e"))
    DataSet::setEnzyme(ELASTASE);
  if (cmd.optionSet("c"))
    DataSet::setEnzyme(CHYMOTRYPSIN);
  if (cmd.optionSet("a"))
    DataSet::setAAFreqencies(true);
  if (cmd.optionSet("b"))
    DataSet::setPTMfeature(true);
  if (cmd.optionSet("x"))
    xv_type=WHOLE;
  if (cmd.optionSet("i")) {
    niter = cmd.getInt("i",0,100000000);
  }
  if (cmd.optionSet("m")) {
    int m = cmd.getInt("m",1,30000); 
    DataSet::setHitsPerSpectrum(m);
    Scores::pi0 = 1 - (1-Scores::pi0)/(double)m;
  }
  if (cmd.optionSet("q")) {
    Scores::pi0 = cmd.getDouble("q",0.0,1.0);
  }
  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v",0,10));
  }
  if (cmd.optionSet("F")) {
    selectionfdr = cmd.getDouble("F",0.0,1.0);
  }
  if (cmd.optionSet("t")) {
    test_fdr = cmd.getDouble("t",0.0,1.0);
  }
  if (cmd.optionSet("S")) {
    seed = cmd.getInt("S",0,20000);
  }
  if (cmd.optionSet("2")) {
    spectrumFile = cmd.options["2"];
  }
  if (cmd.optionSet("M"))
    DescriptionOfCorrect::setIsotopeMass(true);
  if (cmd.optionSet("D")) {
    docFeatures = true;
    DataSet::setCalcDoc(true);
  }


  if (cmd.optionSet("N"))    
    trainNN=true;
  if (cmd.optionSet("hu")) {
    num_hu = cmd.getInt("hu",0,50);
  }

  
  if (cmd.arguments.size()>4) {
      cerr << "Too many arguments given" << endl;
      cmd.help();
  }
  if (cmd.arguments.size()==0) {
      cerr << "No arguments given" << endl;
      cmd.help();
  }
  if (cmd.arguments.size()>0)
     forwardFN = cmd.arguments[0];
  if (cmd.arguments.size()>1)
     decoyFN = cmd.arguments[1];
  if (cmd.arguments.size()>2)
    {
      decoyFN1 = cmd.arguments[2];
      shuffled1_present = true;
    }
  if (cmd.arguments.size()>3)
    {
      decoyFN2 = cmd.arguments[3];
      shuffled2_present = true;
    }
  return true;
}




void Caller::readRetentionTime(string filename) {
  MSToolkit::MSReader r;
  MSToolkit::Spectrum s;

  r.setFilter(MSToolkit::MS2);
  
  char* cstr = new char [filename.size()+1];
  strcpy (cstr, filename.c_str());
  r.readFile(cstr,s);

  while (s.getScanNumber()!=0){
    scan2rt[s.getScanNumber()] = (double) s.getRTime();
    r.readFile(NULL,s);
  }
  delete [] cstr;
}

void Caller::filelessSetup(unsigned int nsets, const unsigned int numFeatures, int* numSpectra, char ** featureNames, double pi0) {

  pCheck = new SanityCheck();
  normal.filelessSetup(numFeatures, numSpectra[0],1);
  shuffled.filelessSetup(numFeatures, numSpectra[1],-1);
  shuffled1_present=false;
  shuffled2_present=false;
  if (nsets > 2)
    {
      shuffled1.filelessSetup(numFeatures, numSpectra[2],-1);
      shuffled1_present = true;
    }
  if (nsets > 3)
    {
      shuffled2.filelessSetup(numFeatures, numSpectra[3],-1);
      shuffled2_present = true;
    }
  Scores::pi0 = pi0;
  for (unsigned int ix=0;ix<numFeatures;ix++){
    string fn = featureNames[ix];
    DataSet::getFeatureNames().insertFeature(fn);
  }
}

void Caller::readFiles(bool &doSingleFile) {
  if (gistInput) {
    pCheck = new SanityCheck();
    normal.readGist(forwardFN,decoyFN,1);
    shuffled.readGist(forwardFN,decoyFN,-1);
  } else if (tabInput) {
    pCheck = new SanityCheck();
    normal.readTab(forwardFN,1);
    shuffled.readTab(forwardFN,-1);
  } else if (!doSingleFile) {
    pCheck = new SqtSanityCheck();
    normal.readFile(forwardFN,1);
    shuffled.readFile(decoyFN,-1);
    if(shuffled1_present)
      shuffled1.readFile(decoyFN1,-1);
    if(shuffled2_present)
      shuffled2.readFile(decoyFN2,-1);   
  } else {
    pCheck = new SqtSanityCheck();
    normal.readFile(forwardFN,decoyWC,false);  
    shuffled.readFile(forwardFN,decoyWC,true);  
  }
  if (spectrumFile.size()>0)
      readRetentionTime(spectrumFile);
}


void Caller::train() {
  train_many_nets();
}


void Caller::fillFeatureSets() {

    if(shuffled2_present)
      Scores::fillFeaturesSplit(trainset,testset,normal,shuffled,shuffled1,shuffled2,trainRatio);
    else if( shuffled1_present)
      Scores::fillFeaturesSplit(trainset,testset,normal,shuffled,shuffled1,trainRatio);     
    else
    Scores::fillFeaturesSplit(trainset,testset,normal,shuffled,trainRatio);     
  thresholdset=trainset;
  
  fullset.fillFeatures(normal,shuffled);

  //if (VERB>0) {
    cerr << "Train set contains " << trainset.posSize() << " positives and " << trainset.negSize() << " negatives, size ratio="
         << trainset.factor << " and pi0=" << trainset.pi0 << endl;
    cerr << "Test set contains " << testset.posSize() << " positives and " << testset.negSize() << " negatives, size ratio="
         << testset.factor << " and pi0=" << testset.pi0 << endl;
    //}
  
  set<DataSet *> all;
  all.insert(normal.getSubsets().begin(),normal.getSubsets().end());
  all.insert(shuffled.getSubsets().begin(),shuffled.getSubsets().end());
  if (shuffled1_present)
    all.insert(shuffled1.getSubsets().begin(),shuffled1.getSubsets().end());
  if (shuffled2_present)
    all.insert(shuffled2.getSubsets().begin(),shuffled2.getSubsets().end());
  pNorm=Normalizer::getNormalizer();
  pNorm->setSet(all,FeatureNames::getNumFeatures(),docFeatures?DescriptionOfCorrect::totalNumRTFeatures():0);
  pNorm->normalizeSet(all);
  if (docFeatures) {
    for (set<DataSet *>::iterator myset=all.begin();myset!=all.end();++myset)
      (*myset)->setRetentionTime(scan2rt);
  }
}




int Caller::preIterationSetup() {
     trainset.createXvalSets(xv_train,xv_test,xval_fold);
}    

int Caller::run() {
  srand(seed);
  if(VERB>0)  cerr << extendedGreeter();
  //File reading
  bool doSingleFile = !decoyWC.empty();
  readFiles(doSingleFile);
  fillFeatureSets();
  preIterationSetup();
  train();
  
  cerr << " Found " << getOverFDR(fullset, net, selectionfdr) << " over q<" << selectionfdr << "\n";
  normal.print(fullset);
  
  return 0;
}

} // qranker namspace

