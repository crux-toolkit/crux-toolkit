#ifndef QRANKER_H_
#define QRANKER_H_
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
#ifdef CRUX
#include "CruxApplication.h"
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#endif
#include "analyze_psms.h"
#include "PepXMLWriter.h"
#include "NeuralNet.h"
#include "DataSet.h"
#include "PSMScores.h"
#include "SQTParser.h"
#include "CruxParser.h"
#include "TabDelimParser.h"
#include "CruxApplication.h"
#include "mass.h"
#include "objects.h"

class QRanker: public CruxApplication
{

public:
  QRanker();
  virtual ~QRanker();

  int run();
  void train_net_sigmoid(PSMScores &set, int interval);
  void train_net_ranking(PSMScores &set, int interval);
  void train_net_hinge(PSMScores &set, int interval);
  void count_pairs(PSMScores &set, int interval);
  void train_many_general_nets();
  void train_many_target_nets();
  void train_many_nets();
    
  int getOverFDR(PSMScores &set, NeuralNet &n, double fdr);
  void getMultiFDR(PSMScores &set, NeuralNet &n, vector<double> &qval);
  void getMultiFDRXCorr(PSMScores &set, vector<double> &qval);
  void printNetResults(vector<int> &scores);
  void write_results();
  void write_results_psm_tab(ofstream &osTarget, ofstream &osDecoy);
  void get_pep_seq(string &pep, string &seq, string &n, string &c);
  void write_results_psm_xml(PepXMLWriter& osTarget, PepXMLWriter& osDecoy);
  void computePEP();

  void write_max_nets(string filename, NeuralNet *max_net);
  void write_unique_peptides(string filename, NeuralNet* max_net);
  void write_num_psm_per_spectrum(NeuralNet* max_net);

  inline void set_fileroot(string &fl){fileroot = fl;}
  inline void set_overwrite_flag(int flag) {overwrite_flag = flag;}
  inline void set_input_dir(string &input_dir) {in_dir = input_dir; d.set_input_dir(input_dir);}
  inline void set_output_dir(string &output_dir){out_dir = output_dir;}
  void print_description();
  int set_command_line_options(int argc, char **argv);
  int crux_set_command_line_options(int argc, char *argv[]);
  

  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  virtual bool needsOutputDirectory();
  virtual COMMAND_T getCommand(); 
  FILE_FORMAT_T check_file_format(string &filePath);
  string file_extension(string filename); 
  void get_protein_id(int pepind, vector<string> &prot); 
  void  print_protein_ids(vector<string> &prots, ofstream &os, int psmind);  

protected:

    Dataset d;
    string res_prefix;

    PSMScores fullset; 
    PSMScores trainset,testset,thresholdset;

    int seed;
    double selectionfdr;
    
    int num_features;
    NeuralNet net;
    int num_hu;
    double mu;
    double weightDecay;
    string cleavage_type; 
    int ind_low;
    int interval;
    int niter;
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
    NeuralNet* nets;

    string in_dir;
    string out_dir;
    int skip_cleanup_flag;
    int overwrite_flag;
    string fileroot;
    int feature_file_flag;
    ostringstream feature_file_name;
    
    string file_format_; 
    
    TabDelimParser pars;
    SQTParser* parser; 
     

};



#endif /*QRANKER_H_*/
