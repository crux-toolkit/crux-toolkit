#ifndef PEPRANKER_H_
#define PEPRANKER_H_
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
#include "PepScores.h"
#include "SQTParser.h"
#include "CruxParser.h"
#include "TabDelimParser.h"
#include "CruxApplication.h"
#include "mass.h"
#include "objects.h"

class PepRanker: public CruxApplication
{

public:
  PepRanker();
  virtual ~PepRanker();

  int run();
  int getOverFDRPSM(PSMScores &set, NeuralNet &n, double fdr);
  double get_peptide_score_xcorr(int pepind);
  void getMultiFDRXCorr(PepScores &set, vector<double> &qvalues);
  double get_peptide_score(int pepind, NeuralNet &n);
  int getOverFDR(PepScores &set, NeuralNet &n, double fdr);
  void getMultiFDR(PepScores &set, NeuralNet &n, vector<double> &qvalues);
  void write_max_nets(string filename, NeuralNet *max_net);
  double get_protein_score(int pepind);
  void calc_gradients(int pepind, int label, int i);
  void train_net_ranking(PepScores &set, int interval);
  void train_many_general_nets();
  void train_many_target_nets();
  void train_many_nets();
  void setup_for_training();  
  void setup_for_reporting_results();
  void report_results_xml_tab();
    
  void get_tab_delim_proteins(string protein_str, vector<string> &proteins);
  void get_pep_seq(string &pep, string &seq, string &n, string &c);
  void write_results_psm_xml(PepXMLWriter& os);
  void computePEP();
  int computePepNSAF();

  void write_results_peptides_tab(ofstream &os);
  void write_results_psm_tab(ofstream &os);
  void report_results_tab();
  void write_results_pep_xml(PepXMLWriter& xmlfile);
  void write_results_peptides_xml(ofstream &os);
  void write_results_psm_xml(ofstream &os);
  void report_results_xml();

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

    PepScores fullset; 
    PepScores trainset,testset,thresholdset;

    PSMScores psmfullset;
    PSMScores psmtrainset, psmtestset;

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
    vector<int> overFDRmulti;
    vector<int> max_overFDR;

    NeuralNet* max_net_gen;
    NeuralNet* max_net_targ;
    NeuralNet *net_clones;
    
    int psm_count;
    vector<int> max_psm_inds;
    vector<int> pepind_to_max_psmind;

    
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
