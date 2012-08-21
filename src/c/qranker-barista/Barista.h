#ifndef BARISTA_H
#define BARISTA_H
#define CRUX
#include <sys/types.h>
#ifndef _MSC_VER
#include <dirent.h>
#endif
#include <iostream>
#include <sstream>
#include <set>
#include <algorithm>
#include <assert.h>
#include <cstdio>
#include <iomanip>
#ifdef CRUX
#include "CruxApplication.h"
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#endif
#include "analyze_psms.h"
#include "PepXMLWriter.h"
#include "DataSet.h"
#include "ProtScores.h"
#include "PSMScores.h"
#include "PepScores.h"
#include "NeuralNet.h"
#include "SQTParser.h"
#include "CruxParser.h"
#include "QRanker.h"
#include "PepRanker.h"
using namespace std;
#include "mass.h"

class Barista : public CruxApplication
{
 public:
  Barista() 
    : verbose(0),
    skip_cleanup_flag(0),
    overwrite_flag(0),
    feature_file_flag(0),
    list_of_files_flag(0),
    in_dir(""), 
    out_dir(""), 
    fileroot(""), 
    seed(0), 
    selectionfdr(0.01), 
    nepochs(20), 
    num_features(0), 
    num_hu(3), 
    mu(0.05),
    weightDecay(0.0), 
    alpha(0.3),
    max_psms_in_prot(0),
    net_clones(0),
    max_fdr(0),
    max_peptides(0),   
    max_fdr_psm(0),
    max_fdr_pep(0){}
  ~Barista(){clear();}
  void clear();
  void print_description();
  int crux_set_command_line_options(int argc, char *argv[]);
  int set_command_line_options(int argc, char *argv[]);
  void setup_for_training(int trn_to_tst);
  void setup_for_reporting_results();
  int run();

  int run_tries();
  int run_tries_multi_task();
  double train_hinge(int protind, int label);
  double train_hinge_psm(int psmind, int label);
  void train_net(double selectionfdr);
  void train_net_multi_task(double selectionfdr, int interval);

  void calc_gradients(int protind, int label);

  int getOverFDRProt(ProtScores &set, NeuralNet &n, double fdr);
  int getOverFDRProt(ProtScores &set, double fdr);

  double get_protein_score(int protind);
  double get_protein_score(int protind, NeuralNet &n);
  double get_protein_score_parsimonious(int protind, NeuralNet &n);
  int getOverFDRProtParsimonious(ProtScores &set, NeuralNet &n, double fdr);
  void computePEP();
  int computeNSAF();
  int computePepNSAF();

  void write_results_prot(string &out_dir, int fdr);
  void report_all_results();
  void get_pep_seq(string &pep, string &seq, string &n, string &c);
  void get_tab_delim_proteins(string protein_str, vector<string> &proteins);
  void write_protein_special_case_xml(ofstream &os, int i);
  void write_results_prot_xml(ofstream &os);
  void write_subset_protein_special_case_xml(ofstream &os, int i);
  void write_results_subset_prot_xml(ofstream &os);
  void write_results_peptides_xml(ofstream &os);
  void write_results_psm_xml(ofstream &os);
  void write_results_pep_xml(PepXMLWriter& xmlfile);
  void report_all_results_xml();
  void write_results_prot_special_case_tab(ofstream &os, int i);
  void write_results_prot_tab(ofstream &os);
  void write_subset_protein_special_case_tab(ofstream &os, int i);
  void write_results_subset_prot_tab(ofstream &os);
  void write_results_peptides_tab(ofstream &os);
  void write_results_psm_tab(ofstream &os);
  void report_all_results_tab();
  void report_all_results_xml_tab();
  void report_prot_fdr_counts(vector<double> &qvals, ofstream &of);
  void report_psm_fdr_counts(vector<double> &qvals, ofstream &of);
  void report_pep_fdr_counts(vector<double> &qvals, ofstream &of);
  void report_all_fdr_counts();
  void clean_up();
  void get_protein_id(int pepind, vector<string> &prots);
  void print_protein_ids(vector<string> &proteins,ofstream &os,int psmind);  

  int getOverFDRPSM(PSMScores &set, NeuralNet &n, double fdr);
  double get_peptide_score(int pepind, NeuralNet &n);
  int getOverFDRPep(PepScores &set, NeuralNet &n, double fdr);

  inline void set_input_dir(string input_dir) {in_dir = input_dir; d.set_input_dir(input_dir);}
  inline void set_output_dir(string output_dir){out_dir = output_dir;}

  /* CruxApplication Methods */
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  virtual bool needsOutputDirectory();
  virtual COMMAND_T getCommand();


  double check_gradients_hinge_one_net(int protind, int label);
  double check_gradients_hinge_clones(int protind, int label);

  FILE_FORMAT_T check_file_format(string filename);
  string file_extension(string str); 
 protected:
  SQTParser* parser; 
  int verbose;
  int skip_cleanup_flag;
  int overwrite_flag;
  int feature_file_flag;
  ostringstream feature_file_name;
  int list_of_files_flag; 

  Dataset d;
  string in_dir;
  string out_dir;
  string fileroot;

  string cleavage_type; 

  int seed;
  double selectionfdr;
  int nepochs;

  NeuralNet net;
  int num_features;
  int num_hu;
  double mu;
  double weightDecay;
  double alpha;


  int max_psms_in_prot;
  NeuralNet *net_clones;
  NeuralNet max_net_prot;
  int max_fdr;
  ProtScores trainset, thresholdset, testset;

  int max_peptides;
  vector<int> max_psm_inds;
  vector<double> max_psm_scores;
  
  //for parsimony counts
  vector<int> used_peptides;
  vector<int> pepind_to_max_psmind;

  PSMScores psmtrainset, psmtestset;
  NeuralNet max_net_psm;
  int max_fdr_psm;
  
  PepScores peptrainset,peptestset;
  NeuralNet max_net_pep;
  int max_fdr_pep;
  
  string file_format_; 
  ofstream fdebug;

  string opt_type;
  QRanker qr;
  PepRanker pr;

};


#endif
