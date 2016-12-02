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

#include "app/CruxApplication.h"
#include "io/carp.h"
#include "util/crux-utils.h"
#include "parameter.h"

#include "io/PepXMLWriter.h"
#include "NeuralNet.h"
#include "DataSetCrux.h"
#include "PSMScores.h"
#include "SQTParser.h"
#include "CruxParser.h"
#include "TabDelimParser.h"
#include "app/CruxApplication.h"
#include "util/mass.h"
#include "model/objects.h"

double* compute_PEP(double* target_scores, ///< scores for target matches
  int num_targets,       ///< size of target_scores
  double* decoy_scores,  ///< scores for decoy matches
  int num_decoys,         ///< size of decoy_scores
  bool ascending ///< are the scores ascending or descending
  );

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
  virtual std::string getName() const;
  virtual std::string getDescription() const;
  virtual std::vector<std::string> getArgs() const;
  virtual std::vector<std::string> getOptions() const;
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;
  virtual bool needsOutputDirectory() const;
  virtual COMMAND_T getCommand() const;
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
