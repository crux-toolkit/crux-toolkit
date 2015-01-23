#ifndef TABDELIMPARSER_H
#define TABDELIMPARSER_H

#include <sys/types.h>
#ifndef _MSC_VER
#include <dirent.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <string>
#include <math.h>
#include <map>
#include <set>
#include <cstring>
#include <stdlib.h>
#include "SpecFeatures.h"
#include "BipartiteGraph.h"
using namespace std;

class TabDelimParser{
 public:

  TabDelimParser();
  ~TabDelimParser();
  void clear();
  inline void set_output_dir(string &dir){out_dir = dir;}
  inline string& get_output_dir(){return out_dir;}
  int run(vector<string> &filenames);

  void get_tokens(string &line, vector<string>&tokens, string &delim);
  void first_pass(ifstream &fin);
  void second_pass(ifstream &fin, int label);
  void allocate_feature_space();
  void extract_psm_features(vector<string> & tokens, double *x);
  void save_data_in_binary(string out_dir);
  void clean_up(string dir);
  /*******************************************************/
  int run_on_xlink(vector<string> &filenames);
  void first_pass_xlink(ifstream &fin);
  void second_pass_xlink(ifstream &fin, int label);
  void allocate_feature_space_xlink();
  void extract_xlink_features(vector<string> & tokens, double *x);
  void save_data_in_binary_xlink(string out_dir);
  void clean_up_xlink(string dir);
 protected:
  //auxiliary variables
  int num_mixed_labels;
  
  //psm info
  int* psmind_to_scan;
  int* psmind_to_charge;
  int* psmind_to_label;
  int* psmind_to_num_pep;
  int* psmind_to_ofst;
  int* psmind_to_pepind;
  double *psmind_to_neutral_mass;
  double *psmind_to_peptide_mass;

  //peptide info
  map<string,int> pep_to_ind;
  map<int,string> ind_to_pep;
  
  //summary of the dataset
  int num_features;
  int num_psm;
  int num_pos_psm;
  int num_neg_psm;
  int num_pep;
  int num_pep_in_all_psms;
  int curr_ofst;

  int psmind;

  //psm feature vector
  double *x;

  //writing out data
  string out_dir;
  ofstream f_psm;

  //final hits per spectrum
  int fhps;
  //decoy prefix
  string decoy_prefix;
  //enzyme
  // enzyme e;
  //max peptide length to be considered
  int max_len;
  //min peptide length to be considered
  int min_len;
  
  /************xlink-specific*************/
  int num_xlink_features;
  map<int,string> psmind_to_peptide1;
  map<int,string> psmind_to_peptide2;
  map<int,string> psmind_to_loc;
  map<int,string> psmind_to_protein1;
  map<int,string> psmind_to_protein2;
  
};

#endif
