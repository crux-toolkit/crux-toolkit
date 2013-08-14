#ifndef SQTPARSER_H
#define SQTPARSER_H
#define CRUX
#include <sys/stat.h>
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
#include "SpecFeatures.h"
#include "BipartiteGraph.h"
#ifdef CRUX
#include "CruxApplication.h"
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#endif
using namespace std;

typedef enum {TRYPSIN_ENZ,CHYMOTRYPSIN_ENZ,ELASTASE_ENZ} enzyme;

class SQTParser{
 public:
  struct sqt_match{
    int scan;
    int charge;
    double precursor_mass;
    int num_sequence_comparisons;
    vector <int> xcorr_rank;
    vector<int> sp_rank;
    vector<double> calc_mass;
    vector<double> delta_cn;
    vector<double> xcorr_score;
    vector<double> sp_score;
    vector <int> num_ions_matched;
    vector<double> num_total_ions;
    vector<string> peptides;
    vector<int> num_proteins_in_match;
    vector<string> proteins;
    vector<int> peptide_pos; 

  };
  /*
 SQTParser() : num_spectra(0),num_psm(0),num_pos_psm(0),num_neg_psm(0),num_features(0),
    num_spec_features(0), num_pep(0),num_pos_pep(0), num_neg_pep(0),num_prot(0),num_pos_prot(0),num_neg_prot(0),
   num_mixed_labels(0),psmind(0){}
  */
  SQTParser();
  virtual ~SQTParser();
  void clear();
  virtual int run();
  void clear_matches();
  void erase_matches();
  inline void set_num_features(int nf) {num_features = nf;}
  inline void set_num_spec_features(int nsf) {num_spec_features = nsf;}
  inline int get_num_features() const {return num_features;}
  inline void set_use_quadratic_features(int use){use_quadratic_features = use;}
  inline int get_use_quadratic_features()const {return use_quadratic_features;}
  int set_output_dir(string &output_dir, int overwrite_flag);
  int set_output_dir(string &output_dir){out_dir = output_dir; return 1;}
  inline string& get_output_dir(){return out_dir;}
  inline void set_input_dir(string input_dir){in_dir = input_dir;}
  inline void set_input_dir_ms2(string input_dir){in_dir_ms2 = input_dir;}
  inline string& get_input_dir(){return in_dir;}
  inline void set_db_name(string database){db_name = database;}
  inline void set_decoy_prefix(string prefix){decoy_prefix = prefix;}
  inline void set_num_hits_per_spectrum(int hits_per_spectrum){fhps = hits_per_spectrum;}
  inline vector<string> get_final_features_header(){return final_features_header_;}
  void set_enzyme(string &enz);
  int is_ending(string &name, const string &ext);
  int is_spectrum_file(string &fname);
  int is_fasta(string &fname);
  int set_database_source(string &db_source);
  virtual int match_file_to_ms2(string &sqt_source, string &prefix);
  int collect_ms2_files(string &ms2_source, string &sqt_source);
  virtual int set_input_sources(string &ms2_source, string &sqt_source);
  void read_list_of_files(string &list, vector<string> &fnames);
  
  /********* for separate database searches ********************************************/
  virtual int match_target_file_to_ms2(string &sqt_source, string &prefix);
  virtual int match_decoy_file_to_ms2(string &sqt_source, string &prefix);
  int collect_ms2_files(string &ms2_source, string &sqt_target_source, string &sqt_decoy_source);
  virtual int set_input_sources(string &ms2_source, string &sqt_target_source, string &sqt_decoy_source);
  /*************************************************************************************/


  void read_sqt_file(ifstream &is, string &decoy_prefix, int final_hits_per_spectrum, enzyme enz, bool decoy);
  int parse_sqt_spectrum_matches(ifstream &is, sqt_match &m);
  void read_S_line(ifstream &is, sqt_match &m);
  void read_M_line(ifstream &is, sqt_match &m);
  void digest_database(ifstream &f_db, enzyme e);
  int cntEnzConstraints(string& seq,enzyme enz);
  
  static int cntEnz(const string& peptide, enzyme enz);
  static double isTryptic(const char n,const char c);
  static double isChymoTryptic(const char n,const char c);
  static double isElastasic(const char n,const char c);
  static double isEnz(const char n,const char c, enzyme enz);
  void extract_psm_features(sqt_match &m, enzyme enz, double *x, int i);
  void extract_psm_features(sqt_match &m, enzyme enz, double *x, int i, int hits_read);
  void extract_features(sqt_match &m, int hits_read, int final_hits,enzyme enz);
  void add_matches_to_tables(sqt_match &m, string &decoy_prefix, int hits_read, int final_hits, bool decoy);
  void allocate_feature_space();
  void fill_graphs_and_save_data(string &out_dir);
  
  void open_files(string &out_dir);
  void close_files();
  void clean_up(string dir);
  int check_file(ostringstream &fname);
  int check_input_dir(string &in_dir);
  
  //spec features generator
  SpecFeaturesGenerator sfg;
  
  void write_features_header();
  void add_quadratic_features_header(); 
  //void write_quad_features_header();
  
 protected:
  //auxiliary variables
  int database_exists;
  int num_prot_not_found_in_db;
  int num_pos_prot_not_found_in_db;
  int num_neg_prot_not_found_in_db;
   
  vector<string> features_header_; 
  vector<string> spec_features_header_3_;
  vector<string> spec_features_header_7_;
  vector<string> final_features_header_; 
  sqt_match m;
  int num_mixed_labels;
  map<int,set<int> > pepind_to_protinds_map;
  map<int,set<int> > protind_to_pepinds_map;
  map<int,set<int> > pepind_to_psminds_map;
  
  //summary of the dataset
  int num_features;
  int num_spec_features;
  int num_total_features;
  int use_quadratic_features;
  int num_spectra;
  int num_psm;
  int num_pos_psm;
  int num_neg_psm;
  int num_pep;
  int num_pos_pep;
  int num_neg_pep;
  int num_prot;
  int num_pos_prot;
  int num_neg_prot;
  
  int num_cur_psm;
  int num_cur_prot;
  int prot_offset;
  
  

  //psm feature vector
  double *x;
  //spec feature vector
  double *xs;
  
  //peptide info
  map<string,int> pep_to_ind;
  map<int,string> ind_to_pep;
  BipartiteGraph pepind_to_protinds;
  BipartiteGraph pepind_to_psminds;
  
  //protein info
  map<string,int> prot_to_ind;
  map<int,string> ind_to_prot;
  BipartiteGraph protind_to_pepinds;
  
  //digested database info
  map<string,int> protein_to_num_all_pep_map;  
  map<int,int> protind_to_num_all_pep_map;
  map<string,int> protein_to_length_map;  
  map<int,int> protind_to_length_map;
  int *protind_to_num_all_pep;
  int *protind_to_length;

  //writing out data
  string in_dir;
  string in_dir_ms2;
  string out_dir;
  string db_name;
  vector<string> sqt_file_names;
  vector<string> ms2_file_names;
  vector<string> db_file_names;

  string cur_fname;
  int cur_fileind;
  
  //files for writing out data
  ofstream f_psm;
  ofstream f_psmind_to_label;
  ofstream f_psmind_to_scan;
  ofstream f_psmind_to_charge;
  ofstream f_psmind_to_precursor_mass;
  ofstream f_psmind_to_pepind;
  ofstream f_pepind_to_label;
  ofstream f_protind_to_label;
  ofstream f_protind_to_num_all_pep;
  ofstream f_protind_to_length;
  ofstream f_fileind_to_fname;
  ofstream f_psmind_to_fileind;
  
  ofstream f_psmind_to_xcorr;
  ofstream f_psmind_to_spscore;
  ofstream f_psmind_to_deltaCn;
  ofstream f_psmind_to_calculated_mass;
  
  ofstream f_psmind_to_sp_rank;//sp rank
  ofstream f_pmsind_to_matches_spectrum; //matches_spectrum  
  ofstream f_psmind_to_xcorr_rank;//xcorr rank 
  ofstream f_psmind_to_by_ions_matched;// b/y ions match  
  ofstream f_psmind_to_by_ions_total;  //b/y ions total   
  ofstream f_psmind_to_peptide_position; //peptide position 
  
  //final hits per spectrum
  int fhps;
  //decoy prefix
  string decoy_prefix;
  //enzyme
  enzyme e;
  //max peptide length to be considered
  int max_len;
  //min peptide length to be considered
  int min_len;

  
  virtual bool read_search_results(string& cur_fname, bool decoy); 
  virtual string get_parser_extension();

  vector<string> spectrumExts_; // spectrum file extensions
  
};

#endif
