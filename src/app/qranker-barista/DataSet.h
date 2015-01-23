#ifndef DATASET_H_
#define DATASET_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <cmath>
#include <map>
#include "BipartiteGraph.h"
using namespace std;



class Dataset
{
 public:

  Dataset();
  ~Dataset();
  
  void load_data_psm_training();
  void clear_data_psm_training();
  void load_labels_psm_training();
  void clear_labels_psm_training();
  void load_data_psm_results();
  void clear_data_psm_results();

  void load_data_prot_training();
  void clear_data_prot_training();
  void load_labels_prot_training();
  void clear_labels_prot_training();
  void load_data_all_results();
  void clear_data_all_results();
  void load_data_pep_training();
  void clear_data_pep_training();
  void load_labels_pep_training();
  void clear_labels_pep_training();
  void load_data_pep_results();
  void clear_data_pep_results();

  inline void set_input_dir(string input_dir){in_dir = input_dir;}
  void normalize_psms();
  int print_features(string &filename);
  
  inline double* psmind2features(int psmind){return (psmind_to_features+num_features*psmind);}
  inline int psmind2label(int psmind){return psmind_to_label[psmind];}
  inline int psmind2scan(int psmind){return psmind_to_scan[psmind];}
  inline int psmind2charge(int psmind){return psmind_to_charge[psmind];}
  inline int psmind2pepind(int psmind){return psmind_to_pepind[psmind];}
  inline double psmind2precursor_mass(int psmind){return psmind_to_precursor_mass[psmind];}
  inline double psmind2xcorr(int psmind){return psmind_to_xcorr[psmind];}
  inline double psmind2spscore(int psmind){return psmind_to_spscore[psmind];}
  inline double psmind2deltaCn(int psmind){return psmind_to_deltaCn[psmind];}
  inline double psmind2peptide_mass(int psmind){return psmind_to_calculated_mass[psmind];}
  string& psmind2fname(int psmind);
  inline string& ind2pep(int ind){return ind_to_pep[ind];}
  inline int get_num_psms(){return num_psms;}
  inline int get_num_features(){return num_features;}

  inline int get_num_peptides(){return num_pep;}
  inline int pepind2label(int pepind){return pepind_to_label[pepind];}
  inline int pepind2num_psm(int pepind){return pepind_to_psminds.get_range_length(pepind);}
  inline int* pepind2psminds(int pepind){return pepind_to_psminds.get_range_indices(pepind);}
  inline int pepind2num_prot(int pepind){return pepind_to_protinds.get_range_length(pepind);}
  inline int* pepind2protinds(int pepind){return pepind_to_protinds.get_range_indices(pepind);}

  inline int get_num_proteins(){return num_prot;}
  inline string& ind2prot(int ind){return ind_to_prot[ind];}
  inline int protind2label(int protind){return protind_to_label[protind];}
  inline int protind2num_pep(int protind){return protind_to_pepinds.get_range_length(protind);}
  inline int* protind2pepinds(int protind){return protind_to_pepinds.get_range_indices(protind);}
  inline int protind2num_all_pep(int protind){return protind_to_num_all_pep[protind];}
  inline int protind2length(int protind){return protind_to_length[protind];}
  //returns false if not subset, true if yes, subset
  inline bool is_prot_subset(int protind1, int protind2){return protind_to_pepinds.is_subset(protind1, protind2);} 
  inline int psmind2sp_rank(int psmind){return psmind_to_sp_rank[psmind];}//Sp rank 
  inline int psmind2xcorr_rank(int psmind){return psmind_to_xcorr_rank[psmind];}//xcorr rank
  inline double psmind2by_ions_matched(int psmind){return psmind_to_by_ions_matched[psmind];}//b/y ions matched 
  inline double psmind2by_ions_total(int psmind){return psmind_to_by_ions_total[psmind];}//b/y ions total 
  inline int psmind2matches_spectrum(int psmind) {return psmind_to_matches_spectrum[psmind];}///<matchs/spectrum
  inline int psmind2peptide_position(int psmind){return psmind_to_peptide_position[psmind];}///<peptide position in protein 
  
  inline void get_features_header(vector<string> str){features_header_.swap(str);}


 protected:
  int num_psms;
  int num_pos_psms;
  int num_neg_psms;
  int num_features;
  double* psmind_to_features;
  int* psmind_to_label;
  int *psmind_to_pepind;
  int *psmind_to_scan;
  int *psmind_to_charge;
  double *psmind_to_xcorr;
  double *psmind_to_spscore;
  double *psmind_to_deltaCn;
  double *psmind_to_calculated_mass;
  double *psmind_to_precursor_mass;
  int *psmind_to_fileind;
  int num_pep;
  int num_pos_pep;
  int num_neg_pep;
  int *pepind_to_label;
  int num_prot;
  int num_pos_prot;
  int num_neg_prot;
  int *protind_to_label;
  int *protind_to_num_all_pep;
  int* psmind_to_sp_rank;//Sp rank 
  int *protind_to_length;
  int* psmind_to_xcorr_rank;//xcorr rank 
  int *psmind_to_matches_spectrum;//matches/spectrum
  double* psmind_to_by_ions_matched; //b/y ions matched
  double* psmind_to_by_ions_total; //b/y ions total
  int*  psmind_to_peptide_position;//<pepetide position in protein 

  map <int, string> fileind_to_fname;
 
  vector <string> features_header_;
  
  BipartiteGraph pepind_to_psminds;
  BipartiteGraph pepind_to_protinds;
  map <int, string> ind_to_pep;
  

 
  
  BipartiteGraph protind_to_pepinds;
 
  map <int, string> ind_to_prot;

  string in_dir;
};



#endif /*DATASET_H */
