#ifndef PROTSCORES_H_
#define PROTSCORES_H_
#include <vector>
#include <algorithm>
#include <string>
using namespace std;
#include "DataSet.h"

class ProtScoreHolder{
public:
  double score; 
  int protind;
  int label;
  double q;
  double nsaf;
  double PEP;
  //meta group number
  int group_number;
  //flag indicating that it is protein identical to some other protein in terms of peptide composition
  int indistinguishable_prot;
  int subset_prot;
  //inds of identical prots in terms of their peptide composition
  vector<int> indistinguishable_protinds;
  vector<int> subset_protinds;
  vector<int> parent_groups;
  
  int has_complement;
  //intersection of peptides in a group
  vector<int> intersection;
  map<int, vector<int> > protind2complement;
 ProtScoreHolder():score(0.0),protind(0),q(0.0),nsaf(0.0),PEP(0.0),group_number(0),indistinguishable_prot(0),subset_prot(0), has_complement(0){;}
  ~ProtScoreHolder() {;}
};


class ProtScores
{
public:
	ProtScores();
	~ProtScores();
	void clear();
    vector<ProtScoreHolder>::iterator begin() {return scores.begin();}
    vector<ProtScoreHolder>::iterator end() {return scores.end();}    
    static double pi0;
    double factor;

    int calcOverFDR(double fdr);
    void calcMultiOverFDR(vector<double> &fdr, vector<int> &overFDR);
    inline ProtScoreHolder& operator[](int ix){return scores[ix];}
    void static fillProteinsFull(ProtScores& fullset, Dataset &d);
    int static traverse_peptide(Dataset &d, int pepind, int trn_tst, vector<int> &assignment_array);
    int static traverse_protein(Dataset &d, int protind, int trn_tst, vector<int> &assignment_array);
    void static fillProteinsSplit(ProtScores& train,ProtScores& test,Dataset &d, double ratio);
    void make_meta_set(Dataset &d);
    bool is_subset(Dataset &d, int protind1, int protind2);
    bool are_equal(Dataset &d, int protind1, int protind2, vector<int> &overlap, vector<int> &only1, vector<int> &only2);
    void get_pep_seq(string &pep, string &seq);
    int is_equivalent_pep(Dataset &d, int protind1, string &pep2);
    bool is_equivalent(Dataset &d, int protind1, vector<int> &pep_set);
    void get_intersection(vector<int> &set1, vector<int> &set2, vector<int> &result);
    void get_complement(Dataset &d, int protind1, vector<int> &set2, vector<int> &result);

    inline int get_num_subsets(){return scores_subset.size();}
    inline ProtScoreHolder& get_subset_prot(int i){return scores_subset[i];}

    inline int size(){return scores.size();}
protected:
    int neg,pos,posNow;
    vector<ProtScoreHolder> scores;
    vector<ProtScoreHolder> scores_meta;
    vector<ProtScoreHolder> scores_same;
    vector<ProtScoreHolder> scores_subset;

};



#endif /*PROTSCORES_H_*/
