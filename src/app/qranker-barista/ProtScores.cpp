#include "ProtScores.h"


//inline bool operator<(const ProtScoreHolder &one, const ProtScoreHolder &other) 
//{return (one.score<other.score);}

inline bool operator<(const ProtScoreHolder &one, const ProtScoreHolder &other) 
{return (one.score<other.score) || (one.score == other.score && one.protind < other.protind);}



ProtScores::ProtScores()
{
    factor=1;
    neg=0;
    pos=0;
    posNow=0;
}

ProtScores::~ProtScores()
{
}


void ProtScores::clear()
{
  scores.clear();
  pos=neg=posNow = 0;
}

/**
 * Percentage of target scores that are drawn according to the null.
 */
double ProtScores::pi0 = 0.9;

void ProtScores :: fillProteinsFull(ProtScores& fullset, Dataset &d)
{
  int n = d.get_num_proteins();
  ProtScoreHolder s;
  fullset.scores.resize(n,s);
  int num_pos = 0;
  int num_neg = 0;
  for(int i = 0; i < n; i++)
    {
      fullset[i].protind = i;
      int label = d.protind2label(i);
      fullset[i].label = label;
      if(label == 1)
	num_pos++;
      else
	num_neg++;
    }

  fullset.pos = num_pos;
  fullset.neg = num_neg;
  fullset.factor = (double)num_pos/(double)num_neg;
}


int ProtScores :: traverse_peptide(Dataset &d, int pepind, int trn_tst, vector<int> &assignment_array)
{
  int num_added = 0;
  int num_prot = d.pepind2num_prot(pepind);
  int *protinds = d.pepind2protinds(pepind);
  for (int i = 0; i < num_prot;i++)
    {
      int protind = protinds[i];
      if(assignment_array[protind] == 0)
	{
	  assignment_array[protind] = trn_tst;
	  int num = traverse_protein(d,protind,trn_tst,assignment_array);
	  num_added+=num;
	}
      else
	{
	  if(assignment_array[protind] != trn_tst)
	    cout << "warning: conflict spliting into trn-tst " << assignment_array[protind] << endl;
	}
    }
  return num_added;
}

int ProtScores :: traverse_protein(Dataset &d, int protind, int trn_tst, vector<int> &assignment_array)
{
  int num_pep = d.protind2num_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  int num_added = 0;
  for(int i = 0; i < num_pep; i++)
    {
      int pepind = pepinds[i];
      int num = traverse_peptide(d,pepind,trn_tst, assignment_array);
      num_added += num;
    }
  return num_added+1;
}

void ProtScores :: fillProteinsSplit(ProtScores& train,ProtScores& test,Dataset &d, double ratio)
{
  int n = d.get_num_proteins();
  vector<int> assignment_array;
  assignment_array.resize(n,0);
  int num_trn = 0;
  int num_tst = 0;
  int trn_tst = 0;
  int max_num = 0;
  for(int i = 0; i < n; i++)
    {
      if(assignment_array[i] == 0)
	{
	  if (num_trn >= ratio*num_tst) 
	    trn_tst = 2;
	 else
	    trn_tst = 1;
	  assignment_array[i] = trn_tst;
	  int num = traverse_protein(d, i, trn_tst, assignment_array);
	  if(num > max_num)
	    {
	      max_num = num;
	      //cout << "max_num " << max_num << endl;
	    }
	
	  if (trn_tst == 1) 
	    num_trn += num;
	  else
	    num_tst += num;
	} 
    } 
  //cout << "total proteins " << num_trn+num_tst << " trainset size " << num_trn << " testeset size " << num_tst << " ratio " << (double)num_trn/(double)num_tst << endl;
  //cout << "total proteins " << num_trn+num_tst << " trainset size " << num_trn << " testeset size " << num_tst << endl;

  int trn= 0;
  int tst= 0;
  int num_pos_trn = 0;
  int num_neg_trn = 0;
  int num_pos_tst = 0;
  int num_neg_tst = 0;

  ProtScoreHolder s;
  train.scores.resize(num_trn,s);
  test.scores.resize(num_tst,s);

  for(unsigned int i = 0; i < assignment_array.size();i++)
    {
      if(assignment_array[i] == 1)
	{
	  train[trn].protind = i;
	  int label = d.protind2label(i);
	  train[trn].label = label;
	  if(label == 1)
	    num_pos_trn++;
	  else
	    num_neg_trn++;
	  trn++;
	}
      else if (assignment_array[i] == 2)
	{
	  test[tst].protind = i;
	  int label = d.protind2label(i);
	  test[tst].label = label;
	  if(label == 1)
	    num_pos_tst++;
	  else
	    num_neg_tst++;
	  tst++;
	}
      else
	cout << "protein_not assigned!\n";
    }

    train.pos=num_pos_trn;
    test.pos=num_pos_tst;
    train.neg=num_neg_trn;
    test.neg=num_neg_tst;
    train.factor = (double)train.pos/(double)train.neg;
    test.factor = (double)train.pos/(double)train.neg;

}


int ProtScores::calcOverFDR(double fdr) {
  
  vector<ProtScoreHolder>::iterator it = scores.begin();

  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());
  
  /*
  ProtScoreHolder s;
  for(int i = scores.size()-1; i > 0; i--)
    {
      if(scores[i].score == scores[i-1].score)
	{
	  if(scores[i].protind < scores[i-1].protind)
	    {
	      s = scores[i];
	      scores[i] = scores[i-1];
	      scores[i-1] = scores[i];
	    }
	}
    }
  */

  int positives=0,nulls=0;
  double efp=0.0,q;
  posNow = 0;
  register unsigned int ix=0;
  for(it=scores.begin();it!=scores.end();it++) {
    if (it->label!=-1)
      positives++;
    if (it->label==-1) {
      nulls++;
      efp=pi0*nulls*factor;
    }
    if (positives)
      q=efp/(double)positives;
    else
      q=pi0;
    if (q>pi0)
      q=pi0;
    it->q=q;
    if (fdr>=q)
      posNow = positives;
  }

  for (ix=scores.size();--ix;) {
    if (scores[ix-1].q > scores[ix].q)
      scores[ix-1].q = scores[ix].q;  
  }
  return posNow;
}



void ProtScores :: get_pep_seq(string &pep, string &seq)
{
  string tmp;
  size_t pos = pep.find('.');
  if(pos > 0 && pos != string::npos)
    tmp = pep.substr(pos+1,pep.size());
  pos = tmp.rfind('.');
  if(pos > 0 && pos != string::npos)
    tmp = tmp.substr(0,pos);
  seq = tmp;
}

int ProtScores :: is_equivalent_pep(Dataset &d, int protind1, string &pep2)
{

  int num_pep1 = d.protind2num_pep(protind1);
  int *pepinds1 = d.protind2pepinds(protind1);

  string p2;
  get_pep_seq(pep2,p2);
  
  for(int j = 0; j < num_pep1; j++)
    {
      string pep1 = d.ind2pep(pepinds1[j]);
      
      int flag = 1;
      if(pep1.size() != pep2.size())
	flag = 0;
      
      string p1;
      get_pep_seq(pep1, p1);
      
      unsigned int ix1 = 0; unsigned int ix2 = 0; 
      while((flag == 1) && (ix1 < p1.size()) && (ix2 < p2.size()))
	{
	  char c1 = p1.at(ix1);
	  char c2 = p2.at(ix2);
	  
	  if(c1 != c2)
	    {
	      flag = 0;
	      
	      if(c1 == 'L' && c2 =='I')
		flag = 1;
	      
	      if(c1 == 'I' && c2 =='L')
		flag = 1;
	      if(c1 == 'T' && c2 =='S')
		flag = 1;
	      if(c1 == 'S' && c2 =='T')
		flag = 1;
	    }
	  ix1++; ix2++;
	}
      if(flag == 1)
	return pepinds1[j];
	
    }
  return 0;
    
}

bool ProtScores :: is_equivalent(Dataset &d, int protind1, vector<int> &pep_set)
{
  
  for(unsigned int i = 0; i < pep_set.size(); i++)
    {
      string pep2 = d.ind2pep(pep_set[i]);
      if(is_equivalent_pep(d, protind1, pep2) == 0)
	return false;
    }
  return true;
}



bool ProtScores :: is_subset(Dataset &d, int protind1, int protind2)
{
  int num_pep1 = d.protind2num_pep(protind1);
  int *pepinds1 = d.protind2pepinds(protind1);
  int num_pep2 = d.protind2num_pep(protind2);
  int *pepinds2 = d.protind2pepinds(protind2);

  vector<int> only2;

  int ix1 = 0;
  int ix2 = 0;
  
  while(ix1 < num_pep1 && ix2 < num_pep2)
    {
      if(pepinds1[ix1] == pepinds2[ix2])
	{
	  ix1++; ix2++;
	}
      else if(pepinds1[ix1] < pepinds2[ix2])
	ix1++;
      else
	{
	  only2.push_back(pepinds2[ix2]); ix2++;
	}
    }
  while(ix2 < num_pep2)
    {
      only2.push_back(pepinds2[ix2]);
      ix2++;
    }
  if(only2.size() == 0)
    return true;
  else 
    return is_equivalent(d, protind1, only2);
}

bool ProtScores :: are_equal(Dataset &d, int protind1, int protind2, vector<int> &overlap, vector<int> &only1, vector<int> &only2)
{
  int num_pep1 = d.protind2num_pep(protind1);
  int *pepinds1 = d.protind2pepinds(protind1);
  int num_pep2 = d.protind2num_pep(protind2);
  int *pepinds2 = d.protind2pepinds(protind2);
  
  int ix1 = 0;
  int ix2 = 0;
  overlap.clear();
  only1.clear();
  only2.clear();
  
  while(ix1 < num_pep1 && ix2 < num_pep2)
    {
      if(pepinds1[ix1] == pepinds2[ix2])
	{
	  overlap.push_back(pepinds1[ix1]);
	  ix1++; ix2++;
	  
	}
      else if(pepinds1[ix1] < pepinds2[ix2])
	{
	  only1.push_back(pepinds1[ix1]);
	  ix1++;
	}
      else
	{
	  only2.push_back(pepinds2[ix2]); ix2++;
	}
    }
  while(ix1 < num_pep1)
    {
      only1.push_back(pepinds1[ix1]);
      ix1++;
    }
  while(ix2 < num_pep2)
    {
      only2.push_back(pepinds2[ix1]);
      ix2++;
    }

  if(only1.size() == 0 && only2.size() == 0)
    return true;
  else
    //return (is_equivalent(d, protind1, only2) && is_equivalent(d, protind2, only1)); 
    return false;
}

void ProtScores :: get_intersection(vector<int> &set1, vector<int> &set2, vector<int> &result)
{
  result.clear();
  unsigned int ix1 = 0;
  unsigned int ix2 = 0;
  while(ix1 < set1.size() && ix2 < set2.size())
    {
      if(set1[ix1] == set2[ix2])
	{
	  result.push_back(set1[ix1]);
	  ix1++; ix2++;
	  
	}
      else if(set1[ix1] < set2[ix2])
	{
	  ix1++;
	}
      else
	{
	  ix2++;
	}
    }
}


void ProtScores :: get_complement(Dataset &d, int protind1, vector<int> &set2, vector<int> &result)
{
  int num_pep1 = d.protind2num_pep(protind1);
  int *pepinds1 = d.protind2pepinds(protind1);

  result.clear();
  int ix1 = 0;
  unsigned int ix2 = 0;
  while(ix1 < num_pep1 && ix2 < set2.size())
    {
      if(pepinds1[ix1] == set2[ix2])
	{
	  ix1++; ix2++;
	}
      else if(pepinds1[ix1] < set2[ix2])
	{
	  result.push_back(pepinds1[ix1]);
	  ix1++;
	}
      else
	{
	  ix2++;
	}
    }
  while(ix1 < num_pep1)
    {
      result.push_back(pepinds1[ix1]);
      ix1++;
    }
}


void ProtScores :: make_meta_set(Dataset &d)
{
  map<int,int>protind2pos_in_array;
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      int protind = scores[i].protind;
      protind2pos_in_array[protind] = i;
    }

  vector<int> processed;
  processed.resize(d.get_num_proteins(),0);

  int group_num = 0;
  for(unsigned int i = 0; i < scores.size();i++)
    {
      int protind1 = scores[i].protind;
      if(processed[protind1] == 1)
	continue;
      processed[protind1] = 1;
      group_num++;
      scores[i].group_number = group_num;
      
      int num_pep1 = d.protind2num_pep(protind1);
      int *pepinds = d.protind2pepinds(protind1);
      int flag = 0;
      for(int j = 0; j < num_pep1; j++)
	{
	  int pepind = pepinds[j];
	  int num_prot = d.pepind2num_prot(pepind);
	  int *protinds = d.pepind2protinds(pepind);
	  
	  for(int k = 0; k < num_prot; k++)
	    {
	      int protind2 = protinds[k];
	      int pos_in_array = protind2pos_in_array[protind2];
	      assert(scores[pos_in_array].protind == protind2);
	      if(protind2 == protind1)
		{
		  flag = 1;
		  continue;
		}
	      else if (processed[protind2] == 1)
		continue;
	      else
		{
		  if(is_subset(d,protind1, protind2) && is_subset(d,protind2,protind1)) 
		    {
		      processed[protind2] = 1;
		      (scores[i].indistinguishable_protinds).push_back(protind2);
		      scores[pos_in_array].group_number = group_num;
		      scores[pos_in_array].indistinguishable_prot = 1;
		    }
		}
	    }
	}
      assert(flag == 1);
    }
  
  cout << "number of meta groups " << group_num <<endl;
  
  //do some checks
  map<int, set<int> > protind2group;
  map<int, set<int> > group2protind;
  int count = 0;
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      int protind = scores[i].protind;
      int group = scores[i].group_number;
      //protind2group[protind].insert(group);
      //group2protind[group].insert(protind);
      if(scores[i].indistinguishable_prot == 0)
	{
	  count++;
	  vector<int> same_prot_list = scores[i].indistinguishable_protinds;
	  if(same_prot_list.size() > 0)
	    {
	      vector<int> same_prot_list = scores[i].indistinguishable_protinds;
	      
	      for(unsigned int k = 0; k < same_prot_list.size(); k++)
		{
		  int protind2 = same_prot_list[k];
		  int pos_in_array = protind2pos_in_array[protind2];
		  assert(scores[pos_in_array].protind == protind2);
		  assert(scores[pos_in_array].indistinguishable_prot == 1);
		  assert(scores[pos_in_array].group_number == group);
		
		  vector<int> overlap;
		  vector<int> only1;
		  vector<int> only2;
		  if(!are_equal(d, protind, protind2, overlap, only1, only2))
		    scores[i].has_complement = 1;
		}
	    }
	}
    }
    
  //make sure that each protein belongs only to a single group
  //for(map<int, set<int> > :: iterator it = protind2group.begin(); it != protind2group.end(); it++)
  //assert((it->second).size() == 1);
  
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      if(scores[i].indistinguishable_prot == 1)
	continue;
      
      int protind = scores[i].protind;
      //int group = scores[i].group_number;
      vector<int> same_prot_list = scores[i].indistinguishable_protinds;
      
      //do a couple more checks
      //if(same_prot_list.size() > 0)
      //{
      //  set<int> protinds_in_group = group2protind[group];
      //  assert(same_prot_list.size()+1 == protinds_in_group.size());
      //}

      if(scores[i].has_complement == 1)
	{
	  vector<int> intersection;
	  int num_pep = d.protind2num_pep(protind);
	  int *pepinds = d.protind2pepinds(protind);
	  for(int k = 0; k < num_pep; k++)
	    intersection.push_back(pepinds[k]);
	  vector<int> complement;
	  map<int, vector<int> > protind2complement;
	  for(unsigned int k = 0; k < same_prot_list.size(); k++)
	    {
	      vector<int> overlap;
	      vector<int> only1;
	      vector<int> only2;
	      int protind2 = same_prot_list[k];
	      are_equal(d, protind, protind2, overlap, only1, only2);
	      vector<int> result;
	      get_intersection(intersection, overlap, result);
	      intersection = result;
	    }
	  
	  for(unsigned int k = 0; k < same_prot_list.size(); k++)
	    {
	      int protind2 = same_prot_list[k];
	      get_complement(d, protind2, intersection, complement);
	      protind2complement[protind2] = complement;
	    }
	
	  get_complement(d, protind, intersection, complement);
	  protind2complement[protind] = complement;
	  
	  scores[i].intersection = intersection;
	  scores[i].protind2complement = protind2complement;
	}
    }
  

  //now find the subsets
  for(unsigned int i = 0; i < scores.size();i++)
    {
      //if this is a protein that is indistinguishable from some other, don't consider it
      if(scores[i].indistinguishable_prot == 1)
	continue;
      
      int protind1 = scores[i].protind;
      int num_pep1 = d.protind2num_pep(protind1);
      int *pepinds = d.protind2pepinds(protind1);
      
      int flag = 0;
      for(int j = 0; j < num_pep1; j++)
	{
	  int pepind = pepinds[j];
	  int num_prot = d.pepind2num_prot(pepind);
	  int *protinds = d.pepind2protinds(pepind);
	  
	  for(int k = 0; k < num_prot; k++)
	    {
	      int protind2 = protinds[k];
	      int pos_in_array = protind2pos_in_array[protind2];
	      assert(protind2 == scores[pos_in_array].protind);
	      int group = scores[pos_in_array].group_number;
	      int found = 0;
	      for(unsigned int k = 0; k < scores[i].parent_groups.size(); k++)
		{
		  if(scores[i].parent_groups[k] == group)
		    {
		      found = 1; 
		      break;
		    }
		}
	      if(found == 1)
		continue;
	      else if(protind2 == protind1)
		{
		  flag = 1;
		  continue;
		}
	      else if (scores[pos_in_array].indistinguishable_prot == 1)
		continue;
	      else
		{
		  // is protind1 subset of protind2
		  if(is_subset(d,protind2,protind1))
		    {
		      (scores[i].parent_groups).push_back(group);
		      scores[i].subset_prot=1;
		      (scores[pos_in_array].subset_protinds).push_back(protind1);
		    }
		
		}
	    }
	}
      assert(flag == 1);
    }

  
  //do some checking
  int num_unique = 0;
  int num_same = 0;
  int num_subset = 0;
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      if(scores[i].indistinguishable_prot == 1)
	num_same++;
      else if(scores[i].subset_prot == 1)
	num_subset++;
      else
	{
	  num_unique++;
	
	  //int protind = scores[i].protind;
	  //int group = scores[i].group_number;
	  //if(scores[i].subset_protinds.size() > 0)
	  //{
	  //  for(int j = 0; j < scores[i].subset_protinds.size(); j++)
	  //{
	  //  int protind2 = scores[i].subset_protinds[j];
	  //  assert(is_subset(d,protind, protind2));
	  //  int pos_in_array = protind2pos_in_array[protind2];
	  //  assert(scores[pos_in_array].subset_prot == 1);
	  //  int flag = 0;
	  //  for(int k = 0; k < scores[pos_in_array].parent_groups.size(); k++)
	  //    {
	  //      if(scores[pos_in_array].parent_groups[k] == group)
	  //	flag = 1;
	  //    }
	  //  assert(flag == 1);
	  //}
	}
      
    }
  
  cout << "number of subset groups " << num_subset << endl;

  ProtScoreHolder s;
  scores_meta.resize(num_unique,s);
  scores_same.resize(num_same,s);
  scores_subset.resize(num_subset,s);
  
  int idx_unique = 0;
  int idx_same = 0;
  int idx_subset = 0;

  for(unsigned int i = 0; i < scores.size(); i++)
    {
      if(scores[i].indistinguishable_prot == 1)
	{
	  scores_same[idx_same] = scores[i];
	  idx_same++;
	}
      else if(scores[i].subset_prot == 1)
	{
	  scores_subset[idx_subset] = scores[i];
	  idx_subset++;
	}
      else
	{
	  scores_meta[idx_unique] = scores[i];
	  idx_unique++;
	}
    }

  assert((idx_unique+idx_subset) == group_num);
  assert((idx_unique+idx_subset+idx_same) == (int)scores.size());

  scores.clear();
  scores = scores_meta;
  scores_meta.clear();

}



/*
bool ProtScores :: is_subset(Dataset &d, int protind1, int protind2)
{
  int num_pep1 = d.protind2num_pep(protind1);
  int *pepinds1 = d.protind2pepinds(protind1);
  int num_pep2 = d.protind2num_pep(protind2);
  int *pepinds2 = d.protind2pepinds(protind2);
  
  int ix1 = 0;
  int ix2 = 0;
  vector<int> overlap;
  vector<int> only1;
  vector<int> only2;
  
  while(ix1 < num_pep1 && ix2 < num_pep2)
    {
      if(pepinds1[ix1] == pepinds2[ix2])
	{
	  ix1++; ix2++;
	}
      else if(pepinds1[ix1] < pepinds2[ix2])
	ix1++;
      else
	{
	  only2.push_back(pepinds2[ix2]); ix2++;
	}
    }
  if(only2.size() == 0)
    return true;
  for(unsigned int i = 0; i < only2.size(); i++)
    {
      string pep2 = d.ind2pep(only2[i]);
      string p2;
      size_t pos = pep2.find('.');
      if(pos > 0 && pos != string::npos)
	p2 = pep2.substr(pos+1,pep2.size());
      pos = p2.find('.');
      if(pos > 0 && pos != string::npos)
	p2 = p2.substr(0,pos);
      
      for(int j = 0; j < num_pep1; j++)
	{
	  string pep1 = d.ind2pep(pepinds1[j]);
	  string p1;
	  pos = pep1.find('.');
	  if(pos > 0 && pos != string::npos)
	    p1 = pep1.substr(pos+1,pep1.size());
	  pos = p1.find('.');
	  if(pos > 0 && pos != string::npos)
	    p1 = p1.substr(0,pos);
	  
	  unsigned int ix1 = 0; unsigned int ix2 = 0;
	  while(ix1 < p1.size() && ix2 < p2.size())
	    {
	      char c1 = p1.at(ix1);
	      char c2 = p2.at(ix2);
	      if(c1 != c2)
		{
		  if(c1 == 'L' || c1 =='I')
		    {
		      if(c2 != 'L' && c2 !='I')
			return false; 
		    }
		  else if (c1 == 'T' || c1 =='S')
		    {
		      if(c2 != 'T' && c2 !='S')
			return false;
		    }
		  else
		    return false;
		}
	      ix1++; ix2++;
	    }
  	}
    }
  return true;
}



void ProtScores :: make_meta_set(Dataset &d)
{
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      int protind = scores[i].protind;
      int num_pep = d.protind2num_pep(protind);
      scores[i].score = (double) num_pep;
    }
  sort(scores.begin(),scores.end());
  reverse(scores.begin(),scores.end());

  vector<int> used_prot;
  used_prot.resize(d.get_num_proteins(),0);

  vector<int> processed;
  processed.resize(d.get_num_proteins(),0);

  for(unsigned int i = 0; i < scores.size();i++)
    {
      int protind1 = scores[i].protind;
      if(processed[protind1] == 1)
	continue;
      processed[protind1] = 1;

      int num_pep = d.protind2num_pep(protind1);
      int *pepinds = d.protind2pepinds(protind1);

      int flag = 0;
      for(int j = 0; j < num_pep; j++)
	{
	  int pepind = pepinds[j];
	  int num_prot = d.pepind2num_prot(pepind);
	  int *protinds = d.pepind2protinds(pepind);
	  
	  for(int k = 0; k < num_prot; k++)
	    {
	      int protind2 = protinds[k];
	      if(protind2 == protind1)
		flag = 1;
	      else if (used_prot[protind2] == 1)
		continue;
	      else if (processed[protind2] == 1)
		continue;
	      else
		{
		  //if(d.is_prot_subset(protind1, protind2)) 
		  if(is_subset(d,protind1, protind2)) 
		    {
		      processed[protind2] = 1;
		      used_prot[protind2] = 1;
		      (scores[i].subset_protinds).push_back(protind2);
		    }
		}
	    }
	}
      assert(flag == 1);
    }
  
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      int protind = scores[i].protind;
      if(used_prot[protind] == 1)
	{
	  //if(scores[i].subset_protinds.size()!=0)
	  //cout << " here " << scores[i].subset_protinds.size() << " " << i  << " " << scores[i].subset_protinds[0]<< endl;
	  scores[i].in_meta_prot = 1;
	}
    }
  
  for(unsigned int i = 0; i < scores.size(); i++)
    {
      for(unsigned int k = 0; k < scores[i].subset_protinds.size(); k++)
	{
	  int protind = scores[i].subset_protinds[k];
	  assert(used_prot[protind] == 1);
	  assert(d.is_prot_subset(scores[i].protind, protind));
	    
	}
    }
  
  if(0)
    {
      ProtScoreHolder s;
      //mix up the examples
      int n = scores.size();
      for(int i = 0; i < n; i++)
	{
	  int p1 = (int)((double)rand()/RAND_MAX*(n-1)); 
	  int p2 = (int)((double)rand()/RAND_MAX*(n-1)); 
	  s = scores[p1];
	  scores[p1] = scores[p2];
	  scores[p2] = s;
	}
    }
   
}
*/
