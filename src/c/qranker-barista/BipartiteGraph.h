#ifndef BIPARTITEGRAPH_H
#define BIPARTITEGRAPH_H
#include<iostream>
#include<fstream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
using namespace std;

class BipartiteGraph{
 public:
  struct Range {int p; int len;};
 BipartiteGraph():nranges(0),nindices(0), ranges((Range*)0), indices((int*)0){}
  ~BipartiteGraph(){clear();}
  void clear(){delete[] ranges; ranges = 0; delete[] indices; indices = 0; nranges=0; nindices=0;}
  void create_bipartite_graph(map<int, set<int> > data);

  bool is_index_in_range(int index, int range);
  //assumes that inds are sorted
  bool is_subset(int r1, int r2);
  int get_nranges() const {return nranges;}
  int get_nindices() const {return nindices;}
  int get_range_length(int r){return ranges[r].len;}
  int* get_range_indices(int r){return (indices+ranges[r].p);}

  void save(ofstream &os);
  void load(ifstream &is);
 private:
  int nranges; //how many ranges
  int nindices; //size of the index array
  Range* ranges;
  int* indices;
};


#endif //BIPARTITEGRAPH_H
