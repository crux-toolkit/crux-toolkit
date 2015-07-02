#include "BipartiteGraph.h"

void BipartiteGraph::create_bipartite_graph(map<int, set<int> > data)
{
  int n = data.size();
  //allocate the range array
  ranges = new Range[n];
  nranges = 0;
  nindices = 0;
  int max_len = 0;
  for(map <int, set<int> >::iterator it = data.begin(); it != data.end(); it++)
    {
      int len = 0;
      for(set<int>::iterator itp = (it->second).begin(); itp != (it->second).end(); itp++)
	{
	  len++;
	  nindices++;
	}
      if(len > max_len)
	max_len = len;
      int ind = it->first;
      ranges[ind].len = len;
      nranges++;
    }
  indices = new int[nindices];

  vector<int> t;
  t.reserve(max_len);
  set<int> s;
  int pos = 0;
  for (int i = 0; i < nranges; i++)
    {
      s = data[i];
      t.erase(t.begin(),t.end());
      for(set<int>::iterator itp = s.begin(); itp != s.end(); itp++)
	t.push_back(*itp);
      sort(t.begin(),t.end());
      assert((int)t.size() == ranges[i].len);
      ranges[i].p = pos;
      for (unsigned int k = 0; k < t.size(); k++)
	{
	  indices[pos] = t[k];
	  pos++;
	}
    }
}


void BipartiteGraph::save(ofstream &os)
{
  os.write((char*)(&nranges),sizeof(int));
  os.write((char*)(&nindices),sizeof(int));
  os.write((char*)ranges,sizeof(Range)*nranges);
  os.write((char*)indices,sizeof(int)*nindices);
}

void BipartiteGraph::load(ifstream &is)
{
  is.read((char*)(&nranges),sizeof(int));
  is.read((char*)(&nindices),sizeof(int));
  ranges = new Range[nranges];
  indices = new int[nindices];
  is.read((char*)ranges,sizeof(Range)*nranges);
  is.read((char*)indices,sizeof(int)*nindices);

}


bool BipartiteGraph::is_index_in_range(int index, int range)
{
  assert(range < nranges);
  for(int i = 0; i < ranges[range].len; i++)
    {
      int ofst = ranges[range].p+i;
      if(indices[ofst] == index)
	return true;
    }
  return false;
}

bool BipartiteGraph :: is_subset(int r1, int r2)
{
  assert(r1 < nranges);
  assert(r2 < nranges);
  int len1 = ranges[r1].len;
  int len2 = ranges[r2].len;
  int *inds1 = indices+ranges[r1].p;
  int *inds2 = indices+ranges[r2].p;
  int i1 = 0;
  int i2 = 0;
  
  while(i1 < len1 && i2 < len2)
    {
      if(inds1[i1] == inds2[i2])
	{
	  i1++;
	  i2++;
	}
      else if(inds1[i1] < inds2[i2])
	i1++;
      else
	return false;
    }
  return true;
}

