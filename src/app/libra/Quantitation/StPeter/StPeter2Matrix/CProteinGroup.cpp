#include "CProteinGroup.h"

using namespace std;

//Add a new protein name to the group.
//Checks if protein name is already in the group.
void CProteinGroup::add(string prot){
  size_t i;
  for(i=0;i<proteins.size();i++){
    if(proteins[i].compare(prot)==0) break;
  }
  if(i==proteins.size()) proteins.push_back(prot);
}

//Check if the protein group contains the requested name.
bool CProteinGroup::check(string& prot){
  size_t i;
  for (i = 0; i<proteins.size(); i++){
    if (proteins[i].compare(prot) == 0) return true;
  }
  return false;
}

//Clears out all protein names and sets index and quant values
//to zeroes.
void CProteinGroup::reset(size_t sz){
  validated=false;
  proteins.clear();
  index.clear();
  quant.clear();
  for(size_t i=0;i<sz;i++) {
    index.push_back(0);
    quant.push_back(0);
  }
}
