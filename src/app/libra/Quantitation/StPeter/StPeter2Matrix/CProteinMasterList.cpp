#include "CProteinMasterList.h"

using namespace std;

//Return the protein and its list of indexes
ProtMLValue& CProteinMasterList::operator[](const size_t& index){
  return masterList[index];
}

//Add a new protXML file to the data matrix
void CProteinMasterList::addProtXML(ProtXMLParser& p, string id){
  size_t i,j,k;
  ProtMLValue v;
  
  //add our id (typically protXML file name), and indexes to existing array
  protXMLID.push_back(id);
  for(i=0;i<masterList.size();i++){
    masterList[i].index.push_back(-1); //default index is negative to indicate no protein
  }

  //set up new entry stub
  for (i = 0; i<protXMLID.size();i++) v.index.push_back(-1);

  //iterate over all proteins in the current protXML
  for(i=0;i<p.size();i++){
    j=findProtein(p[i].proteinName);
    if(j==SIZE_MAX){ //add mew protein if protein not in matrix yet
      v.protein=p[i].proteinName;
      v.index.back()=i;

      //after adding the protein, resort the list every 100th entry
      masterList.push_back(v);
      if (masterList.size() % 100 == 0) sort(masterList.begin(), masterList.end(), compareProtein);
    } else { //or simply set the index if protein is in matrix
      masterList[j].index.back()=i;
    }

    //iterater over alternate proteins
    //altertnate protein names get their own entries in the matrix, but the indexes
    //may or may not point to the same location as other proteins in the matrix.
    for(k=0;k<p[i].altProteinNames->size();k++){
      j=findProtein(p[i].altProteinNames->at(k));
      if (j == SIZE_MAX){ //add new protein
        v.protein = p[i].altProteinNames->at(k);
        v.index.back() = i;

        //after adding the protein, resort the list every 100th entry
        masterList.push_back(v);
        if (masterList.size() % 100 == 0) sort(masterList.begin(), masterList.end(), compareProtein);
      } else { //or set the index
        masterList[j].index.back() = i;
      }
    }
  }

}

//Find if a protein is already in the matrix. Use a binary search, however sorting only
//occurs every 100th addition to the array, so must iteratively search the last [up to 99] entries.
size_t CProteinMasterList::findProtein(std::string& prot){
  // Find if peptide already listed by binary search
  size_t sz = masterList.size() - masterList.size() % 100; //buffer of 100 unsorted entries
  size_t lower = 0;
  size_t mid = sz / 2;
  size_t upper = sz;
  int i;

  if (sz>0){
    i = masterList[mid].protein.compare(prot);
    while (i != 0){
      if (lower >= upper) break;
      if (i>0){
        if (mid == 0) break;
        upper = mid - 1;
        mid = (lower + upper) / 2;
      } else {
        lower = mid + 1;
        mid = (lower + upper) / 2;
      }
      if (mid == sz) break;
      i = masterList[mid].protein.compare(prot);
    }

    //match by peptide sequence, so check modifications
    if (i == 0) return mid;
  }

  //check unsorted proteins
  for (mid = sz; mid < masterList.size(); mid++){
    if (masterList[mid].protein.compare(prot)==0) return mid;
  }

  //SIZE_MAX indicates protein not found
  return SIZE_MAX;

}

size_t CProteinMasterList::size(){
  return masterList.size();
}

size_t CProteinMasterList::sizeSamples(){
  return protXMLID.size();
}

bool CProteinMasterList::compareProtein(const ProtMLValue& a, const ProtMLValue& b){
  return (a.protein.compare(b.protein)<0);
}