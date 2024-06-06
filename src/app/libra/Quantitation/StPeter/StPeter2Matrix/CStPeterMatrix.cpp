#include "CStPeterMatrix.h"

using namespace std;

CStPeterMatrix::CStPeterMatrix(){
  outputFileName.clear();
  protFilter.clear();
}

CProteinGroup& CStPeterMatrix::operator[](const size_t& index){
  return matrix[index];
}

//add protXML file names to name array.
void CStPeterMatrix::addName(string name){
  names.push_back(name);
}

//Iterates over all proteins in a protXML file for a specific groupID
void CStPeterMatrix::addProteinsToGroup(CProteinGroup& g, ProtXMLParser& p, size_t index, int groupID){
  size_t i, k;
  for (i = 0; i<p.size(); i++){
    if(p[i].probability<0) continue; //don't count proteins already counted
    if (p[i].stPeter.dSIn == 0 && p[i].stPeter.SIn==0) continue; //don't count proteins that weren't quantified
    if (p[i].groupID == groupID){
      g.add(p[i].proteinName);
      for (k = 0; k<p[i].altProteinNames->size(); k++) g.add(p[i].altProteinNames->at(k));
      if(p[i].stPeter.dSIn!=0){
        g.quant[index]+=pow(2,p[i].stPeter.dSIn);
      } else if (p[i].stPeter.SIn != 0) {
        g.quant[index] += pow(2, p[i].stPeter.SIn);
      }
      p[i].probability-=2; //mark this protein as analyzed
    }
  }
}

void CStPeterMatrix::crossFilter(ProtXMLParser& p){
  int i;
  size_t j,k;
  double prob=p.getProbability(0.01);
  cout << "Cross-filter probability: " << prob << endl;
  if(prob<0.9) {
    cout << "Setting probability to 0.9" << endl;
    prob=0.9;
  }
  for(i=0;i<p.size();i++){
    if(p[i].probability>=prob){
      for(j=0;j<matrix.size();j++){
        if(matrix[j].check(p[i].proteinName)){
          matrix[j].validated=true;
          break;
        }
        for(k=0;k<p[i].altProteinNames->size();k++){
          if (matrix[j].check(p[i].altProteinNames->at(k))){
            matrix[j].validated = true;
            break;
          }
        }
      }
    }
  }

  for(j=0;j<matrix.size();j++){
    if(matrix[j].validated) continue;
    matrix.erase(matrix.begin()+j);
    j--;
  }

}

//If a decoy filter is applied, remove all protein names that qualify as decoys
void CStPeterMatrix::decoyFilter(){
  if(protFilter.size()==0) return; //only filter results if the filter is set.

  size_t i,j;
  for(i=0;i<matrix.size();i++){
    for(j=0;j<matrix[i].proteins.size();j++){
      if(matrix[i].proteins[j].find(protFilter)!=string::npos){
        matrix[i].proteins.erase(matrix[i].proteins.begin()+j);
        j--;
      }
    }

    //if a protein group was all decoys, remove that protein group from the matrix.
    if(matrix[i].proteins.size()==0){
      matrix.erase(matrix.begin()+i);
      i--;
    }
  }
}

//Write the matrix to file. Format is groupID (a number), proteins, and sample values.
bool CStPeterMatrix::exportMatrix(){
  size_t i,j;
  
  if(outputFileName.size()==0) {
    cout << "ERROR: no file name set for output" << endl;
    return false;
  }

  FILE* f=fopen(outputFileName.c_str(),"wt");
  if(f==NULL){
    cout << "ERROR: cannot open " << outputFileName << " for output" << endl;
    return false;
  }
  fprintf(f,"ID\tProteins");
  for(i=0;i<names.size();i++) fprintf(f,"\t%s",names[i].c_str());
  fprintf(f,"\n");
  for(i=0;i<matrix.size();i++){
    fprintf(f,"%d\t",(int)i+1);
    for(j=0;j<matrix[i].proteins.size();j++){
      if(j>0) fprintf(f,",");
      fprintf(f,"%s",matrix[i].proteins[j].c_str());
    }
    for(j=0;j<matrix[i].quant.size();j++){
      if(matrix[i].quant[j]==0) fprintf(f,"\tNA");
      else fprintf(f,"\t%.6lf",log2(matrix[i].quant[j]));
    }
    fprintf(f,"\n");
  }
  fclose(f);
  return true;
}

//Main function of class. Accepts an array of protXML classes and a list of all indexed proteins
//contained in those classes. Then iterates over every protXML class to group overlapping protein
//groups into supergroups and quantifies them, creating a matrix summarizing all protXML files.
void CStPeterMatrix::extractStPeter(vector<ProtXMLParser>& p, CProteinMasterList& ml){
  size_t i, j, k;
  CProteinGroup group;

  //iterate over all protXML
  for (i = 0; i<p.size(); i++){

    //iterate over all proteins in the protXML
    for (j = 0; j<p[i].size(); j++){
      if (p[i][j].stPeter.dSIn == 0 && p[i][j].stPeter.SIn == 0) continue; //skip any proteins that haven't been quantified
      if (p[i][j].probability<0) continue; //skip any proteins already processed into a group

      //create new group and add all other proteins with the same groupID
      group.reset(p.size());
      addProteinsToGroup(group, p[i], i, p[i][j].groupID);
      group.index[i] = group.proteins.size();

      //iterate over other samples to add all their proteins that
      //match to any protein in this group
      for(k=i+1;k<p.size();k++){
        matchGroups(p[k],group,ml,k);
        group.index[k] = group.proteins.size();
      }

      //If the other protXML files extended the group
      //revisit in each protXML file to add the new entries.
      //Iterate over each until satisfied (i.e. no new proteins).
      while (true){
        bool bFinished = true;
        for (k = 0; k<group.index.size(); k++){
          if (group.index[k]<group.proteins.size()){
            bFinished = false;
            //search for new proteins
            matchGroups(p[k],group,ml,k);
            group.index[k]=group.proteins.size();
          }
        }
        if (bFinished) break;
      }

      //Group is now complete. Add group to matrix.
      matrix.push_back(group);
    }
  }
}

//Iterates over a protXML to find all proteins that match the current set of groupIDs in a matrix group.
void CStPeterMatrix::matchGroups(ProtXMLParser& p, CProteinGroup& g, CProteinMasterList& ml, size_t mlIndex){
  size_t i,j;
  int index;

  vector<int> finishedIDs; //keep track of the groupIDs already searched so they aren't repeated

  //iterate over the newly added proteins
  for(i=g.index[mlIndex];i<g.proteins.size();i++){

    //find the index of the protein in the protXML of interest
    j=ml.findProtein(g.proteins[i]);
    if(j==SIZE_MAX) continue;
    index=ml[j].index[mlIndex];
    if(index==-1) continue;

    //see if we've already processed the group of this protein
    for(j=0;j<finishedIDs.size();j++){
      if(p[index].groupID==finishedIDs[j]) break;
    }
    if(j==finishedIDs.size()) {
      //if not, add this protein to the group (and any other proteins with the same groupID).
      addProteinsToGroup(g, p, mlIndex, p[index].groupID);
      finishedIDs.push_back(p[index].groupID);
    }
  }
}

void CStPeterMatrix::setFilter(const char* filter){
  protFilter = filter;
}

void CStPeterMatrix::setFilter(std::string filter){
  setFilter(filter.c_str());
}

void CStPeterMatrix::setOutput(const char* fname){
  outputFileName=fname;
}

void CStPeterMatrix::setOutput(std::string fname){
  setOutput(fname.c_str());
}

size_t CStPeterMatrix::size(){
  return matrix.size();
}

