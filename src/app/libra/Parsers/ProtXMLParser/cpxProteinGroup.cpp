#include "cpxProteinGroup.h"

using namespace std;

cpxProteinGroup::cpxProteinGroup(){
  groupNumber = 0;
  probability = 0;
  protein = new vector<cpxProtein>;
}

cpxProteinGroup::cpxProteinGroup(const cpxProteinGroup& c){
  groupNumber = c.groupNumber;
  probability = c.probability;
  protein = new vector<cpxProtein>;
  for (size_t i = 0; i<c.protein->size(); i++) protein->push_back(c.protein->at(i));
}

cpxProteinGroup::~cpxProteinGroup(){
  delete protein;
}

cpxProtein& cpxProteinGroup::operator[](const size_t index){
  return protein->at(index);
}

cpxProteinGroup& cpxProteinGroup::operator=(const cpxProteinGroup& c){
  if (this != &c){
    groupNumber = c.groupNumber;
    probability = c.probability;
    delete protein;
    protein = new vector<cpxProtein>;
    for (size_t i = 0; i<c.protein->size(); i++) protein->push_back(c.protein->at(i));
  }
  return *this;
}

size_t cpxProteinGroup::size(){
  return protein->size();
}
