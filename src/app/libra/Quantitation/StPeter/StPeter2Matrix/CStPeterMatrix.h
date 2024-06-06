#ifndef _CSTPETERMATRIX_H
#define _CSTPETERMATRIX_H

#include "CProteinGroup.h"
#include "CProteinMasterList.h"
#include "Quantitation/StPeter/ProtXMLParser.h"

//This is the matrix class. The matrix is just a vector of ProteinGroups.
//The class also holds parameters and support functions for manipulating
//and exporting the matrix.

class CStPeterMatrix{
public:
  CStPeterMatrix();
  CProteinGroup& operator[](const size_t& index);
  
  void addName(std::string name);
  void crossFilter(ProtXMLParser& p);
  void decoyFilter();
  void extractStPeter(vector<ProtXMLParser>& p, CProteinMasterList& ml);
  bool exportMatrix();

  void setFilter(const char* filter);
  void setFilter(std::string filter);

  void setOutput(const char* fname);
  void setOutput(std::string fname);

  size_t size();

private:

  void addProteinsToGroup(CProteinGroup& g, ProtXMLParser& p, size_t index, int groupID);
  void matchGroups(ProtXMLParser& p, CProteinGroup& g, CProteinMasterList& ml, size_t mlIndex);

  std::string outputFileName;
  std::string protFilter;
  std::vector<CProteinGroup> matrix;
  std::vector<std::string> names;
};

#endif
