#ifndef _CPROTEINMASTERLIST_H
#define _CPROTEINMASTERLIST_H

#include "Quantitation/StPeter/ProtXMLParser.h"
#include <algorithm>
#include <string>
#include <vector>

//simple matrix whose values consist of a protein identifier,
//and the protein indexes in each corresponding protXML file.
typedef struct ProtMLValue{
  std::string protein;
  std::vector<int> index;
} ProtMLValue;

class CProteinMasterList{
public:
  
  ProtMLValue& operator[](const size_t& index);

  void addProtXML(ProtXMLParser& p, std::string id);
  size_t findProtein(std::string& prot);

  size_t size();
  size_t sizeSamples();

private:

  std::vector<ProtMLValue> masterList;
  std::vector<std::string> protXMLID;
  size_t protXMLCount;

  static bool compareProtein(const ProtMLValue& a, const ProtMLValue& b);
};

#endif
