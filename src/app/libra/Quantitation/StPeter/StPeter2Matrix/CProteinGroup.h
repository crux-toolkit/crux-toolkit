#ifndef _CPROTEINGROUP_H
#define _CPROTEINGROUP_H

#include <string>
#include <vector>

//Simple class of three arrays.
//First array holds a set of protein identifiers contained in the group.
//Second and third arrays of length equal to the number of protXML files.
//The index array holds an internal processing index for efficiency - not a result.
//The quant array holds the actual quantification values of the protein group for each protXML.

class CProteinGroup {
public:

  void add(std::string prot);
  bool check(std::string& prot);
  void reset(size_t sz);

  bool validated;
  std::vector<std::string> proteins;
  std::vector<size_t> index;
  std::vector<double> quant;

private:

};

#endif