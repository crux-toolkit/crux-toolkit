#ifndef _CPEPXMLPEPTIDE_H
#define _CPEPXMLPEPTIDE_H

#include "PepXMLStructs.h"
#include <string>
#include <vector>

class CPepXMLPeptide {
public:

  CPepXMLPeptide();
  CPepXMLPeptide(const CPepXMLPeptide& p);
  ~CPepXMLPeptide();

  CPepXMLPeptide& operator=(const CPepXMLPeptide& p);
  
  void addMod(PepXMLPepMod& m);
  void clear();
  size_t sizeMod();
  size_t sizeOf();

  char label;  //0=alpha, 1=beta : for cross-links
  char nextAA;
  char prevAA;

  double calcPepNeutMass;
  double complementMass;

  std::string modifiedPeptide;
  std::string peptide;

  std::vector<PepXMLPepMod>*   mods;
  std::vector<size_t>*         proteins;
  std::vector<PepXMLPepScore>* xlScores;
  
};

#endif
