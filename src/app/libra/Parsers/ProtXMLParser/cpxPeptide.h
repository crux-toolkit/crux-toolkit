/*
Copyright 2017, Michael R. Hoopmann, Institute for Systems Biology
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef _CPXPEPTIDE_H
#define _CPXPEPTIDE_H

#include "cpxIndPeptide.h"
#include "cpxModInfo.h"
#include <string>
#include <vector>

class cpxPeptide{
public:

  //Constructors & Destructors
  cpxPeptide();
  cpxPeptide(const cpxPeptide& c);
  ~cpxPeptide();

  //Operators
  cpxIndPeptide& operator[](const size_t index);
  cpxPeptide& operator=(const cpxPeptide& c);
  //bool operator==(const CProteinAmbiguityGroup& c);

  //Data members
  int charge;
  double expTotInstances;
  double fpkmAdjustedProbability;
  double initialProbability;
  bool isContributingEvidence;
  bool isNondegenerateEvidence;
  cpxModInfo modificationInfo;
  int nEnzymaticTermini;
  int nInstances;
  double nspAdjustedProbability;
  double nSiblingPeptides;
  int nSiblingPeptidesBin;
  std::string peptideSequence;
  double weight;
  std::vector<cpxIndPeptide>* indistinguishablePeptide;

  //Functions
  size_t size();
  //void writeOut(FILE* f, int tabs = -1);

private:
};

#endif