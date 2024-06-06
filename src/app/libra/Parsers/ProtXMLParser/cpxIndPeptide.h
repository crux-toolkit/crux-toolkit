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

#ifndef _CPXINDPEPTIDE_H
#define _CPXINDPEPTIDE_H

#include "cpxModInfo.h"
#include <string>
#include <vector>

class cpxIndPeptide{
public:

  //Constructors & Destructors
  cpxIndPeptide();
  cpxIndPeptide(const cpxIndPeptide& c);
  ~cpxIndPeptide();

  //Operators
  cpxIndPeptide& operator=(const cpxIndPeptide& c);
  //bool operator==(const CProteinAmbiguityGroup& c);

  //Data members
  int charge;
  double calcNeutralPepMass;
  std::string peptideSequence;
  cpxModInfo modificationInfo;

  //Functions
  //void writeOut(FILE* f, int tabs = -1);

private:
};

#endif