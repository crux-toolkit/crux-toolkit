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

#ifndef _CPXMODINFO_H
#define _CPXMODINFO_H

#include <string>
#include <vector>

typedef struct spxMod{
  int pos;
  double mass;
} spxMod;

class cpxModInfo{
public:

  //Constructors & Destructors
  cpxModInfo();
  cpxModInfo(const cpxModInfo& c);
  ~cpxModInfo();

  //Operators
  cpxModInfo& operator=(const cpxModInfo& c);
  //bool operator==(const CProteinAmbiguityGroup& c);

  //Data members
  std::string modifiedPeptide;
  std::vector<spxMod>* modAminoAcidMass;
  double nTermMod;
  double cTermMod;

  //Functions
  //void writeOut(FILE* f, int tabs = -1);

private:
};

#endif