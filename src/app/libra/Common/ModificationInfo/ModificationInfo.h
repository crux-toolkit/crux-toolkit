#ifndef MODIFICATION_H
#define MODIFICATION_H

/*

Program       : ModificationInfo                                                    
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

Primary data object holding all mixture distributions for each precursor ion charge

Copyright (C) 2003 Andrew Keller

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
akeller@systemsbiology.org

*/

#include <stdio.h>

#include "Common/sysdepend.h"
#include "Common/Array.h"
#include "Parsers/Parser/Tag.h"
#include "Parsers/Algorithm2XML/SearchParams/SearchParams.h"
#include "Quantitation/Option.h"
#include "Common/ResidueMass/ResidueMass.h"
class Modification; // forward reference

class ModificationInfo {
 public:

  ModificationInfo(double nterm, double cterm);
  ModificationInfo(const char* embedded_pep, Array<Modification*>* mods);
  ModificationInfo(Array<Tag*>* tags);
  ModificationInfo(const ModificationInfo* modinfo);
  ~ModificationInfo();
  //  double getModifiedPeptideMass(char* peptide, Boolean monoisotopic);
  void enterAAMod(int pos, double mass);
  double getNtermModMass() const;
  double getCtermModMass() const;
  Boolean isModified() const;
  int getNumModAAs() const;
  int getModAAPos(int index) const;
  double getModAAMass(int index) const;
  Array<Tag*>* getModificationInfoInfoTags(char* std_peptide_name) const;
  
  char* getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error) const;
  char* getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error, const char* starttag, const char* endtag, bool out_stat) const;
  char* getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error, const char* starttag, const char* endtag) const;
  char* getStandardModifiedPeptide(const char* peptide, Array<StaticModificationCount>* static_consts, double error, bool out_stat) const;
  Boolean isModifiedResidue(int pos) const;
  double getModifiedResidueMass(int pos) const;
  Boolean equivalentModification(const ModificationInfo* modinfo, double error, const char* peptide, const char* quant_labels) const;
  char* getModifiedPeptide() const;
  Boolean setQueryString(char* query, int maxlength);

 protected:

  double mod_nterm_mass_;
  double mod_cterm_mass_;
  Array<int>* mod_aa_positions_;
  Array<double>* mod_aa_masses_;
  char* mod_peptide_;
};


#endif
