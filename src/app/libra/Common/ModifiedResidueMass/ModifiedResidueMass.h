#ifndef MOD_RESIDUE_H
#define MOD_RESIDUE_H

#include <stdlib.h>
#include <string.h>

#include "Common/constants.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Common/ResidueMass/ResidueMass.h"

class ModifiedResidueMass {

 public:

  ModifiedResidueMass();

  static double getUnmodifiedPeptideMass(char* pep, Boolean monoisotopic);
  static double getModifiedPeptideMass(ModificationInfo* modinfo, char* pep, Boolean monoisotopic);


 protected:


};


#endif
