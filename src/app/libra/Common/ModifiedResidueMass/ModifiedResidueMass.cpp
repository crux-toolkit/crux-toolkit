
#include "ModifiedResidueMass.h"


double ModifiedResidueMass::getUnmodifiedPeptideMass(char* pep, Boolean monoisotopic) {
  return ResidueMass::getProteinMass(pep, monoisotopic);

  /*
  double total = ResidueMass::getMass('n', monoisotopic) + getMass('c', monoisotopic);
  if(pep != NULL) {
    for(int k = 0; k < strlen(pep); k++)
      total += ResidueMass::getMass(pep[k], monoisotopic);
  }
  return total;
  */
}

double ModifiedResidueMass::getModifiedPeptideMass(ModificationInfo* modinfo, char* pep, Boolean monoisotopic) {
  if(modinfo == NULL)
    return getUnmodifiedPeptideMass(pep, monoisotopic);

  double mass = 0.0;
  if(modinfo->getNtermModMass() > 0.0)
    mass += modinfo->getNtermModMass();
  else
    mass += ResidueMass::getMass('n', monoisotopic);
  if(modinfo->getCtermModMass() > 0.0)
    mass += modinfo->getCtermModMass();
  else
    mass += ResidueMass::getMass('c', monoisotopic);
  int index = 0;
  if(pep != NULL) {
    for(int k = 0; k < (int) strlen(pep); k++)
      if(index < modinfo->getNumModAAs() && modinfo->getModAAPos(index) == k + 1)
	mass += modinfo->getModAAMass(index++); 
      else
	mass += ResidueMass::getMass(pep[k], monoisotopic);
  }
  return mass;

}

