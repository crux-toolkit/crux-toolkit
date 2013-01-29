#ifndef PROTEINMATCH_H_
#define PROTEINMATCH_H_

#include "match_objects.h"
#include "AbstractMatch.h"

#include <map>
#include <string>
#include <vector>

class ProteinMatch : public AbstractMatch {

 protected:
  std::vector<PeptideMatch*> peptide_matches_;
  Crux::Protein* protein_;
  
 public:
  ProteinMatch();
  ProteinMatch(
    Crux::Protein* protein
  );
  
  ~ProteinMatch();

  PeptideMatchIterator peptideMatchBegin();
  PeptideMatchIterator peptideMatchEnd();

  void addPeptideMatch(
    PeptideMatch*
  );

  std::string getId();

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
