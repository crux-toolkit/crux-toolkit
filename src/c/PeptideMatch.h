#ifndef PEPTIDEMATCH_H_
#define PEPTIDEMATCH_H_

#include "AbstractMatch.h"
#include "match_objects.h"

class PeptideMatch : public AbstractMatch {

 protected:  
  std::vector<ProteinMatch*> protein_matches_;
  std::map<ProteinMatch*, PeptideSrc*> protein_match_to_peptide_src_;
  std::vector<SpectrumMatch*> spectrum_matches_;
  Crux::Peptide* peptide_;

 public:
  PeptideMatch();
  PeptideMatch(
    Crux::Peptide* peptide
  );

  ~PeptideMatch();

  void setPeptide(
    Crux::Peptide* peptide
  );

  Crux::Peptide* getPeptide();
   

  void addProteinMatch(
    ProteinMatch* protein_match, 
    PeptideSrc* src 
  );

  PeptideSrc* getSrc(
    ProteinMatch* protein_match
  );

  void addSpectrumMatch(
    SpectrumMatch* spectrum_match
  );
  
  SpectrumMatchIterator spectrumMatchBegin();
  SpectrumMatchIterator spectrumMatchEnd();

  ProteinMatchIterator proteinMatchBegin();
  ProteinMatchIterator proteinMatchEnd();
 
};

#endif
#ifndef PEPTIDEMATCH_H_
#define PEPTIDEMATCH_H_

#include "AbstractMatch.h"
#include "match_objects.h"

class PeptideMatch : public AbstractMatch {

 protected:  
  std::vector<ProteinMatch*> protein_matches_;
  std::map<ProteinMatch*, PeptideSrc*> protein_match_to_peptide_src_;
  std::vector<SpectrumMatch*> spectrum_matches_;
  Crux::Peptide* peptide_;

 public:
  PeptideMatch();
  PeptideMatch(
    Crux::Peptide* peptide
  );

  ~PeptideMatch();

  void setPeptide(
    Crux::Peptide* peptide
  );

  Crux::Peptide* getPeptide();
   

  void addProteinMatch(
    ProteinMatch* protein_match, 
    PeptideSrc* src 
  );

  PeptideSrc* getSrc(
    ProteinMatch* protein_match
  );

  void addSpectrumMatch(
    SpectrumMatch* spectrum_match
  );
  
  SpectrumMatchIterator spectrumMatchBegin();
  SpectrumMatchIterator spectrumMatchEnd();

  ProteinMatchIterator proteinMatchBegin();
  ProteinMatchIterator proteinMatchEnd();
 
};

#endif
