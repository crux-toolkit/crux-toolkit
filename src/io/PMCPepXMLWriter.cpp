#include "PMCPepXMLWriter.h"

using namespace Crux;

/**
 * Writes the data in a ProteinMatchCollection to the currently open file
 */
void PMCPepXMLWriter::write(
  ProteinMatchCollection* collection ///< collection to be written
) {
  if (!file_) {
    carp(CARP_FATAL, "No file open to write to.");
  }

  writeHeader();
  writePSMs(collection);
  writeFooter();
}

/**
 * Writes the PSMs in a ProteinMatchCollection to the currently open file
 */
void PMCPepXMLWriter::writePSMs(
  ProteinMatchCollection* collection ///< collection to be written
) {

  const map<pair<int, int>, int>& spectrum_counts = collection->getMatchesSpectrum();

  // iterate over matches
  for (SpectrumMatchIterator spec_iter = collection->spectrumMatchBegin();
       spec_iter != collection->spectrumMatchEnd();
       ++spec_iter) {
    SpectrumMatch* spec_match = *spec_iter;
    PeptideMatch* pep_match = spec_match->getPeptideMatch();

    vector<string> protein_names;
    vector<string> protein_descriptions;
    for (ProteinMatchIterator prot_iter = pep_match->proteinMatchBegin();
         prot_iter != pep_match->proteinMatchEnd();
         ++prot_iter) {
      ProteinMatch* prot_match = *prot_iter;
      Protein* protein = prot_match->getProtein();
      const char* tmp;
      tmp = protein->getIdPointer();
      protein_names.push_back((tmp != NULL) ? tmp : "");
      tmp = protein->getAnnotationPointer();
      protein_descriptions.push_back((tmp != NULL) ? tmp : "");
    }

    // populate scores
    double scores[NUMBER_SCORER_TYPES] = { 0 };
    bool scores_computed[NUMBER_SCORER_TYPES] = { false };
    for (ScoreMapIterator score_iter = spec_match->scoresBegin();
         score_iter != spec_match->scoresEnd();
         ++score_iter) {
      scores[score_iter->first] = score_iter->second;
      scores_computed[score_iter->first] = true;
    }

    // populate ranks
    int ranks[NUMBER_SCORER_TYPES] = { 0 };
    for (RankMapIterator rank_iter = spec_match->ranksBegin();
         rank_iter != spec_match->ranksEnd();
         ++rank_iter) {
      ranks[rank_iter->first] = rank_iter->second;
    }

    Spectrum* spectrum = spec_match->getSpectrum();
    SpectrumZState zstate = spec_match->getZState();
    Peptide* peptide = pep_match->getPeptide();
    int spec_scan = spectrum->getFirstScan();
    string spec_filename = spectrum->getFilename();
    FLOAT_T spec_neutral_mass = zstate.getNeutralMass();
    int spec_charge = zstate.getCharge();

    // get sequence and modified sequence
    char* seq;
    seq = peptide->getSequence();
    string seq_str(seq);
    free(seq);
    MODIFIED_AA_T* mod_seq = peptide->getModifiedAASequence();
    seq =
      modified_aa_string_to_string_with_masses(mod_seq, peptide->getLength(),
      get_mass_format_type_parameter("mod-mass-format"));
    string mod_seq_str(seq);
    free(seq);
    free(mod_seq);
    // get flanking aas
    char* flanking = peptide->getFlankingAAs();
    string flanking_str(flanking);
    free(flanking);
                                                             
    FLOAT_T peptide_mass = peptide->getPeptideMass();

    // write psm
    map<pair<int, int>, int>::const_iterator lookup =
      spectrum_counts.find(make_pair(spec_scan, spec_charge));
    writePSM(spec_scan, spec_filename.c_str(),
             spec_neutral_mass, spec_charge,
             ranks, seq_str.c_str(), mod_seq_str.c_str(), peptide_mass,
             protein_names.size(), flanking_str.c_str(),
             protein_names, protein_descriptions,
             scores_computed, scores,
             (lookup != spectrum_counts.end()) ? lookup->second : 0);
  }

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

