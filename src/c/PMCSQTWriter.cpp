#include "PMCSQTWriter.h"

/**
 * Writes the data in a ProteinMatchCollection to the currently open file
 */
void PMCSQTWriter::write(
  ProteinMatchCollection* collection, ///< collection to be written
  string database, ///< the database name
  int top_match ///< the top matches to output
) {
  if (!file_) {
    carp(CARP_FATAL, "No file open to write to.");
  }

  // count proteins
  int num_proteins = 0;
  for (ProteinMatchIterator iter = collection->proteinMatchBegin();
       iter != collection->proteinMatchEnd();
       ++iter) {
    ++num_proteins;
  }

  writeHeader(database, num_proteins);
  writePSMs(collection, top_match);
}

/**
 * Writes the PSMs in a ProteinMatchCollection to the currently open file
 */
void PMCSQTWriter::writePSMs(
  ProteinMatchCollection* collection, ///< collection to be written
  int top_match ///< the top matches to output
) {

  const map<pair<int, int>, int>& spectrum_counts = collection->getMatchesSpectrum();

  // iterate over matches
  for (SpectrumMatchIterator spec_iter = collection->spectrumMatchBegin();
       spec_iter != collection->spectrumMatchEnd();
       ++spec_iter) {
    SpectrumMatch* spec_match = *spec_iter;
    PeptideMatch* pep_match = spec_match->getPeptideMatch();

    // write spectrum
    Crux::Spectrum* spectrum = spec_match->getSpectrum();
    SpectrumZState z_state = spec_match->getZState();
    map<pair<int, int>, int>::const_iterator lookup =
      spectrum_counts.find(make_pair(spectrum->getFirstScan(), z_state.getCharge()));
    writeSpectrum(spectrum, z_state, (lookup != spectrum_counts.end()) ?
                  lookup->second : 0);

    FLOAT_T xcorr_score = -1.0;
    int xcorr_rank = -1;
    FLOAT_T sp_score = -1.0;
    int sp_rank = -1;
    FLOAT_T delta_cn = -1.0;
    FLOAT_T by_ions_matched = -1;
    FLOAT_T by_ions_total = -1;

    // xcorr
    if (spec_match->hasScore(XCORR)) {
      xcorr_score = spec_match->getScore(XCORR);
    } else {
      carp_once(CARP_WARNING, "Missing XCorr value writing SQT.");
    }
    if (spec_match->hasRank(XCORR)) {
      xcorr_rank = spec_match->getRank(XCORR);
      if (xcorr_rank > top_match) {
        continue;
      }
    } else {
      carp_once(CARP_WARNING, "Missing XCorr rank writing SQT.");
    }

    // sp
    if (spec_match->hasScore(SP)) {
      sp_score = spec_match->getScore(SP);
    } else {
      carp_once(CARP_WARNING, "Missing SP value writing SQT.");
    }
    if (spec_match->hasRank(SP)) {
      sp_rank = spec_match->getRank(SP);
    } else {
      carp_once(CARP_WARNING, "Missing SP rank writing SQT.");
    }

    // deltacn
    if (spec_match->hasScore(DELTA_CN)) {
      delta_cn = spec_match->getScore(DELTA_CN);
    } else {
      carp_once(CARP_WARNING, "Missing DeltaCN value writing SQT.");
    }

    // b/y ions
    if (spec_match->hasScore(BY_IONS_MATCHED)) {
      by_ions_matched = spec_match->getScore(BY_IONS_MATCHED);
    } else {
      carp_once(CARP_WARNING, "Missing B/Y ions matched value writing SQT.");
    }
    if (spec_match->hasScore(BY_IONS_TOTAL)) {
      by_ions_total = spec_match->getScore(BY_IONS_TOTAL);
    } else {
      carp_once(CARP_WARNING, "Missing B/Y ions total value writing SQT.");
    }

    bool is_decoy = false; // TODO not sure how to determine

    // write psm
    writePSM(pep_match->getPeptide(),
             xcorr_score, xcorr_rank,
             sp_score, sp_rank,
             delta_cn,
             by_ions_matched, by_ions_total,
             is_decoy);
  }
 
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

