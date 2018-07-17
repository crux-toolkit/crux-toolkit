/**
* \file AssignConfidenceApplication.h
* AUTHOR: Attila Kertesz-Farkas
* CREATE DATE: May 01, 2015
****************************************************************************/

#include "AssignConfidenceApplication.h"
#include "ComputeQValues.h"
#include "io/MatchCollectionParser.h"
#include "PosteriorEstimator.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"

#include <map>
#include <utility>

using namespace std;
using namespace Crux;

static const int MAX_PSMS = 10000000;
// 14th decimal place
static const double EPSILON = 0.00000000000001;

/**
* \returns a blank ComputeQValues object
*/
AssignConfidenceApplication::AssignConfidenceApplication():
  spectrum_flag_(NULL), iteration_cnt_(0) {
}

/**
* Destructor
*/
AssignConfidenceApplication::~AssignConfidenceApplication() {
}

/**
 * Generate keys when building a hash on PSMs.
 * http://stackoverflow.com/questions/98153/whats-the-best-hashing-algorithm-to-use-on-a-stl-string-when-using-hash-map/
 */
int stringToIndex(string myString) {
  int returnValue = 0;
  for(string::iterator it = myString.begin(); it != myString.end(); ++it) {
    returnValue = (returnValue * 101) + (int)*it;
  }
  return(returnValue);
}

/**
* main method for ComputeQValues
*/

int AssignConfidenceApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("target input"));
}

int AssignConfidenceApplication::main(const vector<string>& input_files) {
  // Prepare the output files if not in Cascade Search
  if (spectrum_flag_ == NULL) {
    output_ = new OutputFiles(this);
  }
  
  ESTIMATION_METHOD_T estimation_method;
  string method_param = Params::GetString("estimation-method");
  carp(CARP_INFO, "Estimation method = %s.", method_param.c_str());
  if (method_param == "tdc") {
    estimation_method = TDC_METHOD;
  } else if (method_param == "mix-max") {
    estimation_method = MIXMAX_METHOD;
  } else if (method_param == "peptide-level") {
    estimation_method = PEPTIDE_LEVEL_METHOD;
  } else {
    carp(CARP_FATAL, "The estimation method \"%s\" is not supported.", method_param.c_str());
  }
    
  // Perform the analysis.
  if (input_files.empty()) {
    carp(CARP_FATAL, "No input files found.");
  }
  if (estimation_method != PEPTIDE_LEVEL_METHOD &&
      Params::GetBool("combine-modified-peptides")) {
    carp(CARP_WARNING, "The \"combine-modified-peptides\" option is ignored when estimation-method is not peptide-level.");
  }
  bool sidak = Params::GetBool("sidak");

  int top_match = 1;
  if (estimation_method == PEPTIDE_LEVEL_METHOD) {
    top_match = MAX_PSMS+1;
  }

  SCORER_TYPE_T score_type = INVALID_SCORER_TYPE;
  SCORER_TYPE_T derived_score_type = INVALID_SCORER_TYPE;

  string score_param = Params::GetString("score");
  if (!score_param.empty()) {
    switch (get_column_idx(score_param.c_str())) {
    case SP_SCORE_COL:
      score_type = SP;    ///< SEQUEST preliminary score
      break;
    case XCORR_SCORE_COL:
      score_type = XCORR;   ///< SEQUEST primary score
      break;
    case EVALUE_COL:
      score_type = EVALUE;  ///< Comet e-value
      break;
    case PVALUE_COL:        ///< Search-for-xlinks p-value
      score_type = LOGP_BONF_WEIBULL_XCORR;
      break;
    case PERCOLATOR_SCORE_COL:
      score_type = PERCOLATOR_SCORE;
      break;
    case PERCOLATOR_QVALUE_COL:
      score_type = PERCOLATOR_QVALUE;
      break;
    case PERCOLATOR_PEP_COL:
      score_type = PERCOLATOR_PEP;
      break;
    case QRANKER_SCORE_COL:
      score_type = QRANKER_SCORE;
      break;
    case QRANKER_QVALUE_COL:
      score_type = QRANKER_QVALUE;
      break;
    case QRANKER_PEP_COL:
      score_type = QRANKER_PEP;
      break;
    case BARISTA_SCORE_COL:
      score_type = BARISTA_SCORE;
      break;
    case BARISTA_QVALUE_COL:
      score_type = BARISTA_QVALUE;
      break;
    case EXACT_PVALUE_COL:
      score_type = TIDE_SEARCH_EXACT_PVAL;
      break;
    case REFACTORED_SCORE_COL:
      score_type = TIDE_SEARCH_REFACTORED_XCORR;
      break;
    case RESIDUE_PVALUE_COL:
      score_type = RESIDUE_EVIDENCE_PVAL;
      break;
    case RESIDUE_EVIDENCE_COL:
      score_type = RESIDUE_EVIDENCE_SCORE;
      break;
    case BOTH_PVALUE_COL:
      score_type = BOTH_PVALUE;
      break;
    case ELUTION_WINDOW_COL:
      score_type = TIDE_SEARCH_EXACT_SMOOTHED;
      break;      
    default:
      carp(CARP_FATAL, "The PSM feature \"%s\" is not supported.", score_param.c_str());
    }
  }

  if (sidak && score_type != TIDE_SEARCH_EXACT_PVAL) {
    carp(CARP_WARNING, "Sidak adjustment may not be compatible with score: %s", score_param.c_str());
  }

  // Create two match collections, for targets and decoys.
  MatchCollection* target_matches = new MatchCollection();
  map<int, MatchCollection*> decoy_matches; // key is decoy index

  bool ascending, distinct_matches;
  MatchCollectionParser parser;
  map<string, FLOAT_T> BestPeptideScore;

  bool avgTdc = estimation_method == TDC_METHOD;
  for (vector<string>::const_iterator iter = input_files.begin(); iter != input_files.end(); ++iter) {
    string target_path = *iter;
    string decoy_path = *iter;

    if (target_path.find("decoy") != string::npos) {
      carp(CARP_FATAL, "%s appears to be a decoy file. Only target or concatenated files "
                       "should be given to assign-confidence because it automatically searches for "
                       "corresponding decoy files.", target_path.c_str());
    }

    check_target_decoy_files(target_path, decoy_path);

    if (!FileUtils::Exists(target_path)) {
      carp(CARP_FATAL, "Target file %s not found", target_path.c_str());
    } else if (!FileUtils::Exists(decoy_path)) {
      if (estimation_method == MIXMAX_METHOD) {
        carp(CARP_FATAL, "Cannot find file %s. Decoy file from separate target-decoy search is "
                         "required for mix-max q-value calculation", decoy_path.c_str());
      }
      carp(CARP_DEBUG, "Decoy file %s not found", decoy_path.c_str());
      decoy_path = "";
    }

    MatchCollection* match_collection = parser.create(target_path, Params::GetString("protein-database"));
    distinct_matches = match_collection->getHasDistinctMatches();
    if (!match_collection->hasDecoyIndexes()) {
      avgTdc = false;
    }

    carp(CARP_INFO, "Found %d PSMs in %s.", match_collection->getMatchTotal(), target_path.c_str());

    // If necessary, automatically identify the score type.
    // The score type that is used is the first one found
    // in the list below.
    if (score_type == INVALID_SCORER_TYPE) {
      vector<SCORER_TYPE_T> scoreTypes;
      scoreTypes.push_back(XCORR);
      scoreTypes.push_back(EVALUE);
      scoreTypes.push_back(BOTH_PVALUE);
      scoreTypes.push_back(RESIDUE_EVIDENCE_PVAL);
      scoreTypes.push_back(TIDE_SEARCH_EXACT_PVAL);
      scoreTypes.push_back(TIDE_SEARCH_EXACT_SMOOTHED);
      scoreTypes.push_back(LOGP_BONF_WEIBULL_XCORR);
      scoreTypes.push_back(PERCOLATOR_SCORE);
      for (vector<SCORER_TYPE_T>::const_iterator i = scoreTypes.begin(); i != scoreTypes.end(); i++) {
        if (match_collection->getScoredType(*i)) {
          score_type = *i;
          carp(CARP_INFO, "Automatically detected score type: %s", scorer_type_to_string(score_type));
          break;
        }
      }
      if (score_type == INVALID_SCORER_TYPE) {
        carp(CARP_FATAL, "Could not detect score type. Specify the score type using the \"score\" parameter.");
      }
    }
    switch (getDirection(score_type)) {
      case -1:
        ascending = false;
        break;
      case 1:
        ascending = true;
        break;
      default:
        carp(CARP_FATAL, "Cannot infer sort order for score %s.", scorer_type_to_string(score_type));
    }
    carp(CARP_INFO, "Score type=%s, sorting in %s order",
         scorer_type_to_string(score_type), ascending ? "ascending" : "descending");

    if (!match_collection->getScoredType(score_type)) {
      const char* score_str = scorer_type_to_string(score_type);
      carp(CARP_FATAL, "The PSM feature \"%s\" was not found in file \"%s\".", score_str, target_path.c_str());
    }

    // Find and keep the best score for each peptide.
    if (estimation_method == PEPTIDE_LEVEL_METHOD) {
      peptide_level_filtering(match_collection, &BestPeptideScore, score_type, ascending);
      carp(CARP_INFO, "%d distinct target peptides.", BestPeptideScore.size());
    }

    target_matches->setScoredType(score_type, match_collection->getScoredType(score_type));
    target_matches->setScoredType(EVALUE, match_collection->getScoredType(EVALUE));
    target_matches->setScoredType(DELTA_CN, match_collection->getScoredType(DELTA_CN));
    target_matches->setScoredType(SP, match_collection->getScoredType(SP));
    target_matches->setScoredType(BY_IONS_MATCHED, match_collection->getScoredType(BY_IONS_MATCHED));
    target_matches->setScoredType(BY_IONS_TOTAL, match_collection->getScoredType(BY_IONS_TOTAL));
    target_matches->setScoredType(SIDAK_ADJUSTED, sidak);
    for (map<int, MatchCollection*>::iterator i = decoy_matches.begin(); i != decoy_matches.end(); i++) {
      i->second->setScoredType(SIDAK_ADJUSTED, sidak);
    }

    // Counters just to let the user know what's up.
    int num_target_rank_skipped = 0;
    int num_decoy_rank_skipped = 0;
    int num_target_peptide_skipped = 0;
    int num_decoy_peptide_skipped = 0;
    
    if (decoy_path != "") {
      MatchCollection* temp_collection = parser.create(decoy_path, Params::GetString("protein-database"));
      carp(CARP_INFO, "Found %d PSMs in %s.", temp_collection->getMatchTotal(), decoy_path.c_str());

      if (temp_collection->hasDecoyIndexes()) {
        avgTdc = true;
        MatchIterator* decoy_iter = new MatchIterator(temp_collection);
        while (decoy_iter->hasNext()) {
          Crux::Match* decoy_match = decoy_iter->next();
          match_collection->addMatch(decoy_match);
        }
        delete decoy_iter;
      } else {
        // Mark decoy matches
        // key = (filename, scan number, charge, rank); value = index
        map<boost::tuple <int, int, int, int>, int> pairidx;
        int cnt = 0;
        MatchIterator* temp_iter = new MatchIterator(temp_collection);
        while (temp_iter->hasNext()) {
          Crux::Match* decoy_match = temp_iter->next();
          cnt++;

          // Only use top-ranked matches.
          if (decoy_match->getRank(XCORR) > top_match) {
            num_decoy_rank_skipped++;
            continue;
          }

          decoy_match->setNullPeptide(true);
          switch (estimation_method) {
          case MIXMAX_METHOD:
            {
              // Put match directly in the final set of decoys, because no TDC.
              int decoyIndex = decoy_match->decoyIndex();
              map<int, MatchCollection*>::iterator lookup = decoy_matches.find(decoyIndex);
              if (lookup == decoy_matches.end()) {
                decoy_matches[decoyIndex] = new MatchCollection();
                lookup = decoy_matches.find(decoyIndex);
              }
              lookup->second->addMatch(decoy_match);
            }
            break;
          case TDC_METHOD:
          case PEPTIDE_LEVEL_METHOD:
            {
              int fileIndex = stringToIndex(decoy_match->getSpectrum()->getFullFilename());
              int scanid = decoy_match->getSpectrum()->getFirstScan();
              int charge = decoy_match->getCharge();
              int rank = decoy_match->getRank(XCORR);
              boost::tuple<int, int, int, int> t(fileIndex, scanid, charge, rank);

              // If the PSM is already there, that means there was a tie
              // for top-ranked decoys.  In that case, there is no need to
              // store a pointer to the second one.
              if (pairidx[t] == 0) {
                pairidx[t] = cnt;
              }
            }
            break;
          case NUMBER_METHOD_TYPES:
          case INVALID_METHOD:
            carp(CARP_FATAL, "No estimation method specified.");
          }
        }
        delete temp_iter;

        // Find and keep the best score for each decoy peptide.
        if (estimation_method == PEPTIDE_LEVEL_METHOD) {
          peptide_level_filtering(temp_collection, &BestPeptideScore, score_type, ascending);
          carp(CARP_INFO, "%d distinct target+decoy peptides.", BestPeptideScore.size());
        }

        if (estimation_method != MIXMAX_METHOD) {
          int numCompetitions = 0;
          int numLostDecoys = 0;
          int numTies = 0;
          MatchCollection* tdc_collection = new MatchCollection();
          tdc_collection->setScoredType(score_type, true);
          MatchIterator* target_iter = new MatchIterator(match_collection);
          MatchIterator* decoy_iter = new MatchIterator(temp_collection);
          while (target_iter->hasNext()) {
            Crux::Match* target_match = target_iter->next();

            // Only use top-ranked matches.
            if (target_match->getRank(XCORR) > top_match) {
              num_target_rank_skipped++;
              continue;
            }

            // Retrieve the index of the corresponding decoy PSM.
            int fileIndex = stringToIndex(target_match->getSpectrum()->getFullFilename());
            int scanid = target_match->getSpectrum()->getFirstScan();
            int charge = target_match->getCharge();
            int rank = target_match->getRank(XCORR);
            int decoy_idx = pairidx[boost::tuple <int, int, int, int>(fileIndex, scanid, charge, rank)];
            if (decoy_idx == 0) {
              carp(CARP_DEBUG, "Failed to find decoy for file=%s scan=%d charge=%d rank=%d.",
                   target_match->getSpectrum()->getFullFilename(), scanid, charge, rank);
              numLostDecoys++;
            }

            if (estimation_method == PEPTIDE_LEVEL_METHOD) {
              if (decoy_idx == 0) {
                int numCandidates = target_match->getTargetExperimentSize();
                target_match->setTargetExperimentSize(numCandidates);
                tdc_collection->addMatch(target_match);
              } else {
                Crux::Match* decoy_match = decoy_iter->getMatch(decoy_idx - 1);
                int numCandidates = target_match->getTargetExperimentSize() + decoy_match->getTargetExperimentSize();
                target_match->setTargetExperimentSize(numCandidates);
                decoy_match->setTargetExperimentSize(numCandidates);
                tdc_collection->addMatch(target_match);
                tdc_collection->addMatch(decoy_match);
              }
            } else {
              if (decoy_idx == 0) {
                tdc_collection->addMatch(target_match);
                continue;
              }
              Crux::Match* decoy_match = decoy_iter->getMatch(decoy_idx - 1);
              int numCandidates = target_match->getTargetExperimentSize() + decoy_match->getTargetExperimentSize();
              target_match->setTargetExperimentSize(numCandidates);
              decoy_match->setTargetExperimentSize(numCandidates);

              // This is where the target-decoy competition happens.
              carp(CARP_DEBUG, "TDC: Comparing target (%d, +%d) with score %g to decoy (%d, +%d) with score %g.",
                   target_match->getSpectrum()->getFirstScan(), target_match->getCharge(), target_match->getScore(score_type),
                   decoy_match->getSpectrum()->getFirstScan(), decoy_match->getCharge(), decoy_match->getScore(score_type));

              FLOAT_T score_difference = target_match->getScore(score_type) - decoy_match->getScore(score_type);
              numCompetitions++;
              // Randomly break ties.
              if (fabs(score_difference) < 1e-10) {
                numTies++;
                score_difference += 0.5 - ((double)myrandom() / UNIFORM_INT_DISTRIBUTION_MAX);
              }
              if (ascending) { // smaller scores are better
                score_difference *= -1.0;
              }
              if (score_difference >= 0.0) {
                tdc_collection->addMatch(target_match);
              } else {
                tdc_collection->addMatch(decoy_match);
              }
            }
          }
          delete target_iter;
          delete decoy_iter;
          delete match_collection;
          match_collection = tdc_collection;
          carp(CARP_INFO, "%d tdc_collection", match_collection->getMatchTotal());
          if (numCompetitions > 0) {
            carp(CARP_INFO, "Randomly broke %d ties in %d target-decoy competitions.", numTies, numCompetitions);
          }
          if (numLostDecoys > 0) {
            carp(CARP_INFO, "Failed to find %d decoys.", numLostDecoys);
          }
        }
      }
      delete temp_collection;
    }

    // Iterate, gathering matches into one or two collections.
    MatchIterator* match_iterator = new MatchIterator(match_collection, score_type, false);
    while (match_iterator->hasNext()) {
      Match* match = match_iterator->next();
      bool is_decoy = match->getNullPeptide();

      // Only use top-ranked matches.
      if (score_type == BOTH_PVALUE) {
        if (match->getRank(BOTH_PVALUE) > top_match) {
          if (is_decoy) {
            num_decoy_rank_skipped++;
          } else {
            num_target_rank_skipped++;
          }
          continue;
        }
      } else if (score_type == RESIDUE_EVIDENCE_PVAL) {
        if (match->getRank(RESIDUE_EVIDENCE_PVAL) > top_match) {
          if (is_decoy) {
            num_decoy_rank_skipped++;
          } else {
            num_target_rank_skipped++;
          }
          continue;
        }
      } else {
        if (match->getRank(XCORR) > top_match) {
          if (is_decoy) {
            num_decoy_rank_skipped++;
          } else {
            num_target_rank_skipped++;
          }
          continue;
        }
      }

      // Find and keep the best score for each decoy peptide.
      if (estimation_method == PEPTIDE_LEVEL_METHOD) {
        FLOAT_T score = match->getScore(score_type);
        string peptideStr = getPeptideSeq(match);

        FLOAT_T bestScore;
        try {
          bestScore = BestPeptideScore.at(peptideStr);
          if (bestScore != score) {  //not the best scoring peptide
            if (is_decoy) {
              num_decoy_peptide_skipped++;
            } else {
              num_target_peptide_skipped++;              
            }
            continue;
          } else {
            BestPeptideScore.at(peptideStr) += ascending ? -1.0 : 1.0;  //make sure only one best scoring peptide reported.
          }
        } catch (const out_of_range& oor) {
          carp(CARP_DEBUG, "Error in peptide-level filtering");
        }
      }

      // Do the Sidak correction.
      if (sidak) {
        if (match->getRank(XCORR) > 1) {
          carp_once(CARP_WARNING, "Sidak correction is not defined for non-top-matches. Further warnings are not shown.");
        }
        double sidak_adjustment = 1.0 - pow(1.0 - match->getScore(score_type), match->getTargetExperimentSize());
        match->setScore(SIDAK_ADJUSTED, sidak_adjustment);
      }

      // Add this match to one of the collections.
      if (is_decoy) {
        int decoyIndex = match->decoyIndex();
        map<int, MatchCollection*>::iterator lookup = decoy_matches.find(decoyIndex);
        if (lookup == decoy_matches.end()) {
          decoy_matches[decoyIndex] = new MatchCollection();
          lookup = decoy_matches.find(decoyIndex);
        }
        lookup->second->addMatch(match);
      } else {
        target_matches->addMatch(match);
      }
      Match::freeMatch(match);
    }
    delete match_iterator;
    delete match_collection;
    if (num_decoy_rank_skipped + num_target_rank_skipped > 0) {
      carp(CARP_INFO, "Skipped %d target and %d decoy PSMs with rank > %d.",
           num_target_rank_skipped, num_decoy_rank_skipped, top_match);
    }
    if (num_target_peptide_skipped + num_decoy_peptide_skipped > 0) {
      carp(CARP_INFO, "Skipped %d target and %d decoy PSMs due to peptide-level filtering.",
           num_target_peptide_skipped, num_decoy_peptide_skipped);
    }
  }

  if (sidak) {
    score_type = SIDAK_ADJUSTED;
  }

  target_matches->setScoredType(score_type, true);
  for (map<int, MatchCollection*>::iterator i = decoy_matches.begin(); i != decoy_matches.end(); i++) {
    i->second->setScoredType(score_type, true);
  }

  // get from the input files which columns to print in the output files
  if (iteration_cnt_ == 0) {
    vector<bool> cols_to_print(NUMBER_MATCH_COLUMNS);
    cols_to_print[FILE_COL] = Params::GetBool("file-column");
    cols_to_print[SCAN_COL] = true;
    cols_to_print[CHARGE_COL] = true;
    cols_to_print[SPECTRUM_PRECURSOR_MZ_COL] = true;
    cols_to_print[SPECTRUM_NEUTRAL_MASS_COL] = true;
    cols_to_print[PEPTIDE_MASS_COL] = true;
    cols_to_print[DELTA_CN_COL] = target_matches->getScoredType(DELTA_CN);
    cols_to_print[SP_SCORE_COL] = target_matches->getScoredType(SP);
    cols_to_print[SP_RANK_COL] = target_matches->getScoredType(SP);

    if (score_type == BOTH_PVALUE) {
      cols_to_print[BOTH_PVALUE_COL] = true;
      cols_to_print[BOTH_PVALUE_RANK] = true;
    } else {
      cols_to_print[XCORR_SCORE_COL] = !target_matches->getScoredType(TIDE_SEARCH_EXACT_PVAL);
      cols_to_print[XCORR_RANK_COL] = true;
      cols_to_print[EVALUE_COL] = target_matches->getScoredType(EVALUE);
      cols_to_print[EXACT_PVALUE_COL] = target_matches->getScoredType(TIDE_SEARCH_EXACT_PVAL);
      cols_to_print[PVALUE_COL] = target_matches->getScoredType(LOGP_BONF_WEIBULL_XCORR);
      cols_to_print[SIDAK_ADJUSTED_COL] = sidak;
      if (target_matches->getScoredType(TIDE_SEARCH_EXACT_PVAL)) {
        cols_to_print[REFACTORED_SCORE_COL] = true;
      }
    }

    cols_to_print[BY_IONS_MATCHED_COL] = target_matches->getScoredType(BY_IONS_MATCHED);
    cols_to_print[BY_IONS_TOTAL_COL] = target_matches->getScoredType(BY_IONS_TOTAL);

    if (distinct_matches) {
      cols_to_print[DISTINCT_MATCHES_SPECTRUM_COL] = true;
    } else {
      cols_to_print[MATCHES_SPECTRUM_COL] = true;
    }

    switch (estimation_method) {
    case TDC_METHOD:
    case PEPTIDE_LEVEL_METHOD: // FIXME: Make a peptide-level q-value column. --WSN 4 Feb 2016
      cols_to_print[QVALUE_TDC_COL] = true;
      break;
    case MIXMAX_METHOD:
      cols_to_print[QVALUE_MIXMAX_COL] = true;
      break;
    case NUMBER_METHOD_TYPES:
    case INVALID_METHOD:
      carp(CARP_FATAL, "No estimation method specified.");
    }
    cols_to_print[SEQUENCE_COL] = true;
    cols_to_print[CLEAVAGE_TYPE_COL] = true;
    cols_to_print[PROTEIN_ID_COL] = true;
    cols_to_print[FLANKING_AA_COL] = true;
    if (spectrum_flag_ != NULL) {
      cols_to_print[INDEX_NAME_COL] = true;
    }

    output_->writeHeaders(cols_to_print);
  }
  switch (estimation_method) {
  case TDC_METHOD:
  case PEPTIDE_LEVEL_METHOD:
    derived_score_type = QVALUE_TDC;
    break;
  case MIXMAX_METHOD:
    derived_score_type = QVALUE_MIXMAX;
    break;
  case NUMBER_METHOD_TYPES:
  case INVALID_METHOD:
    carp(CARP_FATAL, "No estimation method specified.");
  }

  vector<FLOAT_T> target_scores;

  // Compute q-values.
  vector<FLOAT_T> qvalues;
  switch (estimation_method) {
  case TDC_METHOD:
    if (avgTdc) {
      carp(CARP_INFO, "Using a-TDC (%d decoy sets).", decoy_matches.size());
      vector< vector<FLOAT_T> > decoy_scores;
      AtdcScoreSet::getScores(target_matches, decoy_matches, score_type, ascending, target_scores, decoy_scores);
      qvalues = AtdcScoreSet(target_scores, decoy_scores, ascending).fdps();
      convert_fdr_to_qvalue(qvalues);
      break;
    }
  case PEPTIDE_LEVEL_METHOD:
    {
      target_scores = target_matches->extractScores(score_type);
      vector<FLOAT_T> decoy_scores;
      for (map<int, MatchCollection*>::const_iterator i = decoy_matches.begin(); i != decoy_matches.end(); i++) {
        vector<FLOAT_T> curScores = i->second->extractScores(score_type);
        copy(curScores.begin(), curScores.end(), back_inserter(decoy_scores));
      }
      carp(CARP_INFO, "There are %d target and %d decoy PSMs for q-value computation.",
           target_scores.size(), decoy_scores.size());
      qvalues = compute_decoy_qvalues_tdc(target_scores, decoy_scores, ascending, 1.0);
    }
    break;
  case MIXMAX_METHOD:
    {
      target_scores = target_matches->extractScores(score_type);
      vector<FLOAT_T> decoy_scores;
      for (map<int, MatchCollection*>::const_iterator i = decoy_matches.begin(); i != decoy_matches.end(); i++) {
        vector<FLOAT_T> curScores = i->second->extractScores(score_type);
        copy(curScores.begin(), curScores.end(), back_inserter(decoy_scores));
      }
      carp(CARP_INFO, "There are %d target and %d decoy PSMs for q-value computation.",
           target_scores.size(), decoy_scores.size());
      qvalues = compute_decoy_qvalues_mixmax(target_scores, decoy_scores, ascending, Params::GetDouble("pi-zero"));
    }
    break;
  default:
      carp(CARP_FATAL, "No estimation method specified.");
  }

  unsigned int fdr1 = 0;
  unsigned int fdr5 = 0;
  unsigned int fdr10 = 0;
  for (vector<FLOAT_T>::const_iterator i = qvalues.begin(); i != qvalues.end(); i++) {
    if (*i < 0.01) ++fdr1;
    if (*i < 0.05) ++fdr5;
    if (*i < 0.10) ++fdr10;
  }
  carp(CARP_INFO, "Number of PSMs at 1%% FDR = %d.", fdr1);
  carp(CARP_INFO, "Number of PSMs at 5%% FDR = %d.", fdr5);
  carp(CARP_INFO, "Number of PSMs at 10%% FDR = %d.", fdr10);

  // Store p-values to q-values as a hash, and then assign them.
  map<FLOAT_T, FLOAT_T> qvalue_hash = store_arrays_as_hash(target_scores, qvalues);

  target_matches->assignQValues(&qvalue_hash, score_type, derived_score_type);

  // Store targets by score.
  target_matches->sort(score_type);
  if (spectrum_flag_ == NULL) {
    output_->writeMatches(target_matches);
    output_->writeFooters();
    delete output_;
  } else {
    accepted_psms_ = 0;
    //print out accepted matches in Cascade Search
    FLOAT_T qValueThreshold = Params::GetDouble("q-value-threshold");
    if (is_final_) {
      qValueThreshold = 1.0;
    }

    MatchCollection accepted_matches;
    MatchIterator match_iterator(target_matches, score_type, false);
    while (match_iterator.hasNext()) {
      Match* match = match_iterator.next();
      if (match->getScore(QVALUE_TDC) > qValueThreshold) {
        break;
      }
      spectrum_flag_->insert(make_pair(pair<string, unsigned int>(
        match->getSpectrum()->getFullFilename(),
        match->getSpectrum()->getFirstScan() * 10 + match->getCharge()), true));

      match->setDatabaseIndexName(index_name_);

      accepted_matches.addMatch(match);
      ++accepted_psms_;
    }
    carp(CARP_INFO, "Accepted %d PSMs at q < %g.", accepted_matches.getMatchTotal(), qValueThreshold);
    output_->writeMatches(&accepted_matches);
  }
  delete target_matches;
  for (map<int, MatchCollection*>::const_iterator i = decoy_matches.begin(); i != decoy_matches.end(); i++) {
    delete i->second;
  }

  return 0;
} // Main


/**
* Find the best-scoring match for each peptide in a given collection.
* Only consider the top-ranked PSM per spectrum.
*
* Results are stored in the given match collection.
*/
void AssignConfidenceApplication::identify_best_psm_per_peptide(
  MatchCollection* all_matches,
  SCORER_TYPE_T score_type
) {
  /* Instantiate a hash table.  key = peptide; value = maximal xcorr
     for that peptide. */
  map<string, FLOAT_T> best_score_per_peptide;

  // Store in the hash the best score per peptide.
  MatchIterator* match_iterator 
    = new MatchIterator(all_matches, score_type, false);
  while(match_iterator->hasNext()) {
    Match* match = match_iterator->next();

    // Skip matches that are not top-ranked.
    if (match->getRank(score_type) == 1) {
      char *peptide = match->getModSequenceStrWithSymbols();
      FLOAT_T this_score = match->getScore(score_type);

      map<string, FLOAT_T>::iterator map_position 
        = best_score_per_peptide.find(peptide);

      if (map_position == best_score_per_peptide.end()) {
        best_score_per_peptide[peptide] = this_score;
      } else {
        // FIXME: Need a generic compare operator for score_type.
        if (map_position->second < this_score) {
          best_score_per_peptide[peptide] = this_score;
        }
      }
      free(peptide);
    }
  }
  delete match_iterator;

  // Set the best_per_peptide Boolean in the match, based on the hash.
  match_iterator = new MatchIterator(all_matches, score_type, false);
  while(match_iterator->hasNext()) {
    Match* match = match_iterator->next();

     // Skip matches that are not top-ranked.
    if (match->getRank(score_type) == 1) {
      char* peptide = match->getModSequenceStrWithSymbols();
      FLOAT_T this_score = match->getScore(score_type);

      map<string, FLOAT_T>::iterator map_position 
        = best_score_per_peptide.find(peptide);

      if (map_position->second == this_score) {
        match->setBestPerPeptide();
        
        // Prevent ties from causing two peptides to be best.
        best_score_per_peptide[peptide] = HUGE_VAL;
      }
      
      free(peptide);
    }
  }
  delete match_iterator;
}


/**
 * The q-value is defined as the minimum FDR at which a given score is
 * deemed significant.  This function takes a list of FDRs and
 * converts them into q-values.  The FDRs should be ordered from
 * lowest to highest, sorted according to the underlying score.
 */
void AssignConfidenceApplication::convert_fdr_to_qvalue(
  vector<FLOAT_T>& qvalues ///< Come in as FDRs, go out as q-values.
) {
  FLOAT_T prev_fdr = qvalues[qvalues.size() - 1];
  for (int idx = qvalues.size() - 2; idx >= 0; idx--) {
    carp(CARP_DETAILED_DEBUG, "fdr[%i] = %.10f", idx, qvalues[idx]);
    FLOAT_T this_fdr = qvalues[idx];
    if (prev_fdr < this_fdr) {
      qvalues[idx] = prev_fdr;
    }
    prev_fdr = qvalues[idx];
    carp(CARP_DETAILED_DEBUG, "qvalue[%i] = %.10f", idx, qvalues[idx]);
  }
}

/**
 * Store two parallel arrays of floats in a hash table.
 */
map<FLOAT_T, FLOAT_T> AssignConfidenceApplication::store_arrays_as_hash(
  const vector<FLOAT_T>& keys,
  const vector<FLOAT_T>& values
) {
  map<FLOAT_T, FLOAT_T> return_value;
  for (size_t i = 0; i < keys.size(); i++) {
    carp(CARP_DETAILED_DEBUG, "%g maps to %g", keys[i], values[i]);
    return_value[keys[i]] = values[i];
  }
  return return_value;
}

AssignConfidenceApplication::AtdcScoreSet::AtdcScoreSet(
  const vector<FLOAT_T>& targetScores,
  const vector< vector<FLOAT_T> >& decoyScores,
  bool ascending
) {
  for (vector< vector<FLOAT_T> >::const_iterator i = decoyScores.begin(); i != decoyScores.end(); i++) {
    if (i->size() != targetScores.size()) {
      carp(CARP_FATAL, "a decoy set contained %d scores when %d were expected", i->size(), targetScores.size());
    }
  }

  for (size_t i = 0; i < targetScores.size(); i++) {
    FLOAT_T targetScore = ascending ? targetScores[i] : -targetScores[i];
    if (i > 0 && targetScore < scores_.back().first) {
      carp(CARP_FATAL, "internal error; unsorted scores were passed into AtdcScoreSet");
    }
    scores_.push_back(make_pair(targetScore, vector<FLOAT_T>(decoyScores.size(), 0)));
    for (size_t j = 0; j < decoyScores.size(); j++) {
      scores_.back().second[j] = ascending ? decoyScores[j][i] : -decoyScores[j][i];
    }
  }
}

void AssignConfidenceApplication::AtdcScoreSet::getScores(
  MatchCollection* targets,
  map<int, MatchCollection*> decoys,
  SCORER_TYPE_T scoreType,
  bool ascending,
  vector<FLOAT_T>& outTargetScores,
  vector< vector<FLOAT_T> >& outDecoyScores
) {
  outTargetScores.clear();
  outDecoyScores.clear();

  vector< boost::tuple<FLOAT_T, int, int, int> > targetScores; // score, file, scan, charge
  for (MatchIterator i(targets, scoreType, false); i.hasNext(); ) {
    Match* match = i.next();
    targetScores.push_back(boost::make_tuple(
      match->getScore(scoreType),
      match->getFileIndex(),
      match->getSpectrum()->getFirstScan(),
      match->getCharge()));
  }

  if (ascending) {
    sort(targetScores.begin(), targetScores.end(), sortScoresAsc);
  } else {
    sort(targetScores.begin(), targetScores.end(), sortScoresDesc);
  }

  map< boost::tuple<int, int, int>, pair<size_t, size_t> > idxMap; // <file, scan, charge> -> <idx, decoys found>
  for (size_t i = 0; i < targetScores.size(); i++) {
    const boost::tuple<FLOAT_T, int, int, int>& score = targetScores[i];
    outTargetScores.push_back(score.get<0>());
    int file = score.get<1>();
    int scan = score.get<2>();
    int charge = score.get<3>();
    boost::tuple<int, int, int> idxKey = boost::make_tuple(file, scan, charge);
    if (!idxMap.insert(make_pair(idxKey, make_pair(i, 0))).second) {
      carp(CARP_FATAL, "Multiple target scores found for file %d, scan %d, charge %d", file, scan, charge);
    }
  }

  for (map<int, MatchCollection*>::const_iterator i = decoys.begin(); i != decoys.end(); i++) {
    outDecoyScores.push_back(vector<FLOAT_T>(outTargetScores.size(), numeric_limits<FLOAT_T>::quiet_NaN()));
    for (MatchIterator j(i->second, scoreType, false); j.hasNext(); ) {
      Match* match = j.next();
      boost::tuple<int, int, int> idxKey = boost::make_tuple(
        match->getFileIndex(), match->getSpectrum()->getFirstScan(), match->getCharge());
      map< boost::tuple<int, int, int>, pair<size_t, size_t> >::iterator k = idxMap.find(idxKey);
      if (k != idxMap.end()) {
        k->second.second++;
        outDecoyScores.back()[k->second.first] = match->getScore(scoreType);
      }
    }
  }

  for (map< boost::tuple<int, int, int>, pair<size_t, size_t> >::const_iterator i = idxMap.begin();
       i != idxMap.end();
       i++) {
    size_t decoysFound = i->second.second;
    if (decoys.size() > decoysFound) {
      carp(CARP_FATAL, "Missing %d decoy scores for file %d, scan %d, charge %d (expected %d, found %d)",
           decoys.size() - decoysFound, i->first.get<0>(), i->first.get<1>(), i->first.get<2>(),
           decoys.size(), decoysFound);
    } else if (decoysFound > decoys.size()) {
      carp(CARP_FATAL, "Found %d extra decoy scores for file %d, scan %d, charge %d (expectede %d, found %d)",
           decoysFound - decoys.size(), i->first.get<0>(), i->first.get<1>(), i->first.get<2>(),
           decoys.size(), decoysFound);
    }
  }
}

vector<FLOAT_T> AssignConfidenceApplication::AtdcScoreSet::fdps() const {
  const size_t numScores = scores_.size();
  const size_t numDecoySets = scores_.front().second.size();
  const FLOAT_T bc1 = -1/(FLOAT_T)numDecoySets;

  vector<size_t> optIdxCnt(numScores, 0);
  vector<FLOAT_T> sumNtds(numScores, 0);
  vector<FLOAT_T> sumNdds(numScores, 0);
  for (size_t i = 0; i < numDecoySets; i++) {
    vector<int> ntds(numScores, 0);
    vector<int> ndds(numScores, 0);
    for (size_t j = 0; j < numScores; j++) {
      FLOAT_T scoreTarget = scores_[j].first;
      FLOAT_T scoreDecoy = scores_[j].second[i];
      bool targetBetter = scoreTarget != scoreDecoy ? scoreTarget < scoreDecoy : myrandom_limit(2) == 0;
      if (targetBetter) {
        optIdxCnt[j]++;
        histBin(ntds, scoreTarget);
      } else {
        histBin(ndds, scoreDecoy);
      }
    }

    int ntdsTotal = 0, nddsTotal = 0;
    for (size_t j = 0; j < numScores; j++) {
      sumNtds[j] += (ntdsTotal += ntds[j]);
      sumNdds[j] += (nddsTotal += ndds[j]);
    }
  }

  bool oldATDC = Params::GetBool("use-old-atdc");
  if (oldATDC) {
    for (size_t i = 0; i < numScores; i++) {
      sumNtds[i] /= numDecoySets;
      sumNdds[i] /= numDecoySets;
    }
  }
  
  vector<bool> isTargetPsm(numScores, true);
  int numTargetPsms = 0;
  vector< vector<int> > optIdxCntMat(numDecoySets + 1, vector<int>(numScores, 0));
  vector<size_t> currentOptIdx(numDecoySets + 1, 0);
  size_t lowestOptIdx = numDecoySets + 1;
  for (size_t i = 0; i < numScores; i++) {
    size_t idx = optIdxCnt[i];
    size_t idx2;
    if (oldATDC) {
      idx2 = currentOptIdx[idx]++;
    } else {
      idx2 = ++currentOptIdx[idx];
    }
    optIdxCntMat[idx][idx2] = i;
    if (idx < lowestOptIdx) {
      lowestOptIdx = idx;
    }
    int denominator;
    if (oldATDC) {
      denominator = 1;
    } else {
      denominator = numDecoySets;
    }
    if ((FLOAT_T)numTargetPsms <= (sumNtds[i] / denominator) - 0.5) {
      numTargetPsms++;
      continue;
    }
    if (oldATDC) {
      isTargetPsm[optIdxCntMat[lowestOptIdx][--currentOptIdx[lowestOptIdx]]] = false;
    } else {
      isTargetPsm[optIdxCntMat[lowestOptIdx][currentOptIdx[lowestOptIdx]--]] = false;
    }
    for ( ; lowestOptIdx < numDecoySets + 1 && currentOptIdx[lowestOptIdx] == 0; lowestOptIdx++);
  }

  if (oldATDC) {

    int targetPsmsTotal = 0;
    for (size_t i = 0; i < numScores; i++) {
      if (isTargetPsm[i]) {
        targetPsmsTotal++;
      }
      sumNdds[i] = (1 + sumNdds[i]) / max(1, targetPsmsTotal);
    }

  } else {
  
    if (bc1 >= 0) {
      int targetPsmsTotal = 0;
      for (size_t i = 0; i < numScores; i++) {
        if (isTargetPsm[i]) {
          targetPsmsTotal++;
        }
        sumNdds[i] = (bc1 + sumNdds[i] / numDecoySets) / max(1, targetPsmsTotal);
      }
    } else {
      vector<FLOAT_T> ndds(numScores, 0);
      FLOAT_T prev = 0;
      for (size_t i = 0; i < numScores; i++) {
        FLOAT_T current = sumNdds[i];
        ndds[i] = current - prev;
        prev = current;
      }
      int targetPsmsTotal = 0;
      for (size_t i = 0; i < numScores; i++) {
        if (isTargetPsm[i]) {
          targetPsmsTotal++;
        }
        sumNdds[i] = ((min((FLOAT_T)numDecoySets, max((FLOAT_T)1, ndds[i])) + sumNdds[i]) / numDecoySets) /
          max(1, targetPsmsTotal);
      }
    }
  }
  
  return sumNdds;
}

void AssignConfidenceApplication::AtdcScoreSet::histBin(vector<int>& hist, FLOAT_T x) const {
  // the value x is in the nth bin if edges[n] < x <= edges[n+1]
  // edges are [-inf, <target scores...>]
  for (size_t i = 0; i < scores_.size(); i++) {
    if (x <= scores_[i].first) {
      hist[i]++;
      return;
    }
  }
}

bool AssignConfidenceApplication::AtdcScoreSet::sortScoresAsc(
  const boost::tuple<FLOAT_T, int, int, int>& x,
  const boost::tuple<FLOAT_T, int, int, int>& y
) {
  FLOAT_T scoreX = x.get<0>(), scoreY = y.get<0>();
  return scoreX != scoreY ? scoreX < scoreY : sortScoresByIdentifier(x, y);
}

bool AssignConfidenceApplication::AtdcScoreSet::sortScoresDesc(
  const boost::tuple<FLOAT_T, int, int, int>& x,
  const boost::tuple<FLOAT_T, int, int, int>& y
) {
  FLOAT_T scoreX = x.get<0>(), scoreY = y.get<0>();
  return scoreX != scoreY ? scoreX > scoreY : sortScoresByIdentifier(x, y);
}

bool AssignConfidenceApplication::AtdcScoreSet::sortScoresByIdentifier(
  const boost::tuple<FLOAT_T, int, int, int>& x,
  const boost::tuple<FLOAT_T, int, int, int>& y
) {
  int fileX = x.get<1>(), fileY = y.get<1>();
  if (fileX != fileY) {
    return fileX < fileY;
  }
  int scanX = x.get<2>(), scanY = y.get<2>();
  if (scanX != scanY) {
    return scanX < scanY;
  }
  int chargeX = x.get<2>(), chargeY = y.get<2>();
  return chargeX < chargeY;
}

/**
 * \brief Compute q-values from a given set of scores, using a second
 * set of scores as an empirical null.  Sorts the incoming target
 * scores and returns a corresponding list of q-values.
 *
 * This function is only exported to allow unit testing.
 */
vector<FLOAT_T> AssignConfidenceApplication::compute_decoy_qvalues_tdc(
  vector<FLOAT_T>& target_scores,
  vector<FLOAT_T>& decoy_scores,
  bool ascending,
  FLOAT_T pi_zero
) {
  if (target_scores.empty() || decoy_scores.empty()) {
    carp(CARP_FATAL, "Cannot compute q-values (0 targets, 0 nulls).");
  }
  carp(CARP_DEBUG, "Computing decoy q-values with %d targets and %d decoys.",
       target_scores.size(), decoy_scores.size());

  // Sort both sets of scores.
  if (ascending) {
    sort(target_scores.begin(), target_scores.end(), Match::ScoreLess);
    sort(decoy_scores.begin(), decoy_scores.end(), Match::ScoreLess);
  } else {
    sort(target_scores.begin(), target_scores.end(), Match::ScoreGreater);
    sort(decoy_scores.begin(), decoy_scores.end(), Match::ScoreGreater);
  }

  // Compute false discovery rate for each target score.
  vector<FLOAT_T> qvalues;
  qvalues.reserve(target_scores.size());
  int decoy_idx = 0;
  for (int target_idx = 0; target_idx < target_scores.size(); target_idx++) {
    double target_score = target_scores[target_idx];

    // Find the index of the first decoy score greater than this target score.
    if (ascending) {
      while (decoy_idx < decoy_scores.size() &&
             Match::ScoreLess(decoy_scores[decoy_idx], target_score)) {
        carp(CARP_DEBUG, "Decoy score %g.", decoy_scores[decoy_idx]);
        decoy_idx++;
      }
    } else {   
      while (decoy_idx < decoy_scores.size() &&
             Match::ScoreGreater(decoy_scores[decoy_idx], target_score)) {
        carp(CARP_DEBUG, "Decoy score %g.", decoy_scores[decoy_idx]);
        decoy_idx++;
      }
    }

    // FDR = (#decoys + 1)/ #targets
    FLOAT_T fdr = ((FLOAT_T)(decoy_idx + 1)/(FLOAT_T)(target_idx + 1));
    if (fdr > 1.0) {
      fdr = 1.0;
    }
    carp(CARP_DEBUG, "FDR for score %g = min(1,%d/%d) = %g",
         target_score, decoy_idx, target_idx + 1, fdr);
    qvalues.push_back(fdr);
    // assign fdr for this and any following with same score
    while (target_idx < target_scores.size() - 1 && target_scores[target_idx + 1] == target_score) {
      target_idx++;
      qvalues.push_back(fdr);
    }
  }
  
  // Convert the FDRs into q-values.
  convert_fdr_to_qvalue(qvalues);

  return qvalues;
}

/**
 * \brief Compute q-values using mix-max procedure. This part is a
 * reimplementation of Uri Keich's code written in R.
 *
 */
vector<FLOAT_T> AssignConfidenceApplication::compute_decoy_qvalues_mixmax(
  vector<FLOAT_T>& target_scores,
  vector<FLOAT_T>& decoy_scores,
  bool ascending,
  FLOAT_T pi_zero
) {
  if (target_scores.empty() || decoy_scores.empty()) {
    carp(CARP_FATAL, "Cannot compute q-values (0 targets, 0 decoys).");
  }
  size_t num_targets = target_scores.size();
  size_t num_decoys = decoy_scores.size();
  if (target_scores.size() != decoy_scores.size()) {
    carp(CARP_WARNING, "The mix-max procedure is not well behaved when # targets (%d) != # of decoys (%d).",
         num_targets, num_decoys);
  }
  //estimate pi0 from data if it is not given.
  if (pi_zero == 1.0) {
    vector< pair<double, bool> > score_labels =
      ComputeQValues::getScoreVector(target_scores, decoy_scores, ascending);
    pi_zero = ComputeQValues::estimatePi0(score_labels);

    carp(CARP_INFO, "Estimated pi_zero = %f", pi_zero);
  }
  // continue with mix-max procedure
  for (int target_idx = 0; target_idx < num_targets; ++target_idx) {
    carp(CARP_DEBUG, "target_scores[%d]=%lf decoy_scores[%d]=%lf",
         target_idx, target_scores[target_idx],
         target_idx, decoy_scores[target_idx]);
  }

  //Sort decoy and target stores
  if (ascending) {
    sort(target_scores.begin(), target_scores.end(), greater<FLOAT_T>());
    sort(decoy_scores.begin(), decoy_scores.end(), greater<FLOAT_T>());
  } else {
    sort(target_scores.begin(), target_scores.end());
    sort(decoy_scores.begin(), decoy_scores.end());
  }

  //histogram of the target scores.
  vector<double> h_w_le_z(num_decoys + 1, 0); //histogram for N_{w<=z}
  vector<double> h_z_le_z(num_decoys + 1, 0); //histogram for N_{z<=z}

  int idx = 0;
  for (int i = 0; i < num_decoys; ++i) {
    while (idx < num_targets && (ascending ?
      decoy_scores[i] <= target_scores[idx] :
      decoy_scores[i] >= target_scores[idx])) {
      ++idx;
    }
    h_w_le_z[i] = (double)idx;
  }
  idx = 0;
  for (int i = 0; i < num_decoys; ++i) {
    while (idx < num_targets && (ascending ?
      decoy_scores[i] <= decoy_scores[idx] :
      decoy_scores[i] >= decoy_scores[idx])) {
      ++idx;
    }
    h_z_le_z[i] = (double)idx;
  }
  h_w_le_z[num_decoys] = (double)(num_targets);
  h_z_le_z[num_decoys] = (double)(num_decoys);

  vector<FLOAT_T> fdrmod(num_targets, 0);
  double estPx_lt_zj = 0.0;
  double E_f1_mod_run_tot = 0.0;
  int j = num_decoys-1;
  int k = num_targets-1;
  int n_z_ge_w = 0;
  int n_w_ge_w = 0;
  double qvalue;
  double cnt_z, cnt_w;
  double prev_fdr = -1;

  for (int i = num_targets - 1; i >= 0; --i) {
    while (j >= 0 && (ascending ? 
      decoy_scores[j] <= target_scores[i] : 
      decoy_scores[j] >= target_scores[i])) {
      cnt_w = h_w_le_z[j + 1];
      cnt_z = h_z_le_z[j + 1];
      estPx_lt_zj = (double)(cnt_w - pi_zero*cnt_z) / ((1.0 - pi_zero)*cnt_z);
      estPx_lt_zj = estPx_lt_zj > 1 ? 1 : estPx_lt_zj;
      estPx_lt_zj = estPx_lt_zj < 0 ? 0 : estPx_lt_zj;
      E_f1_mod_run_tot += estPx_lt_zj * ((1.0 - pi_zero));
      ++n_z_ge_w;
      --j;
    }
    while (k >= 0 && (ascending ?
      target_scores[k] <= target_scores[i] :
      target_scores[k] >= target_scores[i])) {
      ++n_w_ge_w;
      --k;
    }
    qvalue = ((double)n_z_ge_w * pi_zero + E_f1_mod_run_tot) / (double)(n_w_ge_w);
    fdrmod[i] = qvalue > 1.0 ? 1.0 : qvalue;
    
    //convert qvalues to fdr
    if (prev_fdr > fdrmod[i]) {
      fdrmod[i] = prev_fdr;
    }

    prev_fdr = fdrmod[i];
  }
  // Convert the FDRs into q-values.
  return fdrmod;
}

void AssignConfidenceApplication::peptide_level_filtering(
  MatchCollection* match_collection,
  map<string, FLOAT_T>* BestPeptideScore, 
  SCORER_TYPE_T score_type,
  bool ascending) {

    MatchIterator* temp_iter = new MatchIterator(match_collection);

    while (temp_iter->hasNext()) {
      Crux::Match* match = temp_iter->next();
      FLOAT_T score = match->getScore(score_type);
      string peptideStr = getPeptideSeq(match);

      FLOAT_T bestScore = 0.0;
      try {
        bestScore = BestPeptideScore->at(peptideStr);
      } catch (const out_of_range& oor) {
        BestPeptideScore->insert(make_pair(peptideStr, score));
        continue;
      }
      if ((ascending && bestScore > score) || (!ascending && score > bestScore)) {
        BestPeptideScore->at(peptideStr) = score;
      }
    }
    delete temp_iter;
}

string AssignConfidenceApplication::getPeptideSeq(Match* match) {
  string peptideSeq;
  if (Params::GetBool("combine-modified-peptides")) {
    peptideSeq = match->getPeptide()->getSequence();
  } else {
    peptideSeq = match->getPeptide()->getModifiedSequenceWithMasses();
  }
  if (Params::GetBool("combine-charge-states")) {
    peptideSeq += StringUtils::ToString(match->getCharge());
  }
  return peptideSeq;
}

map<pair<string, unsigned int>, bool>* AssignConfidenceApplication::getSpectrumFlag() {
  return spectrum_flag_;
}

void AssignConfidenceApplication::setSpectrumFlag(map<pair<string, unsigned int>, bool>* spectrum_flag) {
  spectrum_flag_ = spectrum_flag;
}

void AssignConfidenceApplication::setIterationCnt(unsigned int iteration_cnt) {
  iteration_cnt_ = iteration_cnt;
}

void AssignConfidenceApplication::setOutput(OutputFiles *output) {
  output_ = output;
}

void AssignConfidenceApplication::setIndexName(string index_name) {
  index_name_ = index_name;
}

void AssignConfidenceApplication::setFinalIteration(bool is_final) {
  is_final_ = is_final;
}

unsigned int AssignConfidenceApplication::getAcceptedPSMs() {
  return accepted_psms_;
}

int AssignConfidenceApplication::getDirection(SCORER_TYPE_T scoreType) {
  // Return -1 for higher score better, 1 for lower score better, or 0 for unknown
  switch (scoreType) {
    case SP:
    case XCORR:
    case LOGP_BONF_WEIBULL_XCORR: // negative log p-values
    case TIDE_SEARCH_REFACTORED_XCORR:
    case RESIDUE_EVIDENCE_SCORE:
    case PERCOLATOR_SCORE:
      // higher score better, ascending = false
      return -1;
    case EVALUE:
    case TIDE_SEARCH_EXACT_PVAL:
    case TIDE_SEARCH_EXACT_SMOOTHED:
    case RESIDUE_EVIDENCE_PVAL:
    case BOTH_PVALUE:
      // lower score better, ascending = true
      return 1;
    default:
      carp(CARP_FATAL, "Unknown score type %s", scorer_type_to_string(scoreType));
      // unknown score type
      return 0;
  }
}

/**
* \returns the command name for ComputeQValues
*/
string AssignConfidenceApplication::getName() const {
  return "assign-confidence";
}

/**
* \returns the description for ComputeQValues
*/
string AssignConfidenceApplication::getDescription() const {
  return
    "[[nohtml:Assign two types of statistical confidence measures (q-values "
    "and posterior error probabilities) to each PSM in a given set.]]"
    "[[html:<p>Given target and decoy scores, estimate a q-value for each "
    "target score. The q-value is analogous to a p-value but incorporates "
    "false discovery rate multiple testing correction. The q-value associated "
    "with a score threshold T is defined as the minimal false discovery rate "
    "(FDR) at which a score of T is deemed significant. In this setting, the "
    "q-value accounts for the fact that we are analyzing a large collection of "
    "scores. For confidence estimation afficionados, please note that this "
    "definition of \"q-value\" is independent of the notion of \"positive FDR\" "
    "as defined in (Storey <em>Annals of Statistics</em> 31:2013-2015:2003).</p>"
    "<p>To estimate FDRs, <code>assign-confidence</code> uses one of two "
    "different procedures. Both require that the input contain both target and "
    "decoy scores. The default, target-decoy competition (TDC) procedure is "
    "described in this article:</p><blockquote>Josh E. Elias and Steve P. Gygi. "
    "\"Target-decoy search strategy for increased confidence in large-scale "
    "protein identifications by mass spectrometry.\" <em>Nature Methods</em>. "
    "4(3):207-14, 2007.</blockquote><p>Note that <code>assign-confidence</code> "
    "implements a variant of the protocol proposed by Elias and Gygi: rather "
    "than reporting a list that contains both targets and decoys, <code>"
    "assign-confidence</code> reports only the targets. The FDR estimate is "
    "adjusted accordingly (by dividing by 2).</p><p>The alternative, <em>"
    "mix-max</em> procedure is described in this article:</p><blockquote>Uri "
    "Keich, Attila Kertesz-Farkas and William Stafford Noble. "
    "<a href=\"http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00081\">"
    "\"An improved false discovery rate estimation procedure for shotgun proteomics.\"</a> "
    "<i>Journal of Proteome Research</i>. 14(8):3148-3161, 2015.</blockquote>"
    "<p>Note that the mix-max procedure requires as "
    "input calibrated scores, such as Comet E-values or p-values produced "
    "using Tide-s <code>exact-p-value</code> option.</p>"
    "<p>The mix-max procedure requires that scores "
    "are reported from separate target and decoy searches. Thus, this approach "
    "is incompatible with a search that is run using the <code>--concat T"
    "</code> option to <code>tide-search</code> or the <code>--decoy_search 2"
    "</code> option to <code>comet</code>. On the other hand, the TDC "
    "procedure can take as input "
    "searches conducted in either mode (concatenated or separate). If given "
    "separate search results and asked to do TDC estimation, <code>"
    "assign-confidence</code> will carry out the target-decoy competition as "
    "part of the confidence estimation procedure.</p><p>In each case, the "
    "estimated FDRs are converted to q-values by sorting the scores then "
    "taking, for each score, the minimum of the current FDR and all of the FDRs "
    "below it in the ranked list.</p><p>A primer on multiple testing correction "
    "can be found here:</p><blockquote>William Stafford Noble. <a href=\""
    "http://www.nature.com/nbt/journal/v27/n12/full/nbt1209-1135.html\">\"How "
    "does multiple testing correction work?\"</a> <em>Nature Biotechnology</em>. "
    "27(12):1135-1137, 2009.</blockquote>]]";
}

/**
* \returns the command arguments
*/
vector<string> AssignConfidenceApplication::getArgs() const {
  string arr[] = {
    "target input+"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
* \returns the command options
*/
vector<string> AssignConfidenceApplication::getOptions() const {
  string arr[] = {
    "estimation-method",
    "decoy-prefix",
    "score",
    "sidak",
    "top-match-in",
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
    "list-of-files",
    "combine-charge-states",
    "combine-modified-peptides",
    "fileroot"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
* \returns the command outputs
*/
vector< pair<string, string> >  AssignConfidenceApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("assign-confidence.target.txt",
    "a <a href=\"../file-formats/txt-format.html\">tab-delimited text file</a> that contains the "
    "targets, sorted by score. The file will contain one new column, named "
    "\"&lt;method&gt; q-value\", where &lt;method&gt; is either \"tdc\" or \"mix-max\"."));
  outputs.push_back(make_pair("assign-confidence.log.txt",
    "a log file containing a copy of all messages that were printed to stderr."));
  outputs.push_back(make_pair("assign-confidence.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  return outputs;
}

/**
* \returns the filestem for ComputeQValues
*/
string AssignConfidenceApplication::getFileStem() const {
  return "assign-confidence";
}

COMMAND_T AssignConfidenceApplication::getCommand() const {
  return QVALUE_COMMAND; // FIXME: Rename this. --WSN 4 Feb 2016
}

/**
* \returns whether the application needs the output directory or not.
*/
bool AssignConfidenceApplication::needsOutputDirectory() const {
  return true;
}

void AssignConfidenceApplication::processParams() {
  if (Params::GetString("estimation-method") != "peptide-level") {
    if (Params::GetInt("top-match") != 1 && !Params::IsDefault("top-match")) {
      carp(CARP_WARNING, "assign-confidence can work with top-match = 1 only.");
    }
    Params::Set("top-match", 1);
  }
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

 
