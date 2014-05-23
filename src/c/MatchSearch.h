/**
 * \file MatchSearch.h 
 * AUTHOR: Barbara Frewen
 * DATE: Aug 19, 2008
 * \brief Object for running search-for-matches
 *****************************************************************************/
#ifndef MATCHSEARCH_H
#define MATCHSEARCH_H

#include "CruxApplication.h"
#include "OutputFiles.h"


#include <string>

class MatchSearch : public CruxApplication {

 protected:
  /* Private function definitions */

  /**
   * \brief Look at matches and search parameters to determine if a
   * sufficient number PSMs have been found.  Returns TRUE if the
   * maximum number of modifications per peptide have been considered.
   * In the future, implement and option and test for a minimum score.
   * \returns TRUE if no more PSMs need be searched.
   */
  bool isSearchComplete(
    MatchCollection* matches, ///< matches to consider
    int mods_per_peptide); ///< number of mods per peptide searched so far

  /**
   * \brief Search the database OR index with up to num_peptide_mods from
   * the list for matches to the spectrum. 
   * Scored PSMs are added to the match_collection, possibly truncating
   * the collection and deleting existing matches in the collection.
   * After searching with each peptide mod, assess if there exists a
   * "good enough" match and end the search if there is, returning the
   * number of peptide mods that were searched.
   * \return The number of peptide mods searched.
   */
  int searchPepMods(
    MatchCollection* match_collection, ///< store PSMs here
    bool is_decoy,   ///< generate decoy peptides from index/db
    Index* index,       ///< index to use for generating peptides
    Database* database, ///< db to use for generating peptides
    Crux::Spectrum* spectrum, ///< spectrum to search
    SpectrumZState& zstate, ///< seach spectrum at this z-state
    PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to apply
    int num_peptide_mods, ///< how many p_mods to use from the list
    bool compute_sp,  ///< compute sp scores
    bool store_scores///< save all scores for p-value estimation
  );

  /**
   * Print the target and decoy match collections to their respective
   * target and decoy files.
   *
   * Three possibilities: 1. combine the target and all decoy
   * collections and print to target file.  2. print targets to target
   * file and combine all decoys and print to one decoy file.  3. print
   * each collection to a separate file.
   * Possible side effectos: Collections may be merged and re-ranked.
   */
  void printSpectrumMatches(
    OutputFiles& output_files, ///< files to print to     
    MatchCollection* target_psms, ///< target psms to print
    vector<MatchCollection*>& decoy_psms, ///< decoy psms to print
    Crux::Spectrum* spectrum, ///< spectrum for all psms
    bool combine_target_decoy, ///< print target and decoys to same file
    int num_decoy_files ///< number of decoy files to print to
   );

  // TODO this should be in match_collection
  /**
   * Search the given database or index using shuffled peptides and the
   * spectrum/charge in the target psm match collection.  Add those
   * scores to the target psm match collection for use in weibull
   * parameter estimation but do not save the matches.  Repeat the
   * search with all peptide mods in the list.
   */
  void addDecoyScores(
    MatchCollection* target_psms, ///< add scores to these matches
    Crux::Spectrum* spectrum, ///< spectrum to score
    SpectrumZState& zstate, ///< charge and mass to use for spectrum
    Index* index, ///< search this index if not null
    Database* database, ///< search this database if not null
    PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to search
    int num_peptide_mods ///< number of mods in the above array
  );

  

 public:

  /**
   * \returns a blank MatchSearch object
   */
  MatchSearch();

  /**
   * Destructor
   */
  virtual ~MatchSearch();

  /**
   * main method for MatchSearch
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for MatchSearch
   */
  virtual std::string getName();

  /**
   * \returns the file stem of the application, default getName.
   */ 
  virtual std::string getFileStem();

  /**
   * \returns the description for MatchSearch
   */
  virtual std::string getDescription();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  /**
   * \returns whether the application needs the output directory or
   * not. (default false). 
   */
  virtual bool needsOutputDirectory();

  virtual bool hidden();

};



#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
