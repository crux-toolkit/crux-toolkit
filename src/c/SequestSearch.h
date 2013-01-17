/**
 * \file SequestSearch.h
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Oct 2, 2009
 * PROJECT: crux
 * \brief The crux search routine that emulates SEQUEST.
 *
 * Scores all candidate peptides with Sp, deletes all but the 500
 * top-scoring candidates, scores remaining 500 with xcorr, sorts
 * results by xcorr and returns the top 5 plus the match with the best
 * Sp score.  Writes results to .sqt, .txt, and .pep.xml files.  Does not
 * compute p-values.
 *****************************************************************************/
#ifndef SEQUESTSEARCH_H
#define SEQUESTSEARCH_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"
#include "OutputFiles.h"


#include <string>
#include <vector>

class SequestSearch: public CruxApplication {

 protected:
  /* Private function definitions */

  /**
   * \brief Pring the target and decoy match collections to their
   * respective target and decoy files.
   *
   * Three possibilities: 1. combine the target and all decoy
   * collections and print to target file.  2. print targets to target
   * file and combine all decoys and print to one decoy file.  3. print
   * each collection to a separate file.
   * Possible side effectos: Collections may be merged and re-ranked.
   */
  void printMatches(
    OutputFiles& output_files,       ///< files to print to
    MatchCollection* target_psms, ///< target psms to print
    std::vector<MatchCollection*>& decoy_psms,///< decoy psms to print
    Crux::Spectrum* spectrum,            ///< all matches are to this spec
    bool combine_target_decoy,  ///< merge targets and decoys?
    int num_decoy_files              ///< merge decoys?
  );

 public:
  /**
   * \returns a blank SequestSearch object
   */
  SequestSearch();
  
  /**
   * Destructor
   */
  ~SequestSearch();

  /**
   * main method for SequestSearch
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for SequestSearch
   */
  virtual std::string getName();

  /**
   * \returns the description for SequestSearch
   */
  virtual std::string getDescription();

  /**
   * \returns the file stem of the application, default getName.
   */
  virtual std::string getFileStem();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();

 /**
  *hide SequestSearch 
  */
  virtual bool hidden(); 

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
