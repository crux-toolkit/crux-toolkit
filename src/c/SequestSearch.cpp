/**
 * \file SequestSearch.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Oct 2, 2009
 * PROJECT: crux
 * \brief The crux search routine that emulates SEQUEST.
 *
 * Scores all candidate peptides with Sp, deletes all but the 500
 * top-scoring candidates, scores remaining 500 with xcorr, sorts
 * results by xcorr and returns the top 5 plus the match with the best
 * Sp score.  Writes results to .sqt, .txt, and .csm files.  Does not
 * compute p-values.
 *****************************************************************************/
#include "SequestSearch.h"

#include "FilteredSpectrumChargeIterator.h"
#include "SearchProgress.h"
#include "SpectrumCollectionFactory.h"
#include "ModifiedPeptidesIterator.h"

using namespace std;
using namespace Crux;

/**
 * \returns a blank SequestSearch object
 */
SequestSearch::SequestSearch() {

}

/**
 * Destructor
 */
SequestSearch::~SequestSearch() {
}

/**
 * \brief The starting point for the crux sequest-search command.
 *
 * After parsing command line and opening input and output files,
 * iterates over every spectrum at each charge state.  Generates
 * matches for each and prints.  Matches to decoy spectra can also be
 * generated.  If num-decoys-per-target is n, there are n decoy match
 * collections which may be merged for printing.  
 * If the user gave the command line
 *
 * crux sequest-search [options] ms2file proteins
 * 
 * then this is passed everything after the 'crux' token
 */
int SequestSearch::main(int argc,   ///< number of cmd line tokens
                        char** argv)///< array of cmd line tokens
{

  const char* option_list[] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "mzid-output",
    "pin-output",
    "pepxml-output",
    "txt-output",
    "spectrum-parser",
    "spectrum-min-mz",
    "spectrum-max-mz",
    "spectrum-charge",
    "output-dir",
    "scan-number",
    "fileroot",
    "decoys",
    "num-decoys-per-target",
    "decoy-location"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  // Define required command line arguments 
  const char* argument_list[] = {"ms2 file", "protein-database"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize(argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  // Get input: protein file
  char* input_file = get_string_parameter("protein-database");

  // Prepare input, fasta or index 
  Index* index = NULL;
  Database* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database); 
  free(input_file);

  carp(CARP_DEBUG, "Found %i proteins", num_proteins);
  if( num_proteins == 0 ){
    carp(CARP_FATAL, "No proteins were found in the protein source.");
  }
  
  // Get input: ms2 file
  const char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  Crux::SpectrumCollection* spectra = SpectrumCollectionFactory::create(ms2_file);

  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!spectra->parse()){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       spectra->getNumSpectra());

  // Prepare output files 
  
  bool combine_target_decoy = get_boolean_parameter("tdc");
  OutputFiles output_files(this); 
  output_files.writeHeaders(num_proteins, combine_target_decoy);

  // get search parameters for match_collection
  int num_decoy_files = get_int_parameter("num-decoy-files");
  bool have_index = (index != NULL);
  int num_decoys_per_target = get_num_decoys(have_index); 

  SearchProgress progress;

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  // create spectrum iterator
  FilteredSpectrumChargeIterator* spectrum_iterator = 
    new FilteredSpectrumChargeIterator(spectra);

  // Perform search on each spectrum
  while(spectrum_iterator->hasNext()){

    SpectrumZState zstate;
    Spectrum* spectrum = 
      spectrum_iterator->next(zstate);

    double mz = spectrum->getPrecursorMz();

    progress.report(spectrum->getFirstScan(), zstate.getCharge());

    // create empty match collections to store results in
    MatchCollection* target_psms = new MatchCollection(false); 
    target_psms->setZState(zstate);

    vector<MatchCollection*> decoy_psm_collections;
    for(int decoy_idx=0; decoy_idx < num_decoys_per_target; decoy_idx++){
      MatchCollection* psms = new MatchCollection(true);
      psms->setZState(zstate);
      decoy_psm_collections.push_back(psms);
    }

    // search with one peptide iterator for each peptide mod
    for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){

      // get peptide mod
      PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];

      // get peptide iterator
     ModifiedPeptidesIterator* peptide_iterator =
        new ModifiedPeptidesIterator(mz,
                                     zstate,
                                     peptide_mod, 
                                     false, // not decoy
                                     index,
                                     database);

      // add matches to targets
      int added = target_psms->addMatches(spectrum,
                              zstate,
                              peptide_iterator,
                              false, // not decoy
                              false, // don't save scores for p-values
                              true,  // do preliminary Sp scoring
                              true   // filter by Sp
                              ); 

      // add matches to each decoy
      for(int decoy_idx = 0; decoy_idx < num_decoys_per_target; decoy_idx++){

        // get new peptide iterator
        delete peptide_iterator;
        peptide_iterator = new ModifiedPeptidesIterator(mz,
                                                        zstate,
                                                        peptide_mod, 
                                                        true,  // is decoy
                                                        index,
                                                        database);
        // add matches
        MatchCollection* cur_decoys = decoy_psm_collections.at(decoy_idx);
        cur_decoys->setTargetExperimentSize(target_psms->getExperimentSize());
        cur_decoys->addMatches(spectrum,
                    zstate,
                    peptide_iterator,
                    true,  // is decoy
                    false, // don't save scores for p-values
                    true,  // do preliminary Sp scoring
                    true   // filter by Sp
                    ); 
      }

      carp(CARP_DEBUG, "Found %d peptides.", added);

      // clean up for next peptide mod
      delete peptide_iterator;

    } // next peptide mod

    // print matches
    int total_matches = target_psms->getMatchTotal();

    if( total_matches == 0 ){
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i.",
           spectrum->getFirstScan(), zstate.getCharge());
      progress.increment(false);

    }else{  
      printMatches(output_files, target_psms, decoy_psm_collections,
                    spectrum, combine_target_decoy, num_decoy_files);
      progress.increment(true);
    }

    // clean up
    delete target_psms;
    vector<MatchCollection*>::iterator it = decoy_psm_collections.begin();
    for(; it < decoy_psm_collections.end(); ++it){
      MatchCollection* psms = *it;
      delete psms;
    }

  } // next spectrum

  output_files.writeFooters();

  // clean up
  delete spectrum_iterator;
  delete spectra;
  for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    free_peptide_mod(peptide_mods[mod_idx]);
  }
  free(peptide_mods);
  Index::free(index);
  Database::freeDatabase(database);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux sequest-search");
  
  return 0;
}

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
void SequestSearch::printMatches(
  OutputFiles& output_files,       ///< files to print to
  MatchCollection* target_psms, ///< target psms to print
  vector<MatchCollection*>& decoy_psms,///< decoy psms to print
  Spectrum* spectrum,            ///< all matches are to this spec
  bool combine_target_decoy,  ///< merge targets and decoys?
  int num_decoy_files              ///< merge decoys?
){ 

  if( combine_target_decoy ){

    // merge all collections
    MatchCollection* all_psms = target_psms;
    for(size_t decoy_idx = 0; decoy_idx < decoy_psms.size(); decoy_idx++){
      MatchCollection::merge(decoy_psms.at(decoy_idx), all_psms);
    }

    // sort and rank
    all_psms->populateMatchRank(SP);
    all_psms->saveTopSpMatch();
    all_psms->populateMatchRank(XCORR);
    vector<MatchCollection*> empty_list;
    output_files.writeMatches(all_psms, // target matches
                              empty_list,     // decoy matches
                              XCORR,    // use XCORR rank for cutoff
                              spectrum); 


  }else{ // targets and decoys in separate files

    if( num_decoy_files == 1 ){ // combine decoys

      // merge decoys
      MatchCollection* merged_decoy_psms = decoy_psms.at(0);
      for(size_t decoy_idx = 1; decoy_idx < decoy_psms.size(); decoy_idx++){
        MatchCollection::merge(decoy_psms.at(decoy_idx),
                                merged_decoy_psms);
      }
      
      // re-sort and rank if we merged multiple collections
      if( decoy_psms.size() > 1 ){
        merged_decoy_psms->populateMatchRank(SP);
        merged_decoy_psms->saveTopSpMatch();
        merged_decoy_psms->populateMatchRank(XCORR);
      }
      vector<MatchCollection*> decoy_list(1, merged_decoy_psms);
      output_files.writeMatches(target_psms, decoy_list, 
                                XCORR, spectrum);
      
    }else{   // write targets and decoys to separate files
      // TODO write a version of OutputFiles::writeMatches that takes
      // a vector of decoy match collections
      int num_decoys = decoy_psms.size();
      vector<MatchCollection*> decoy_psm_array;
      for(int decoy_idx = 0; decoy_idx < num_decoys; decoy_idx++){
        decoy_psm_array.push_back(decoy_psms.at(decoy_idx));
      }

      output_files.writeMatches(target_psms, decoy_psm_array, 
                                XCORR, spectrum);
    }
  }


}



/**
 * \returns the command name for SequestSearch
 */
string SequestSearch::getName() {
  return "sequest-search";
}

/**
 * \returns the description for SequestSearch
 */
string SequestSearch::getDescription() {
  return 
  "Similar to search-for-matches but use Sp as a "
  "preliminary score followed by XCorr.";

}

/**
 * \returns the file stem of the application, default getName.
 */
string SequestSearch::getFileStem() {
  return "sequest";
}

/**
 * \returns the enum of the application, default MISC_COMMAND
 */
COMMAND_T SequestSearch::getCommand() {
  return SEQUEST_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool SequestSearch::needsOutputDirectory() {
  return true;
}

/**
 * hide sequest search 
*/

bool SequestSearch:: hidden(){
  return true; 
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
