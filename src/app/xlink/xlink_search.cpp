
/**
 * \file xlink_search.cpp
 * \brief Object for running search-for-xlinks (new code)
 *****************************************************************************/
#include "SearchForXLinks.h"
#include "XLinkMatch.h"
#include "XLinkMatchCollection.h"
#include "XLinkBondMap.h"
#include "XLinkPeptide.h"
#include "XLinkIonSeriesCache.h"
#include "xlink_compute_qvalues.h"

#include "Weibull.h"

//CRUX INCLUDES
#include "model/objects.h"
#include "model/FilteredSpectrumChargeIterator.h"
#include "io/OutputFiles.h"
#include "io/SpectrumCollectionFactory.h"
#include "util/Params.h"
#include "XLinkDatabase.h"


//C++ Includes
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <ctime>



using namespace std;

void buildArguments(
  vector<string>& args_vec, 
  int &argc, 
  char** &argv) {

  argc = args_vec.size();
  argv = new char*[argc];

  argv[0] = (char*)args_vec[0].c_str();
  carp(CARP_INFO, "argv[0]=%s", argv[0]);
  for (int idx = 1; idx < argc; idx++) {
    argv[idx] = (char*)args_vec[idx].c_str();
    carp(CARP_INFO, "argv[%i]=%s", idx, argv[idx]);
  }

}

void writeTrainingCandidates(XLinkMatchCollection* training_candidates, int scan_num, Weibull& weibull) {
  training_candidates->populateMatchRank(XCORR);
  training_candidates->sort(XCORR);
  string output_dir = Params::GetString("output-dir");

  string output_file = output_dir + "/" +
    "scan." + StringUtils::ToString(scan_num) + 
    ".charge." + StringUtils::ToString(training_candidates->getCharge()) +
    ".training.candidates.txt";
  //cerr<<"writing "<<output_file<<endl;
  DelimitedFileWriter writer(output_file.c_str());

  
  
  writer.setColumnName("scan", 0);
  writer.setColumnName("charge", 1);
  writer.setColumnName("decoy", 2);
  writer.setColumnName("sequence", 3);
  writer.setColumnName("eta", 4);
  writer.setColumnName("beta", 5);
  writer.setColumnName("shift", 6);
  writer.setColumnName("correlation", 7);
  writer.setColumnName("xcorr score", 8);
  writer.setColumnName("p-value", 9);
  writer.setColumnName("p-value ecdf", 10);
  writer.writeHeader();

  for (size_t idx = 0 ;idx < training_candidates->getMatchTotal();idx++) {
    writer.setColumnCurrentRow(0, scan_num);
    writer.setColumnCurrentRow(2, (*training_candidates)[idx]->isDecoy());
    writer.setColumnCurrentRow(1, (*training_candidates)[idx]->getCharge());
    writer.setColumnCurrentRow(3, (*training_candidates)[idx]->getSequenceString());
    writer.setColumnCurrentRow(4, weibull.getEta());
    writer.setColumnCurrentRow(5, weibull.getBeta());
    writer.setColumnCurrentRow(6, weibull.getShift());
    writer.setColumnCurrentRow(7, weibull.getCorrelation());
    writer.setColumnCurrentRow(8, (*training_candidates)[idx]->getScore(XCORR));
    FLOAT_T score = (*training_candidates)[idx]->getScore(XCORR);
    writer.setColumnCurrentRow(9, weibull.getWeibullPValue(score));
    writer.setColumnCurrentRow(10, weibull.getECDFPValue(score));
    writer.writeRow();
  }

}


/**
 * main method for SearchForXLinks that implements the refactored code
 */
int SearchForXLinks::xlinkSearchMain() {
  
  carp(CARP_INFO, "Beginning crux search-for-xlinks.");

  /* Get parameters */
  carp(CARP_DETAILED_INFO, "Getting parameters.");

  vector<string> ms2_files = Params::GetStrings("ms2 file");

  string input_file = Params::GetString("protein fasta file");
  string output_directory = Params::GetString("output-dir");
  int top_match = Params::GetInt("top-match");
  XLinkPeptide::setLinkerMass(Params::GetDouble("link mass"));
  int min_weibull_points = Params::GetInt("min-weibull-points");
  bool compute_pvalues = Params::GetBool("compute-p-values");

  /* Prepare input fasta  */
  carp(CARP_INFO, "Preparing database.");
  XLinkDatabase::initialize();
  Database* database = NULL;
  int num_proteins = 0;//prepare_protein_input(input_file, &database);
  carp(CARP_DETAILED_INFO, "Number of proteins: %d",num_proteins);
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = 0;

  /* Usually for debugging purposes, print out the database of canddiates */
  if (Params::GetBool("xlink-print-db")) {
    carp(CARP_INFO, "Generating and printing xlink database.");
    XLinkDatabase::print();
  }

  SpectrumZState zstate;

  /* Prepare output files */
  carp(CARP_DETAILED_INFO, "Preparing output files.");
  OutputFiles output_files(this);
  output_files.writeHeaders(num_proteins);

  for (vector<string>::iterator ms2_file_iter = ms2_files.begin();
       ms2_file_iter != ms2_files.end(); ++ms2_file_iter) { 
    string ms2_file = *ms2_file_iter;
    
    carp(CARP_INFO, "Loading spectra %s.", ms2_file.c_str());
    Crux::Spectrum* spectrum = NULL;
    Crux::SpectrumCollection* spectra =
      SpectrumCollectionFactory::create(ms2_file);
    spectra->parse();

    FilteredSpectrumChargeIterator* spectrum_iterator =
      new FilteredSpectrumChargeIterator(spectra);


    // main loop over spectra in ms2 file
    int scan_num = 0;
    int skipped_no_candidates = 0;
    int search_count = 0;
    FLOAT_T num_spectra = (FLOAT_T)spectra->getNumSpectra();
    FLOAT_T min_pvalue = 1.0 / num_spectra;
    //class for estimating pvalues.
    Weibull weibull;
  
    // for every observed spectrum 
    carp(CARP_INFO, "Beginning search.");
    int print_interval = Params::GetInt("print-search-progress");

    while (spectrum_iterator->hasNext()) {

      spectrum = spectrum_iterator->next(zstate);
      scan_num = spectrum->getFirstScan();
      
      if (print_interval > 0 && search_count > 0 && search_count % print_interval == 0) {
	carp(CARP_INFO, 
	     "%d spectrum-charge combinations searched, %.0f%% complete",
	     search_count + spectrum_iterator->numSkipped(),
	     (search_count + spectrum_iterator->numSkipped()) / num_spectra * 100);
      }
      search_count++;

      FLOAT_T precursor_mz = spectrum->getPrecursorMz();

      XLinkMatchCollection* target_candidates =
	new XLinkMatchCollection(
				 spectrum,
				 zstate,
				 false,
				 false
				 );

      if (target_candidates->getMatchTotal() < 0) {
	carp(CARP_ERROR, "Scan %d has %d candidates.", scan_num, 
	     target_candidates->getMatchTotal());
      } else if (target_candidates->getMatchTotal() == 0) {
	skipped_no_candidates++;
	carp(CARP_DETAILED_INFO, "Skipping scan %d charge %d mass %lg", 
	     scan_num, 
	     zstate.getCharge(),
	     zstate.getNeutralMass()
	     );
	delete target_candidates;
	//XLink::deleteAllocatedPeptides();
	continue;
      }

      carp(CARP_DETAILED_INFO, "Scan=%d charge=%d mass=%lg candidates=%d", 
	   scan_num, 
	   zstate.getCharge(), 
	   zstate.getNeutralMass(), 
	   target_candidates->getMatchTotal());   

      // Score targets.
      target_candidates->scoreSpectrum(spectrum);

      // Score decoys.
      carp(CARP_DEBUG, "Getting decoy candidates.");
      XLinkMatchCollection* decoy_candidates = new XLinkMatchCollection();
      target_candidates->shuffle(*decoy_candidates);
      carp(CARP_DEBUG, "Scoring decoys.");
      decoy_candidates->scoreSpectrum(spectrum);
      
      if (compute_pvalues) {
	weibull.reset();
	XLinkMatchCollection *target_train_candidates =
	  new XLinkMatchCollection(
				   spectrum,
				   zstate,
				   false,
				   true);
	XLinkMatchCollection *train_candidates =
	  new XLinkMatchCollection(
				   spectrum,
				   zstate,
				   true,
				   true
				   );
      
	for (size_t idx=0;idx < target_train_candidates->getMatchTotal();idx++) {
	  train_candidates->add(target_train_candidates->at(idx), true);
	}
	while(train_candidates->getMatchTotal() < min_weibull_points) {
	  target_train_candidates->shuffle(*train_candidates);
	}
	train_candidates->scoreSpectrum(spectrum);
	for (int idx = 0;idx < train_candidates->getMatchTotal();idx++) {
	  const string& sequence = (*train_candidates)[idx]->getSequenceStringConst();
	  FLOAT_T score = (*train_candidates)[idx]->getScore(XCORR);
	  weibull.addPoint(sequence, score);
	}
	bool write_weibull_points = !weibull.fit();
	
	target_candidates->sort(XCORR);
	
	
	// Calculate pvalues.
	int nprint = min(top_match,target_candidates->getMatchTotal());
	carp(CARP_DEBUG, "Calculating %d target p-values.", nprint);
	for (int idx=0;idx < nprint;idx++) {
	  FLOAT_T score = (*target_candidates)[idx]->getScore(XCORR);
	  (*target_candidates)[idx]->setPValue(weibull.getPValue(score));
	}
	
	nprint = min(top_match, (int)decoy_candidates->getMatchTotal());
	carp(CARP_DEBUG, "Calculating %d decoy p-values.", nprint);
	decoy_candidates->sort(XCORR);
	for (int idx=0;idx < nprint;idx++) {
	  FLOAT_T score = (*decoy_candidates)[idx]->getScore(XCORR);
	  (*decoy_candidates)[idx]->setPValue(weibull.getPValue(score));
	  FLOAT_T wpvalue = weibull.getWeibullPValue(score);
	  FLOAT_T bpvalue = bonferroni_correction(wpvalue, decoy_candidates->getMatchTotal()) * 2.0;
	  if ((wpvalue == 0) || (wpvalue != wpvalue) || (bpvalue  < min_pvalue)) {
	    //If we have a bad fit, 0 or too low pvalue, print out the points.
	    write_weibull_points = true;
	  }
        
	}
      
      
	if (write_weibull_points || Params::GetBool("write-weibull-points")) {
	  writeTrainingCandidates(train_candidates, scan_num, weibull);
	}
	carp(CARP_DEBUG, "Delete train candidates.");
	delete train_candidates;
	carp(CARP_DEBUG, "Delete target train candidates.");
	delete target_train_candidates;
	
      } // if (compute_p_values)
      
      //print out
      target_candidates->setFilePath(ms2_file);
      decoy_candidates->setFilePath(ms2_file);
      vector<MatchCollection*> decoy_vec;
    
      if (Params::GetBool("concat")) {
	for (size_t idx=0;idx < decoy_candidates->getMatchTotal();idx++) {
	  target_candidates->add(decoy_candidates->at(idx), true);
	}
      } else {
	
	decoy_vec.push_back(decoy_candidates);

	if (decoy_candidates->getScoredType(SP) == true) {
	  decoy_candidates->populateMatchRank(SP);
	}
	decoy_candidates->populateMatchRank(XCORR);
	decoy_candidates->sort(XCORR);
      }
      
      carp(CARP_DEBUG, "Ranking.");
      
      if (target_candidates->getScoredType(SP) == true) {
	target_candidates->populateMatchRank(SP);
      }
      target_candidates->populateMatchRank(XCORR);
      target_candidates->sort(XCORR);
      
    
      carp(CARP_DEBUG, "Writing results.");
      output_files.writeMatches(
				(MatchCollection*)target_candidates, 
				decoy_vec,
				XCORR,
				spectrum);

      /* Clean up */
      carp(CARP_DEBUG, "Deleting decoy candidates.");
      delete decoy_candidates;
      carp(CARP_DEBUG, "Deleting target candidates.");
      delete target_candidates;
      XLink::deleteAllocatedPeptides();
    
      //free_spectrum(spectrum);
      
      carp(CARP_DEBUG, "Done with spectrum %d.", scan_num);
      carp(CARP_DEBUG, "=====================================");
    } // get next spectrum

    carp(CARP_INFO, "Skipped %d (%g%%) spectra with 0 candidates.", 
	 skipped_no_candidates, skipped_no_candidates / num_spectra * 100);
    
    carp(CARP_INFO, "Skipped %d (%g%%) spectra with too few peaks.", 
	 spectrum_iterator->numSkipped(), 
	 spectrum_iterator->numSkipped() / num_spectra * 100);

    output_files.writeFooters();

    // clean up


    delete spectrum_iterator;
    delete spectra;
    XLink::deleteAllocatedPeptides();
  }
  for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++) {
    free_peptide_mod(peptide_mods[mod_idx]);
  }
  free(peptide_mods);

  //finalize_weibull();
  //Scorer::finalize();
  XLinkIonSeriesCache::finalize();
  IonSeries::finalize();
  XLinkDatabase::finalize();
  //modifications_finalize(); TODO - Figure where to free the modification cache.

  //Calculate q-values via p-values from weibull fit.
  if (compute_pvalues) {
    carp(CARP_DEBUG, "Computing q-values using p-values.");
    if (Params::GetBool("concat")) {
      carp(CARP_WARNING, "Computing Q-Values with concatenated target-decoy files not implemented yet!");
    } else {
      xlink_compute_qvalues();
    }
  }

  return(0);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
