//TODO - change cout/cerrs to carp.

#include "xhhc.h"
#include "LinkedIonSeries.h"
#include "xhhc_scorer.h"
#include "LinkedPeptide.h"

#include "crux-utils.h"
#include "objects.h"
#include "Scorer.h"
#include "DelimitedFile.h"

#include "xlink_compute_qvalues.h"

#include <math.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <iostream>
#include <fstream>

void getBestBonf(DelimitedFile& matches, int start, int stop, 
  int& best_index, double& best_bonf) {

  int numScans = stop - start + 1;
  //if there is only 1 scan, then return it.
  if (numScans == 1) {
    best_index = start;
    double pvalue = matches.getDouble("p-value", best_index);
    int ntests = matches.getInteger("distinct matches/spectrum", best_index);
    best_bonf = bonferroni_correction(pvalue, ntests);
  } else {
    map<int, pair<int, double> > charge_best_score; 
    map<int, int> charge_ntests;

    for (int match_idx = start;match_idx <= stop;match_idx++) {
      int charge = matches.getInteger("charge", match_idx);
      double pvalue = matches.getDouble("p-value", match_idx);

      if (charge_best_score.find(charge) == charge_best_score.end()) {
        charge_best_score[charge] = make_pair(match_idx, pvalue);
      } else if (pvalue < charge_best_score[charge].second) {
        charge_best_score[charge] = make_pair(match_idx, pvalue);
      }
      if (charge_ntests.find(charge) == charge_ntests.end()) {
        int ntests = matches.getInteger("distinct matches/spectrum", match_idx);
        charge_ntests[charge] = ntests;
      }
    }
  
  
    int best_charge = charge_best_score.begin() -> first;
    double best_charge_pvalue = charge_best_score[best_charge].second;
    double best_charge_pvalue_bonf = 
      bonferroni_correction(best_charge_pvalue, charge_ntests[best_charge]);

    int best_charge_idx = charge_best_score[best_charge].first;

    for (map<int, pair<int, double> >::iterator iter = 
      charge_best_score.begin();
      iter != charge_best_score.end();
      ++iter) {

      double current_charge_pvalue_bonf =
        bonferroni_correction(iter -> second.second, charge_ntests[iter -> first]);
      if (current_charge_pvalue_bonf < best_charge_pvalue_bonf) {
        best_charge = iter -> first;
        best_charge_pvalue = iter -> second.second;
        best_charge_pvalue_bonf = current_charge_pvalue_bonf;
        best_charge_idx = charge_best_score[best_charge].first;  
      }
    
    }

    int ntests_total = 0;
    for (map<int, int>::iterator iter =
      charge_ntests.begin();
      iter != charge_ntests.end();
      ++iter) {

      ntests_total += iter -> second;
    }

    best_bonf = bonferroni_correction(best_charge_pvalue, ntests_total);
    best_index = best_charge_idx;
  }
}

void collapseScans(DelimitedFile& matches_in, DelimitedFile& matches_out) {

  matches_out.clear();
  //matches_out.addColumns(matches_in.getColumnNames());
  vector<string>& column_names = matches_in.getColumnNames();
  for (unsigned int idx=0;idx<column_names.size();idx++) {
    matches_out.addColumn(column_names[idx]);
  }
  
  //make sure scans are together.
  carp(CARP_DETAILED_DEBUG, "Sorting matches by scan");
  matches_in.sortByIntegerColumn("scan");

  matches_out.addColumn("p-value bonf.");

  if (matches_in.numRows() == 0) {
    return;
  }

  int last_scan = matches_in.getInteger("scan", 0);
  int first_row = 0;
  int best_row = 0;
  double best_bonf = 0;
  

  for (unsigned int match_idx = 0;
    match_idx < matches_in.numRows(); 
    match_idx++) {

    int current_scan = matches_in.getInteger("scan", match_idx);
    if (last_scan != current_scan) {
      //process the scans between the first and match_idx-1.
      //find the best row and calculate the bonferroni corrected p-value.
      carp(CARP_DEBUG,"Collaping %d %d %d",last_scan, first_row, match_idx-1);
      getBestBonf(matches_in, first_row, match_idx-1, best_row, best_bonf);  
      carp(CARP_DEBUG,"Copying row best: %d %lf", best_row, best_bonf);
      //update matches out
      int new_row = matches_out.addRow();
      matches_in.copyToRow(matches_out, best_row, new_row);
      carp(CARP_DEBUG,"Adding bonferroni");
      matches_out.setValue<double>("p-value bonf.", new_row, best_bonf);


      //update first_row and last scan.
      first_row = match_idx;
      last_scan = current_scan;
    }
  }

  carp(CARP_DEBUG,"Doing last entry");
  //finish the last entry
  getBestBonf(matches_in, first_row, matches_in.numRows() - 1, best_row, best_bonf);
  carp(CARP_DEBUG,"best: %d %lf", best_row, best_bonf);
  int new_row = matches_out.addRow();
  matches_in.copyToRow(matches_out, best_row, new_row);
  matches_out.setValue("p-value bonf.", new_row, best_bonf);




}


int xlink_compute_qvalues(){

  /* Get Arguments */

  carp(CARP_INFO,"reading targets");
  string output_dir = get_string_parameter_pointer("output-dir");
  string target_filename = "search-for-xlinks.target.txt";

  //Read in targets.
  string target_path = output_dir + "/" + target_filename;
  DelimitedFile target_matches(target_path);

  //cout << target_matches << endl;

  //Read in decoys.

  carp(CARP_INFO,"reading decoys");

  string decoy_filename = "search-for-xlinks.decoy.txt";
  string decoy_path = output_dir + "/" + decoy_filename;
  DelimitedFile decoy_matches(decoy_path);

  //cout << decoy_matches << endl;

  //Collapse to the best target/decoy and calculate bonferonni corrected p-value.

  carp(CARP_INFO,"collapsing targets");
  DelimitedFile target_matches_bonf;
  collapseScans(target_matches, target_matches_bonf);

  carp(CARP_INFO,"collapsing decoys");
  DelimitedFile decoy_matches_bonf;
  collapseScans(decoy_matches, decoy_matches_bonf);

  //Sort by increasing p-value.
  carp(CARP_INFO,"Sorting by p-value bonf.");
  target_matches_bonf.sortByFloatColumn("p-value bonf.");
  decoy_matches_bonf.sortByFloatColumn("p-value bonf."); 

  //test the decoy p-values for accuracy.
  for (unsigned int idx=0;idx < decoy_matches_bonf.numRows();idx++) {
    double calc_pvalue = decoy_matches_bonf.getDouble("p-value bonf.", idx);
    double rank_pvalue = (double)(idx + 1) / (double) decoy_matches_bonf.numRows();

    if ((calc_pvalue > 2 * rank_pvalue) || (calc_pvalue < 0.5 * rank_pvalue)) {
      carp(CARP_WARNING, "inaccurate p-values!");
      carp(CARP_WARNING, "scan:%d charge:%d mass: %g rank:%g calc:%g", 
        decoy_matches_bonf.getInteger("scan", idx),
        decoy_matches_bonf.getInteger("charge", idx),
        decoy_matches_bonf.getDouble("spectrum neutral mass"),
        rank_pvalue,
        calc_pvalue);
    } 
  }
  
  carp(CARP_INFO,"Adding q-value column");
  target_matches_bonf.addColumn("fdr b-h");
  target_matches_bonf.addColumn("fdr decoy");
  target_matches_bonf.addColumn("q-value b-h");
  target_matches_bonf.addColumn("q-value decoy");

  //carp(CARP_INFO,"Calculating q-values");
  unsigned int decoy_idx = 0;
  //cout <<"Number of target rows:"<<target_matches_bonf.numRows()<<endl;
  for (unsigned int target_idx = 0; 
    target_idx < target_matches_bonf.numRows(); 
    target_idx++) {

    //calculate q-value by b-h
    double current_pvalue = target_matches_bonf.getDouble("p-value bonf.", target_idx);
    //cout <<"Current pvalue:"<<current_pvalue<<endl;
    
    double fdr_bh = (double) target_matches_bonf.numRows() /
      (double)(target_idx + 1) *
      current_pvalue;

    //cout <<"q_value_bh:"<<q_value_bh<<endl;

    while ((decoy_idx < decoy_matches_bonf.numRows()) && 
	   (decoy_matches_bonf.getDouble("p-value bonf.", decoy_idx) <= 
	    current_pvalue)) {
      decoy_idx++;
    }

    double fdr_decoy = 0;
    if (decoy_idx != 0) {
      fdr_decoy = (double)decoy_idx / (double)(target_idx + 1);
    }

    

    //cout <<"Setting row "<<target_idx<<":"<<q_value_bh<<":"<<q_value_decoy<<endl;
    target_matches_bonf.setValue("fdr b-h", target_idx, fdr_bh);
    target_matches_bonf.setValue("fdr decoy", target_idx, fdr_decoy);
  }

  double min_fdr = 1.0;
  for (int idx=target_matches_bonf.numRows()-1;idx>=0;idx--) {
    double cur_fdr = target_matches_bonf.getDouble("fdr b-h", idx);
    if (cur_fdr < min_fdr) {
      min_fdr = cur_fdr;
    }
    target_matches_bonf.setValue("q-value b-h", idx, min_fdr); 
  }

  min_fdr = 1.0;
  for (int idx=target_matches_bonf.numRows()-1;idx>=0;idx--) {
    double cur_fdr = target_matches_bonf.getDouble("fdr decoy", idx);
    if (cur_fdr < min_fdr) {
      min_fdr = cur_fdr;
    }
    target_matches_bonf.setValue("q-value decoy", idx, min_fdr); 
  }

  

  //sort back by scans? or q-value?
  //target_matches_bonf.sortByFloatColumn("q-value decoy");

  string result_file = output_dir + "/search-for-xlinks.qvalues..txt";
  target_matches_bonf.saveData(result_file);

  result_file = output_dir + "/search-for-xlinks.pvalues.decoy.txt";
  decoy_matches_bonf.saveData(result_file);

  return 0;
}
