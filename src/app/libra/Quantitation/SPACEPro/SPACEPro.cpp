// SPACEPro.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include "deep_class.h"
#include "scan_reader.h"
#include <algorithm>
#include <list>
#include <iterator>
#include <cstdio>
#include <numeric>

#define ds_version "1.0.2"
#define ds_builddate "April 26 2022"


using namespace std;

my_parameters cmd_input(int argc, char* argv[]);
void info_description(); 
string marquee();

my_parameters cmd_input(int argc, char* argv[]) {

	my_parameters my_params;
	
	my_params.probability = 0.9;
	my_params.filename = argv[1];
	my_params.ret_time = 2;
	my_params.ppm = 10;
	my_params.cleave_loc = "KR";
	my_params.hyphen = "-";
	my_params.anti_cleave_loc = "";
	my_params.enzymeSense = false;

	for (int i = 2; i < argc; i++) {
		if (string(argv[i]) == "--threshold" || string(argv[i]) == "-t") {
			if (i + 1 < argc && argv[i+1][0]!='-') my_params.probability = atof(argv[++i]);
			else {
				cout << "Invalid threshold" << endl;
				info_description();
				exit(10);
			}
		}
    else if (string(argv[i]) == "--enzymeN" || string(argv[i]) == "-n") {
			my_params.enzymeSense=true;
		}
		else if (string(argv[i]) == "--rtime" || string(argv[i]) == "-r") {
			if (i + 1 < argc && argv[i+1][0]!='-') my_params.ret_time = (float)atof(argv[++i]);
			else {
				cout << "Invalid rtime" << endl;
				info_description();
				exit(10);
			}
		}
		else if (string(argv[i]) == "--ppm" || string(argv[i]) == "-p") {
			if (i + 1 < argc && argv[i+1][0]!='-') my_params.ppm = (float)atof(argv[++i]);
			else {
				cout << "Invalid ppm" << endl;
				info_description();
				exit(10);
			}
		}
		else if (string(argv[i]) == "--loc" || string(argv[i]) == "-c") {
			if (i + 1 < argc && argv[i+1][0]!='-') my_params.cleave_loc = string(argv[++i]);
			else {
				cout << "Invalid loc" << endl;
				info_description();
				exit(10);
			}
		}
		else if (string(argv[i]) == "--anti" || string(argv[i]) == "-a") {
			if (i + 1 < argc && argv[i+1][0]!='-') my_params.anti_cleave_loc = string(argv[++i]);
			else {
				cout << "Invalid anti" << endl;
				info_description();
				exit(10);
			}
		}
		else {
			cout << "Invalid parameter" << endl;
			info_description();
			exit(10);
		}
	}

	return my_params;

}


void info_description() {
  cout << "APPLICATION NOTE: PepXML should contain results from a single mzML search." << endl;
	cout << "USAGE: SPACEPro <pepXML> [OPTIONS]" << endl;
	cout << "OPTIONS:" << endl;
	cout << "  -t, --threshold <number>  =  probability threshold. Default = 0.9" << endl;
	cout << "  -r, --rtime <number>      =  +/- retention time (in minutes) to use for precursor ion extraction. Default = 2.0" << endl;
	cout << "  -p, --ppm <number>        =  +/- mass error (in parts-per-million) to use for precursor ion extraction. Default = 10.0" << endl;
	cout << "  -c, --loc <string>        =  amino acids where cleavage occurs. Default = KR" << endl;
	cout << "  -a, --anti <string>       =  amino acids where cleavage rules are ignored. Default = [none]" << endl;
  cout << "  -n, --enzymeN             =  N-terminal sense of enzyme/chemical. If not specified, sense uses C-terminal as default." << endl;
}

string marquee(){
	string s="SPACEPro, copyright Vidur Kailash, Institute for Systems Biology\n";
	s+="Version: ";
	s+=ds_version;
	s+= "\n";
	s+="Build Date: ";
	s+=ds_builddate;
	s+="\n";
	return s;
}


int main(int argc, char* argv[])
{
	string m=marquee();
	cout << m << endl;  //Announce the application to the user
	if(argc<2) {
		info_description();
		return 2;
	}

	//deep_functions my_deep_functions;
	scan_reader sr;
	//protein_reader pr; 

	my_parameters my_params;
	my_params = cmd_input(argc, argv);

	peptide_lists my_peptide_lists; 
	
	metrics my_metrics;

	cout << "XML File: " << my_params.filename << endl; 

  //Parse the pepXML file and return a list of PSMs
  my_peptide_lists.xml_parse(my_params);
	cout << "XML parsed, "  << my_peptide_lists.total_psms_in_xml << " total PSMs parsed." << endl; 
	my_metrics.psm_xml = my_peptide_lists.total_psms_in_xml;

  //Tally the types of PSMs
	my_peptide_lists.miss_cleave(my_params);
	my_peptide_lists.enzymatic_calc(my_params);

	my_metrics.psm_total=(int)my_peptide_lists.all_psm.size();
	my_metrics.psm_enzymatic=0;
	my_metrics.psm_miscleave=0;
	my_metrics.psm_nonspecific=0;
	for (size_t i = 0; i < my_peptide_lists.all_psm.size(); i++) {
		if (my_peptide_lists.all_psm[i].non_enzymatic) my_metrics.psm_nonspecific++;
		if (my_peptide_lists.all_psm[i].miss_cleaves>0) my_metrics.psm_miscleave++;
		if(my_peptide_lists.all_psm[i].miss_cleaves==0 && !my_peptide_lists.all_psm[i].non_enzymatic) my_metrics.psm_enzymatic++;
	}
	//cout << my_metrics.psm_total << " total PSMs above probability threshold." << endl;
	//cout << "  " << my_metrics.psm_enzymatic << " (" << (double)my_metrics.psm_enzymatic/ my_metrics.psm_total*100 << "%) are enzymatic PSMs." << endl;
	//cout << "  " << my_metrics.psm_miscleave << " (" << (double)my_metrics.psm_miscleave / my_metrics.psm_total * 100 << "%) are mis-cleaved PSMs." << endl;
  //cout << "  " << my_metrics.psm_nonspecific << " (" << (double)my_metrics.psm_nonspecific / my_metrics.psm_total * 100 << "%) are nonspecific PSMs." << endl;
  
	sr.getTIC(my_peptide_lists, my_params);
	if(my_peptide_lists.delete_dup()){
		my_metrics.psm_unique=(int)my_peptide_lists.all_psm.size();
		//cout << "\n" << my_metrics.psm_unique << " precursor ions represent the PSMs." << endl;
		cout << "PSM analysis complete." << endl;
	} else {
    //Right now, the function can never return false
	}

	cout << "Extracting precursor ion signals from " << my_params.mzml << " ..." << endl;
	if(!sr.mzml(my_peptide_lists, my_params)){
		cout << " Error reading mzML file. Exiting." << endl;
		return -1;
	} else {
		cout << " Done!" << endl;
	}

	
	if(my_peptide_lists.reader()){
		my_metrics.pep_unique=(int)my_peptide_lists.all_peptides.size();
		my_metrics.pep_count=0;
		my_metrics.pep_total=0;
		my_metrics.pep_enzymatic=0;
		my_metrics.pep_miscleave=0;
		my_metrics.pep_nonspecific=0;
		for(size_t i=0;i<my_peptide_lists.all_peptides.size();i++){
			my_metrics.pep_total+=my_peptide_lists.all_peptides[i].areaXIC;
			if(my_peptide_lists.all_peptides[i].areaXIC>0) my_metrics.pep_count++;
			if(my_peptide_lists.all_peptides[i].miss_cleaves>0) my_metrics.pep_miscleave += my_peptide_lists.all_peptides[i].areaXIC;
			if (my_peptide_lists.all_peptides[i].non_enzymatic) my_metrics.pep_nonspecific += my_peptide_lists.all_peptides[i].areaXIC;
			if (my_peptide_lists.all_peptides[i].miss_cleaves == 0 && !my_peptide_lists.all_peptides[i].non_enzymatic) my_metrics.pep_enzymatic += my_peptide_lists.all_peptides[i].areaXIC;
		}
		//cout << "\n" << my_metrics.pep_count << " total peptides quantified." << endl;
		//cout << "  " << my_metrics.pep_enzymatic/ my_metrics.pep_total*100 << "% of peptide signal is enzymatic." << endl;
		//cout << "  " << my_metrics.pep_miscleave / my_metrics.pep_total * 100 << "% of peptide signal is mis-cleaved." << endl;
		//cout << "  " << my_metrics.pep_nonspecific / my_metrics.pep_total * 100 << "% of peptide signal is nonspecific." << endl;
		cout << "Peptide analysis complete." << endl;
	}	else {
      //MH: Right now, the function can never return false
	}


	if(my_peptide_lists.prot_stats()){
		my_metrics.prot_count=(int)my_peptide_lists.all_proteins.size();
		my_metrics.prot_avg_enzymatic=0;
		my_metrics.prot_avg_miscleave=0;
		my_metrics.prot_avg_nonspecific=0;
		for (size_t i = 0;i<my_peptide_lists.all_proteins.size();i++){
			my_metrics.prot_avg_enzymatic += (my_peptide_lists.all_proteins[i].sumEnz/ my_peptide_lists.all_proteins[i].total*100);
			my_metrics.prot_avg_miscleave += (my_peptide_lists.all_proteins[i].sumMiss / my_peptide_lists.all_proteins[i].total * 100);
			my_metrics.prot_avg_nonspecific += (my_peptide_lists.all_proteins[i].sumNonSp / my_peptide_lists.all_proteins[i].total * 100);
		}
		my_metrics.prot_avg_enzymatic /= my_metrics.prot_count;
		my_metrics.prot_avg_miscleave /= my_metrics.prot_count;
		my_metrics.prot_avg_nonspecific /= my_metrics.prot_count;
		//cout << "\n" << my_metrics.prot_count << " total proteins represented by one or more proteotypic peptides." << endl;
		//cout << "  " << my_metrics.prot_avg_enzymatic << "% average enzymatic percentage among proteins." << endl;
		//cout << "  " << my_metrics.prot_avg_miscleave << "% average mis-cleavage percentage among proteins." << endl;
		//cout << "  " << my_metrics.prot_avg_nonspecific << "% average nonspecific percentage among proteins." << endl;
		cout << "Protein analysis complete." << endl;

		string of;
		of=my_params.filename;
		of+=".sp.json";
		my_peptide_lists.json(of); //MH: I had to move this up here because there are horrid memory leaks that occur otherwise.
	}
	else{}


	//Export report
	cout << "\nReport:\n" << endl;
	my_peptide_lists.calc(my_metrics);
	my_peptide_lists.print(my_metrics,my_params,m);
	
	return 0;

}
