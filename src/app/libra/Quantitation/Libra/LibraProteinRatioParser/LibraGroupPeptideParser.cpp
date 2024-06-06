/*
  Program       : LibraGroupPeptideParser
  Author        : original author: Andy Keller;  additions and maintenance: Nichole King
  Date          : 10.16.05
  SVN info      : $Id: LibraGroupPeptideParser.cpp 8896 2023-03-21 07:00:53Z real_procopio $

  Class to retrieve peptide info from pepXML file, and calculate quantitation.

  Copyright (C) 2005

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "Common/sysdepend.h"
#include "Common/TPPVersion.h" // contains version number, name, revision

#include "LibraGroupPeptideParser.h"
#include "StringConvertor.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <utility>
#include "Util/RACI/RACI.h" // for reading gzipped files with efficient seeks

// overloaded  constructor to handle single xmlfile
LibraGroupPeptideParser::LibraGroupPeptideParser(const char * xmlfile, Array<const char*>* peptides, 
						 double minprob, int norm_channel) : Parser(NULL) {
  parse_all_ = False;
  handleSingleXmlFileAttributes( xmlfile );
  initializeClassAttributes( peptides, minprob, norm_channel, (double)-99.  );

  // calls parse:
  init(NULL);
}

// overloaded  constructor to handle single xmlfile and threshhold
LibraGroupPeptideParser::LibraGroupPeptideParser(const char * xmlfile, Array<const char *>* peptides,
						 double minprob, int norm_channel, double minThreshInten) : Parser(NULL) {
  parse_all_ = False;
  handleSingleXmlFileAttributes( xmlfile );
  initializeClassAttributes( peptides, minprob, norm_channel, minThreshInten );
  init(NULL);
}


// overloaded constructor for array of xml files and threshhold intensity
LibraGroupPeptideParser::LibraGroupPeptideParser(Array<const char *>* xmlfiles, Array<const char *>* peptides,
						 double minprob, int norm_channel, double minThreshInten) : Parser(NULL) {
  parse_all_ = False;
  xmlfiles_ = xmlfiles;
  single_input_ = False;
  initializeClassAttributes( peptides, minprob, norm_channel, minThreshInten );
  init(NULL);
}


// overloaded constructor for array of xmlfiles
LibraGroupPeptideParser::LibraGroupPeptideParser(Array<const char *>* xmlfiles, Array<const char *>* peptides,
						 double minprob, int norm_channel) : Parser(NULL) {
  parse_all_ = False;
  xmlfiles_ = xmlfiles;
  single_input_ = False;
  initializeClassAttributes( peptides, minprob, norm_channel, (double)-99.  );
  init(NULL);
}


LibraGroupPeptideParser::LibraGroupPeptideParser(Array<const char *>* xmlfiles, const char * quantitationFileName,
						 double minprob, int norm_channel, double minThreshInten) : Parser(NULL) {
  parse_all_ = True;
  xmlfiles_ = xmlfiles;
  single_input_ = False;
  quantitation_file_name_ = quantitationFileName;
  initializeClassAttributes( NULL, minprob, norm_channel, minThreshInten );
  init(NULL);
}


void LibraGroupPeptideParser::handleSingleXmlFileAttributes( const char * xmlFile ) {
  xmlfiles_ = new Array<const char *>;
  char* next = new char[strlen(xmlFile)+1];
  strcpy(next, xmlFile);
  xmlfiles_->insertAtEnd(next);
  single_input_ = True;
}


// initialize/define class attributes
void LibraGroupPeptideParser::initializeClassAttributes( Array<const char *>* peptides, double minprob, 
							 int norm_channel, double thresh) {

  peptides_ = peptides;
  min_probability_ = minprob;
  norm_channel_ = norm_channel;
  minIntensity = thresh;
  channel_masses_ = new Array<double>;
  summary_tags_ = new Array<Tag*>;

  ratio_log_sum_ = NULL;
  ratio_sum_ = NULL;
  protein_ratio_array_ = NULL;
  protein_mass_array_ = NULL;
  protein_stdev_array_ = NULL;
  protein_se_array_ = NULL;
  protein_variance_array_ = NULL;
  protein_ratio_array_wrt_refchannel_ = NULL;
  protein_error_array_wrt_refchannel_ = NULL;

  ratio_num_ = 0;
  ratio_num_filtered_ = 0;

  tmpInf = 99.99;
  tmpNoQuantitation = -9.0;
}


LibraGroupPeptideParser::~LibraGroupPeptideParser() {
  clearProteinRatioData();

  if(single_input_) {
    for(int k = 0; k < xmlfiles_->length(); k++)
      if((*xmlfiles_)[k] != NULL)
	delete (char *)(*xmlfiles_)[k]; // cast away const, we allocated it
    delete xmlfiles_;
  }

  for (map <string, vector<psm> >::const_iterator it = all_peptides_map_.begin(); it != all_peptides_map_.end(); it++) {
    for (vector<psm>::const_iterator p = it->second.begin(); p != it->second.end(); p++) {
      delete[] p->intensities;
    }
  }

  delete channel_masses_;
}


void LibraGroupPeptideParser::clearProteinRatioData () {
  if (ratio_log_sum_ != NULL) delete[] ratio_log_sum_;
  if (ratio_sum_ != NULL) delete[] ratio_sum_;
  if (protein_ratio_array_ != NULL) delete[] protein_ratio_array_;
  if (protein_mass_array_ != NULL) delete[] protein_mass_array_;
  if (protein_stdev_array_ != NULL) delete[] protein_stdev_array_;
  if (protein_se_array_ != NULL) delete[] protein_se_array_;
  if (protein_variance_array_ != NULL) delete[] protein_variance_array_;

  if (protein_ratio_array_wrt_refchannel_ != NULL) 
    delete[] protein_ratio_array_wrt_refchannel_;

  if (protein_error_array_wrt_refchannel_ != NULL) 
    delete[] protein_error_array_wrt_refchannel_;
}


// initialize storage classes used in parse and delete last elements if they exist
// @param number of reagent channels
void LibraGroupPeptideParser::initializeStorageWithinParse ( int nChan ) {
  clearProteinRatioData();

  // initialize vars
  ratio_log_sum_ = new double[nChan];

  ratio_sum_ = new double[nChan];

  protein_ratio_array_ = new double[nChan];
  protein_mass_array_ = new double[nChan];
  protein_stdev_array_ = new double[nChan];
  protein_se_array_ = new double[nChan];
  protein_variance_array_ = new double[nChan];
  protein_ratio_array_wrt_refchannel_ = new double[nChan];
  protein_error_array_wrt_refchannel_ = new double[nChan];

  peptide_intensities_map.clear();
  peptide_intensities_prenorm_map.clear();

  peptide_sequence_map.clear();
  peptide_keptFlag_map.clear();

  ratio_num_ = 0;
  ratio_num_filtered_ = 0;

  for(int k = 0; k < nChan; k++) {
    ratio_log_sum_[k] = 0.0;
    ratio_sum_[k] = 0.0;
    protein_ratio_array_[k] = 0.0;
    protein_mass_array_[k] = 0.0;
    protein_stdev_array_[k] = 0.0;
    protein_se_array_[k] = 0.0;
    protein_variance_array_[k] = 0.0;
    protein_ratio_array_wrt_refchannel_[k] = 0.0;
    protein_error_array_wrt_refchannel_[k] = 0.0;
  }
}


// parse pepXML file
// @param xmlfile pepXML file
void LibraGroupPeptideParser::parse(const char * xmlfile) {
  Tag* tag = NULL;
  char *nextline = new char[line_width_];   //line_width_ is 1000000
  char* data = NULL;
  Boolean analyze = False;
  Boolean search_score_found = False;
  double probability=-1;

  char libra_code[500]; // looks like this is no longer used meaningfully - LM
  libra_code[0] = 0;

  //for(int k = 0; k < peptides_->length(); k++)
  //  cout << (*peptides_)[k] << " ";
  // cout << endl;

  int summary_analyze = 0;
  double* next_intensities = NULL;

  //xxxxxxx Think Andy's design was to update is_rejected for each
  // peptide in the interact.xml file, but not doing that here currently

  // whether or not user has specified not to use peptide LIBRA quant
  Boolean rejected = False; 

  cout << "Reading pepXML data..." << std::flush;
  // for each input file:
  for(int k = 0; k < xmlfiles_->length(); k++) {
    int curr_index = 0;
    RACI fin((*xmlfiles_)[k]); // can read gzipped xml

    if(! fin) {
      cout << "LibraGroupPeptideParser: error opening " << (*xmlfiles_)[k] << endl;
      exit(1);
    }

    while(fin.getline(nextline, line_width_)) {
      if(strstr(nextline, "<libra_summary") != NULL || summary_analyze == 1 || analyze || possiblePeptideListMember(nextline)) {
	data = strstr(nextline, "<");

	while(data != NULL) {
	  tag = new Tag(data);

	  if(tag != NULL) {
	    if(tag->isStart() && ! strcmp(tag->getName(), "libra_summary")) { 

	      summary_analyze++;

	      if(strlen(libra_code) == 0)
		strcpy(libra_code, tag->getAttributeValue("channel_code"));

	      else if(strcmp(libra_code, tag->getAttributeValue("channel_code"))) {
		cout << "error: different libra channels used: " << libra_code 
		     << " vs " << tag->getAttributeValue("channel_code") << endl;
		exit(1);
	      }

	      if(summary_analyze == 1 && peptides_ == NULL) {
		char text1[10];
		sprintf(text1, "%d", norm_channel_);
		tag->setAttributeValue("normalization", text1);
	      }
	    } // end if  start tag is libra_summary

	    // store libra_summary attributes:
	    if(summary_analyze == 1) { 
	      summary_tags_->insertAtEnd(tag->copy()); // grab it

	      if(tag->isStart() && ! strcmp(tag->getName(), "fragment_masses")) {
		channel_masses_->insertAtEnd(atof(tag->getAttributeValue("mz")));

	      } else if(tag->isEnd() && ! strcmp(tag->getName(), "libra_summary")) {
		summary_analyze++; // done

		if(norm_channel_ > channel_masses_->length()) {
		  cout << "error: norm channel " << norm_channel_ << 
		    " exceeds actual number: " << channel_masses_->length() << endl;
		  exit(1);
		}

		// initialize storage classes
		next_intensities = new double[channel_masses_->length()];

		for(int k = 0; k < channel_masses_->length(); k++) {
		  next_intensities[k] = 0.0;
		}

		initializeStorageWithinParse ( channel_masses_->length() );
	      }
	      
	      // parse and store <search_hit> contents if hit_rank=1 and peptide is in list
	    } else if(tag->isStart() && ! strcmp(tag->getName(), "search_hit") && 
		      ! strcmp(tag->getAttributeValue("hit_rank"), "1")
		      && peptideListMember(tag->getAttributeValue("peptide"))) {

	      analyze = True;

	      if (dbug) {
		cout << "found a match for " << 
		  tag->getAttributeValue("peptide") << endl;
	      }

	      current_peptide = tag->getAttributeValue("peptide");
	      current_protein = tag->getAttributeValue("protein"); 

	    } else if(analyze) {
	      if(tag->isEnd() && ! strcmp(tag->getName(), "libra_result")) {
		curr_index = 0;

	      } else if(tag->isStart() && ! strcmp(tag->getName(), "libra_result")) {
		rejected = tag->getAttributeValue("is_rejected") != NULL &&
		  ! strcmp(tag->getAttributeValue("is_rejected"), "1");

	      } else if(tag->isStart() && ! strcmp(tag->getName(), "intensity")) {
		//store target_mass
		protein_mass_array_[curr_index] = atof(tag->getAttributeValue("target_mass"));

		// use absolute intensities from interact.xml peptide quantitation:
		next_intensities[curr_index] = atof(tag->getAttributeValue("absolute"));

		if (0 && dbug) {
		  cout << "mass: " << protein_mass_array_[curr_index] << " :: int: "
		       << next_intensities[curr_index] <<endl;
		}

		// if reagent m/z line isn't found, setting m/z to expected 
		// reagent value and setting intensity to zero
		if ( (protein_mass_array_[curr_index] < 1.)) {
		  protein_mass_array_[curr_index] = 
		    (*channel_masses_)[curr_index];

		  next_intensities[curr_index] = 0;
		}

		curr_index++;

	      } else if(tag->isStart() && ! strcmp(tag->getName(), "peptideprophet_result")) {
		probability = atof(tag->getAttributeValue("probability"));

	      } else if(tag->isEnd() && ! strcmp(tag->getName(), "search_hit")) { // process
		// include in protein ratio calculation
		if(! rejected && (min_probability_ == 0.0 || probability >= min_probability_)) { 
		  ratio_num_++; // move to else clause below?

		  if (parse_all_) {
		    psm pep;
		    pep.sequence = current_peptide;
		    pep.protein  = current_protein;
		    pep.num_channels = channel_masses_->length();

		    pep.intensities = new double[channel_masses_->length()];
		    
		    
		    for(int k = 0; k < channel_masses_->length(); k++) {
		      pep.intensities[k] = next_intensities[k];
		    }

		    all_peptides_map_[pep.sequence].push_back(pep);
		  }
		  else {
		    // store peptide intensities and sequence in std::maps, to access later by key number
		    storeAndNormalizeValues (next_intensities, channel_masses_->length() );
		  }

		} // if passes test

		probability = -1.0;
		analyze = False;
		rejected = False;
		current_peptide.erase();
	      } // end process search hit
	    } // end analyzed

	    delete tag;
	  } // end no null

	  data = strstr(data+1, "<");
	} // end next tag
      } // end possible peptide present
    } // end next line

    fin.close();
  } // next inputfile
  cout << "ok" << endl << std::flush;  // done reading pepXML data

  if (!parse_all_) {
    // Calculate basic stats for use in outlier removal
    calculateBasicStats( channel_masses_->length() );

    /*
     * Re-calculate Ratios by not including peptides that deviate from the mean
     * by more than 2 sigma.  
     * If other attributes have been set, such as a minimum threshhold intensity,
     * that's handled here too.
     * Stores kept or removed flag for each peptide.  
     */
    recalculateRatios( channel_masses_->length() );

    /*
     * Appends peptide and protein quantitation to existing output file
     * called quantitation.tsv.  The outfile is created by LibraProteinRatioParser
     * and appended to here.  Started to create a class to handle the File
     * writing, but haven't finished.  It's the QuantitationFile in this package.
     */
    if (current_protein.size() > 0) {
      writePeptideAndProteinQuantitationToOutfile( channel_masses_->length() );
    }
  }

  if (next_intensities != NULL)  delete[] next_intensities;
  delete[] nextline;
}


void LibraGroupPeptideParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "search_hit")) {
    if(tag->isStart() && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) {
      //tag->print();
      filter_ = True;
    }else{
      if(filter_ && tag->isEnd())
        filter_memory_ = True;
    }
  }
}


Array<Tag*>* LibraGroupPeptideParser::getProtXMLSummaryTags(double minpepprob, double minpepwt, double minprotprob) {
  if(summary_tags_ == NULL)
    return summary_tags_;
  char text[200];
  for(int k = 0; k < summary_tags_->length(); k++) {
    if((*summary_tags_)[k]->isStart() && ! strcmp("libra_summary", (*summary_tags_)[k]->getName())) {
      sprintf(text, "%0.1f", minpepprob);
      (*summary_tags_)[k]->setAttributeValue("min_pep_prob", text);
      sprintf(text, "%0.1f", minpepwt);
      (*summary_tags_)[k]->setAttributeValue("min_pep_wt", text);
      sprintf(text, "%0.1f", minprotprob);
      (*summary_tags_)[k]->setAttributeValue("min_prot_prob", text);
      (*summary_tags_)[k]->setAttributeValue("version", szTPPVersionInfo);

      k = summary_tags_->length();
    }
  }
  return summary_tags_;
}


Array<Tag*>* LibraGroupPeptideParser::getProtXMLTags(const char* prot_name, Array<const char *>* peptides) {
  Boolean verbose = False;
  initializeStorageWithinParse ( channel_masses_->length() );
  peptides_ = NULL;
  peptides_ = peptides;

  current_protein = prot_name;

  for(int k = 0; k < peptides_->length(); k++) {
    if (all_peptides_map_.find((*peptides_)[k]) == all_peptides_map_.end()) {
      if(verbose)
	cout << endl << "unmatched: " << (*peptides_)[k] << endl;
    }
    else {
      for (int i = 0; i < all_peptides_map_.at((*peptides_)[k]).size(); i++) {
	if(verbose)
	  cout << "found one psm: " << all_peptides_map_.at((*peptides_)[k])[i].sequence << endl;

	ratio_num_++;
	current_peptide = all_peptides_map_.at((*peptides_)[k])[i].sequence;
	storeAndNormalizeValues (all_peptides_map_.at((*peptides_)[k])[i].intensities, channel_masses_->length() );
      }
    }
  }

  calculateBasicStats( channel_masses_->length() );  // Calculate basic stats for use in outlier removal
  recalculateRatios( channel_masses_->length() );    // Recalculate Ratios by not including peptides that deviate from the mean more than 2-sigma

  if (current_protein.size() > 0) {
    writePeptideAndProteinQuantitationToOutfile( channel_masses_->length() );  // Append peptide and protein quantitation to quantitation.tsv
  }

  return getProtXMLTags();
}

Array<Tag*>* LibraGroupPeptideParser::getProtXMLTags() {
  char text[300];
  Array<Tag*>* output = new Array<Tag*>;
  Tag* next = NULL;

  next = new Tag("libra_result", True, False);

  // sprintf(text, "%d", ratio_num_);
  // replacing number of all peptides with the number passing filters:
  sprintf(text, "%d", ratio_num_filtered_);

  next->setAttributeValue("number", text);
  output->insertAtEnd(next);
 
						      // compute protein ratios for each channel and write to protxml
  for(int k = 0; k < channel_masses_->length(); k++) {
    next = new Tag("intensity", True, True);

    sprintf(text, "%d", k+1);
    next->setAttributeValue("channel", text);

    sprintf(text, "%0.5f", (*channel_masses_)[k]);
    next->setAttributeValue("mz", text);

    double ratio, error;

    if(ratio_num_filtered_ == 0) { // case for no peptides
	ratio = tmpNoQuantitation;
	error = tmpNoQuantitation;

    } else {
      // if user selected normalization:
      if (norm_channel_ > 0) {
	ratio = protein_ratio_array_wrt_refchannel_[k];
	// using 1 sigma errors:
	error = protein_error_array_wrt_refchannel_[k];

      } else {
	  ratio = protein_ratio_array_[k];
	
	// using 1 sigma errors:
	error = protein_se_array_[k];
      }

      if ( ratio_num_filtered_ == 1) {  // case for 1 peptide
	// until have a term for noise due to signal in the final errors, need
	// to use alarmingly high error to indicate infinite standard error
	error = tmpInf;
      }
    }

    // cout << ratio << " +- " << error << endl;

    if (isnan(ratio)) {
      next->setAttributeValue("ratio", "NaN");  // xml-appropriate value
    }
    else if (isinf(ratio)) {
      next->setAttributeValue("ratio", "INF");  // xml-appropriate value
    }
    else {
      sprintf(text, "%0.2f", ratio);
      next->setAttributeValue("ratio", text);
    }
    sprintf(text, "%0.2f", error);
    next->setAttributeValue("error", text);

    output->insertAtEnd(next);

  } // next channel

  next = new Tag("libra_result", False, True);
  output->insertAtEnd(next);

  return output;
}


double LibraGroupPeptideParser::getRatioLogSum(int index) {
  if(index < 0 || index >= channel_masses_->length()) {
    cout << "error: index " << index << " vs " << channel_masses_->length() << endl;
    exit(1);
  }
  return ratio_log_sum_[index];
}


int LibraGroupPeptideParser::getRatioNum() {
  return ratio_num_;
}

// get the number of protein's peptides that pass filters 
// @return ratio_num_filtered
int LibraGroupPeptideParser::getRatioNumFiltered() {
  return ratio_num_filtered_;
}

int LibraGroupPeptideParser::getNumChannels() {
  if(channel_masses_ == NULL)
    return 0;
  return channel_masses_->length();
}


Boolean LibraGroupPeptideParser::peptideListMember(const char* pep) {
  Boolean verbose = False; //strstr(pep, "YEFCTILKK") != NULL;
  if(verbose) {
    cout << "comparing " << pep << " with peptides....";
    for(int k = 0; k < peptides_->length(); k++)
      cout << "=" << (*peptides_)[k] << "=";
    cout << endl;
  }
  if(parse_all_)
    return True;
  if(peptides_ == NULL)
    return False;
  for(int k = 0; k < peptides_->length(); k++)
    if(! strcmp((*peptides_)[k], pep)) {
      if(verbose)
	cout << "returning true" << endl;
      return True;
    }
  if(verbose) {
    cout << "-" << pep << "-" << (*peptides_)[0] << "-" << endl;
    cout << strlen(pep) << " vs " << strlen((*peptides_)[0]) << endl;
  }

  return False;
}


Boolean LibraGroupPeptideParser::possiblePeptideListMember(const char* data) {
  if (parse_all_)
    return True;
  if(peptides_ == NULL || data == NULL)
    return False;
  for(int k = 0; k < peptides_->length(); k++)
    if(strstr(data, (*peptides_)[k]) != NULL)
      return True;
  return False;
}

// store values of peptide that passes filters in parse
// @param intensities pointer to vector of intensities for the spectrum reagent channels
// @param number of reagent channels
void LibraGroupPeptideParser::storeAndNormalizeValues (double* intensities, int nChannels ) {

  double nextratio;
  std::vector<double>* peptide_intensities_vector;   // vector to hold a peptide's channel intensities
  std::vector<double>* peptide_intensities_prenorm_vector;   // vector to hold a peptide's pre-normalized channel intensities

  peptide_intensities_vector = new std::vector<double>;
  peptide_intensities_prenorm_vector = new std::vector<double>;

  /*
   * store peptide sequence for write to outfile later:
   * key = peptide number
   * value = (string) current peptide
   */
  peptide_sequence_map.insert( std::make_pair ( ratio_num_,  current_peptide ) );

  // Normalize intensities by the sum of all channels:
  //sum of intensities of channels 0 through n-1.  used to normalize each channel
  double channel_intensities_sum = 0.;
  int channel;


  
  for(channel = 0; channel < nChannels; channel++) {
    channel_intensities_sum += intensities[channel];
  }

  // store normalized peptide channel intensities in vector:
  for(channel = 0; channel < nChannels; channel++) {

    (*peptide_intensities_vector).push_back( channel_intensities_sum ?
					     ( intensities[channel] / channel_intensities_sum ) : 0.0);

    (*peptide_intensities_prenorm_vector).push_back( channel_intensities_sum ?
						     intensities[channel] : 0.0);
  }

  /*
   * store intensitiy vectors in maps with
   * key = peptide number
   * value = vector<double>* of intensities
   */
  peptide_intensities_map.insert( std::make_pair 
				  (ratio_num_,  (*peptide_intensities_vector) ) );

  peptide_intensities_prenorm_map.insert( std::make_pair 
					  (ratio_num_,  (*peptide_intensities_prenorm_vector) ) );

  // make log10 sum of normalized intensities
  for(channel = 0; channel < nChannels; channel++) {
    nextratio = (*peptide_intensities_vector)[channel];

    if(nextratio > pow(10.0, MAX_LOG))
      nextratio = pow(10.0, MAX_LOG);

    else if(nextratio < pow(10.0, -1 * MAX_LOG))
      nextratio = pow(10.0, -1 * MAX_LOG);

    ratio_log_sum_[channel] += log10(nextratio);

    ratio_sum_[channel] += log10(nextratio);

    if (dbug) {
      cout 
	<< current_protein.c_str() << "  " 
	<< current_peptide.c_str() << "  "
	<< "  peptide number (ratio_num_): " << ratio_num_
	<< "  channel: " << channel+1
	<< "  inten: " << intensities[channel]
	<< "  ratio: " << nextratio
	<< "  channel_intensities sum: " << channel_intensities_sum
	<< "  ratio_log_sum_" << ratio_log_sum_[channel] 
	<< endl;
    }
  } //next channel

  delete peptide_intensities_vector;
  delete peptide_intensities_prenorm_vector;
}


// calculate mean and standard deviation of the mean for subsequent use in outlier removal
// @ number of reagent channels
void LibraGroupPeptideParser::calculateBasicStats( int nChannels ) {
  if (dbug) {
    if (ratio_num_ > 0) {
      cout << "-------DONE w/ 1st pass of all peps in protein (before filtering): "
	   << current_protein.c_str() << "  " 
	   << "  npeps:" << ratio_num_
	   << "   check ratio_log_sum_[1]: " << ratio_log_sum_[1] << "\n" << endl;
    }
  }

  // fill protein_ratio_array_ and initialize 
  //protein_variance_array_[channel]
  int channel;
  for (channel=0; channel < nChannels; channel++) {
    double logmean = ratio_log_sum_[channel] / ratio_num_;
    protein_ratio_array_[channel] = (pow((double) 10,  (double) logmean));
    protein_variance_array_[channel] = 0.0;
  }

  std::map<int, std::vector<double> >::iterator iter;
  std::vector<double> tmpIntensityVector;

  // determine variance
  for (int pep_num = 1; pep_num <= ratio_num_; pep_num++) {
    iter = peptide_intensities_map.find( pep_num );

    // assign peptide reagent intensity vector to tmpIntensityVector
    if (iter != peptide_intensities_map.end())
      tmpIntensityVector = iter->second;

    for (int channel=0; channel < (int) tmpIntensityVector.size(); channel++) {
      protein_variance_array_[channel] +=  
	( (pow( (tmpIntensityVector[channel] - 
		 protein_ratio_array_[channel]), (double) 2. ))/ (ratio_num_ - 1)) ;
    }
  }

  // determine stdev
  for (channel=0; channel < (int) tmpIntensityVector.size(); channel++) {
    if (ratio_num_ == 1) {

      // until have a term for noise due to signal in the final errors, need
      // to use alarmingly high error to indicate infinite standard error
      protein_stdev_array_[channel] = tmpInf;
      protein_variance_array_[channel] = tmpInf;

    } else {
      protein_stdev_array_[channel] = sqrt ( protein_variance_array_[channel] );
    }
  }
}


/**
 * Re-calculate Ratios by not including peptides that deviate from the mean
 * by more than 2 sigma;
 * if other attributes have been set, such as a minimum threshhold intensity,
 * that's handled here too;
 * and stores kept or removed flag for each peptide.  
 * @param nChannels number of reagents speficied in condition xml file
 */
void LibraGroupPeptideParser::recalculateRatios( int nChannels ) {
  int n = nChannels;

  if (dbug) {
    if (ratio_num_ > 0) {
      for (int channel=0; channel < n; channel++) {
	cout << "------(before filtering) channel " << channel+1
	     << current_protein.c_str() << "  " 
	     << "  ratio: " << protein_ratio_array_[channel]
	     << "  variance: " << protein_variance_array_[channel] 
	     << "  stdev: " << protein_stdev_array_[channel] 
	     << "  npeps: " << ratio_num_ << endl;
      }
    }
  }

  // declare and init structures to use in recomputing ratios and errors:

  std::vector<double> peptide_intensities_vector(n);
  std::vector<double> tmp_ratio_log_sum_(n);
  int tmp_ratio_num_ = 0;
  int channel;

  for (channel=0; channel < n; channel++) {
    tmp_ratio_log_sum_[channel] = 0.0;
    peptide_intensities_vector[channel] = 0.0;
  }

  // declare tmp std::map to store intensities that are within 2 sigma
  std::map<int, std::vector<double> > tmp_peptide_intensities_map;
  std::vector<double> tmpIntensityVector;
  std::map<int, std::vector<double> >::iterator iter;

  // For entries more than 2 sigma from the mean, don't use; set removedFlag
  int pep_num;
  for (pep_num = 1; pep_num <= ratio_num_; pep_num++) {
    // flag used for removal; outlier=0 means keep
    int outlier = 0;

    iter = peptide_intensities_map.find( pep_num );

    // assign peptide reagent intensity vector to tmpIntensityVector
    if (iter != peptide_intensities_map.end())
      tmpIntensityVector = iter->second;

    // pre-normalized intensities:
    std::vector<double> tmpVector;
    std::map<int, std::vector<double> >::iterator iter2;

    iter2 =  peptide_intensities_prenorm_map.find( pep_num);

    if (iter2 != peptide_intensities_prenorm_map.end())
      tmpVector = iter2->second;

    double totIntensity = 0.0;
    for (int channel=0; channel < (int) tmpIntensityVector.size(); channel++) {
      totIntensity += tmpIntensityVector[channel];
      // protein_ratio_array_[channel] is the mean value for that channel
      double diff = tmpIntensityVector[channel] - protein_ratio_array_[channel];

      // 2 sigma:
      double max_diff = 2.*( protein_stdev_array_[channel] );

      if (diff < 0.)
	diff = diff * -1.0;

      if (max_diff < 0.)
	max_diff = max_diff * -1.0;


      // if peptide is more than 2 sigma from the mean, mark as outlier      
      if ( ( diff > max_diff ) && (ratio_num_ > 1) )
	outlier = 1;
      
      // if there's a threshhold set, use it:
      if (minIntensity >= -0.00001 ) {
	if (tmpVector[channel] < minIntensity) {
	  outlier = 1;
	}
      }

      /*
      //mark as removed, entries where reagent channel wasn't found.
      // this is handled in parse now by replacing values
      if ( (protein_mass_array_[channel] < 1.)) {
      outlier = 1;
      }
      */

    } // end iteration over channels for this peptide pep_num

    // ignore if no intensity in all channels
    if (totIntensity == 0.0)
      outlier = 1;

    // count revised number of peptides and recompute sum of channels

    if (outlier == 0) {
      tmp_ratio_num_++;

      for (int channel=0; channel < (int) tmpIntensityVector.size(); channel++) {
	if(tmpIntensityVector[channel] < pow(10.0, -1 * MAX_LOG))
	  tmpIntensityVector[channel] = pow(10.0, -1 * MAX_LOG);

	tmp_ratio_log_sum_[channel] += log10( tmpIntensityVector[channel]  );

	// store intensities for this peptide in a  vector
	peptide_intensities_vector[channel] = (tmpIntensityVector[channel]);

	if (dbug) {
	  cout << "(after filtering) rn: " << tmp_ratio_num_
	       << "  ch: " << channel+1
	       << "  nxt_r: " << tmpIntensityVector[channel]
	       << "  r_log_sum_" << tmp_ratio_log_sum_[channel] 
	       << endl;
	}
      }

      // store intensities in map with key = tmp_ratio_num_, value = vector
      tmp_peptide_intensities_map.insert( std::make_pair 
					  (tmp_ratio_num_, peptide_intensities_vector ) );

      if (dbug) {
	cout << tmp_ratio_num_ <<endl;
      }
    } // end if outlier == 0

    std::string kept = "Yes";

    if (outlier == 1)
      kept = "No";
    else if (outlier == 0)
      kept = "Yes";
    else
      kept = "hmmm";

    peptide_keptFlag_map.insert( std::make_pair(pep_num, kept) );

    if (dbug) {
      std::string sequence;
      std::map<int, std::string >::iterator it1;

      it1 = peptide_sequence_map.find( ratio_num_ );
      if ( it1 != peptide_sequence_map.end() )
	sequence = it1->second;

      cout << current_protein.c_str() << ":" <<sequence.c_str()<<":"<< ratio_num_
	   <<":(" << outlier <<"):"<<kept.c_str()<<endl;
    }

  } // end loop over all peptides 

  // if at least one peptide had outlier = 0:
  if (tmp_ratio_num_ > 0) {
    // Recalculate mean

    // store new values in protein_ratio_array_ and 
    // initialize protein_variance_array_[channel]
    for (int channel=0; channel < channel_masses_->length(); channel++) {
      double logmean = tmp_ratio_log_sum_[channel] / tmp_ratio_num_;

      // for nan from negative intensities, set ratio to 0.
      if( !(logmean > 0.) && !(logmean < 0.) )
	protein_ratio_array_[channel] = 0;
      else
	protein_ratio_array_[channel] = (pow((double)10,  logmean));

      protein_variance_array_[channel] = 0.0;
    }
  }

  tmpIntensityVector.clear();

  // determine variance
  for (pep_num = 1; pep_num <= tmp_ratio_num_; pep_num++) {
    iter = tmp_peptide_intensities_map.find( pep_num );

    // assign peptide reagent intensity vector to tmpIntensityVector
    if (iter != tmp_peptide_intensities_map.end())
      tmpIntensityVector = iter->second;


    for (int channel=0; channel < (int) tmpIntensityVector.size(); channel++) {
      protein_variance_array_[channel] +=
	( (sqr( (tmpIntensityVector[channel] - protein_ratio_array_[channel]) ))
	  / (tmp_ratio_num_ - 1)) ;
    }
  }

  // determine stdev and normalize w.r.t. reference channel
  for (channel=0; channel < channel_masses_->length(); channel++) {

    //xxx      if (tmp_ratio_num_ < 2 || (protein_ratio_array_[channel] == 0) ) {
    if (tmp_ratio_num_ < 2) {

      // until have a term for noise due to signal in the final errors, need
      // to use alarmingly high error to indicate infinite standard error
      protein_variance_array_[channel] = tmpInf;
      protein_stdev_array_[channel] = tmpInf;
      protein_se_array_[channel] = tmpInf;
      protein_error_array_wrt_refchannel_[channel] = tmpInf;

    } else {

      protein_stdev_array_[channel] = sqrt ( protein_variance_array_[channel] );
      protein_se_array_[channel] = protein_stdev_array_[channel] / sqrt( (double) tmp_ratio_num_);

      // calc ref channel error as the sum of the standard errors (1 sigma)
      // of channel and ref channel added in quadrature:
      if (norm_channel_ > 0) {
	double ts = (pow(protein_se_array_[channel], (double)2.)) +
	  (pow(protein_se_array_[norm_channel_ - 1], (double)2.));

	protein_error_array_wrt_refchannel_[channel] = sqrt( ts );
      }
    }

    if (norm_channel_ > 0)
      protein_ratio_array_wrt_refchannel_[channel] = 
	protein_ratio_array_[channel] / 
	protein_ratio_array_[norm_channel_ - 1];

  }

  if( tmp_ratio_num_ == 0 ) {
    //else no peptides were kept, update arrays
    for (int channel=0; channel < channel_masses_->length(); channel++) {
      protein_ratio_array_[channel] = tmpNoQuantitation;
      protein_variance_array_[channel] = tmpNoQuantitation;
      protein_stdev_array_[channel] = tmpNoQuantitation;
      protein_se_array_[channel] = tmpNoQuantitation;
      protein_error_array_wrt_refchannel_[channel] = tmpNoQuantitation;
      protein_ratio_array_wrt_refchannel_[channel] = tmpNoQuantitation;
    }
  }

  if (dbug) {
    for (int channel=0; channel < channel_masses_->length(); channel++) {
      cout << "------(after filtering) channel " << channel+1
	   << "  ratio: " << protein_ratio_array_[channel]
	   << "  variance: " << protein_variance_array_[channel] 
	   << "  stdev: " << protein_stdev_array_[channel] 
	   << "  npeps: " << tmp_ratio_num_ 
	   << "  ratio wrt ref: " 
	   << protein_ratio_array_wrt_refchannel_[channel]
	   << endl;
    }
  }

  ratio_num_filtered_ = tmp_ratio_num_;

  //xxxxxxxxxx
  // should be updating a class attribute of lists of peptides to update the peptide element
  // in protXML
}

/*
 * Appends peptide and protein quantitation to existing output file
 * called quantitation.tsv.  The outfile is created by LibraProteinRatioParser
 * and appended to here.  Started to create a class to handle the File
 * writing, but haven't finished.  It's the QuantitationFile in this package.
 * @param number of reagent channels
 */
void LibraGroupPeptideParser::writePeptideAndProteinQuantitationToOutfile( int nChannels ) {

  // formatting info to append to "quantitation.tsv"
  // if want clearer details of format, please see QuantitationFile code

  //local data structures:
  std::vector<std::string> peptide_quantitation_lines;
  std::vector<std::string> protein_quantitation_lines;
  std::string line;


  for (int pep_num = 1; pep_num <= ratio_num_; pep_num++) {

    // expecting peptide line to be of format:
    //prot<tab>pep<tab>nr1<tab>nr2<tab>nr3<tab>nr4<tab>in1<tab>in2<tab>in3<tab>in4<tab>is_rejected

    // retrieve sequence from std::map:
    std::string sequence;
    std::map<int, std::string >::iterator it1;

    it1 = peptide_sequence_map.find( pep_num );

    if ( it1 != peptide_sequence_map.end() )
      sequence = it1->second;

    line = current_protein + std::string("\t") + sequence;

    // retrieve normalized intensities from std::map:
    std::vector<double> intensities;
    std::map<int, std::vector<double> >::iterator it2;

    it2 = peptide_intensities_map.find( pep_num );

    if ( it2 != peptide_intensities_map.end() )
      intensities = it2->second;

    int ii;
    for (ii=0; ii < (int) intensities.size(); ii++) {
      line = line + std::string("\t") + ftos(intensities[ii], 2);
    }

    // retrieve pre-normalized intensities from map:
    intensities.clear();

    it2 = peptide_intensities_prenorm_map.find( pep_num );

    if ( it2 != peptide_intensities_prenorm_map.end() )
      intensities = it2->second;

    for (ii=0; ii < (int) intensities.size(); ii++) {
      line = line + std::string("\t") + ftos(intensities[ii], 2);
    }

    // retrieve rejected flag from std::map:
    std::string flag;
    std::map<int, std::string >::iterator it3;

    it3 = peptide_keptFlag_map.find ( pep_num );

    if ( it3 != peptide_keptFlag_map.end() ) {
      flag = it3->second;
    } else {
      // cerr << "No flags in std::map?" << endl;
    }

    line = line + std::string("\t") + flag;
    peptide_quantitation_lines.push_back( line );
  } // end iteration over peptides

  // Format protein information and store in quantitation file
  // expecting format: protein 2Xn_channels+1 of tabs r1\tr2\tr3\tr4\te1\te2\te3\te4
  // protein name plus tabs for place holders so that 
  // columns won't interfere with peptide ratios's and intensities
  // so user can do math in using Excell if they want
  line.erase();

  line = current_protein + std::string("\t");

  for (int i =1; i <= 2*nChannels; i++) {
    line = line + std::string("\t");
  }

  //tab placeholder for is_rejected column:
  line = line + std::string("\t");

  // get protein quantitation:
  int channel;
  for (channel=0; channel < nChannels; channel++) {
    if (norm_channel_ > 0) {
      line = line + std::string("\t") + ftos(protein_ratio_array_wrt_refchannel_[channel], 2);
    }
    else {
      line = line + std::string("\t") + ftos(protein_ratio_array_[channel], 2);
    }
  }

  // get protein errors:
  for (channel=0; channel < nChannels; channel++) {
    if (norm_channel_ > 0) {
      line = line + std::string("\t") + ftos(protein_error_array_wrt_refchannel_[channel], 2);
    }
    else {
      line = line + std::string("\t") + ftos(protein_se_array_[channel], 2);
    }
  }

  protein_quantitation_lines.push_back(line);

  // Write formatted lines to outfile

  //pQF->reset();

  int successful = 1;
  ofstream out;

  try {
    out.open( quantitation_file_name_, ios::app );
    if (!out) throw 100;

    // write protein lines to outfile:
    int i;
    for (i=0; i < (int) protein_quantitation_lines.size(); i++) {
      out << protein_quantitation_lines[i].c_str() << endl;
    }

    // write peptide lines to outfile:
    for (i=0; i < (int) peptide_quantitation_lines.size(); i++) {
      out << peptide_quantitation_lines[i].c_str() << endl;
    }
    successful = 0;

  } catch (int code) {
    if (code == 99)
      cerr << "Error: quantitation file name not set "<< endl;
    if (code == 100)
      cerr << "Error: could not open " << quantitation_file_name_ << "for appending" << endl;
  }
  out.close();
}
