/*
Program       : ASAPCGIDisplayParser
Authors       : Andrew Keller <akeller@systemsbiology.org>
                Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPCGIDisplayParser.cpp 8030 2020-02-18 19:25:48Z mhoopmann $

Displays ASAPRatio protein information from ProteinProphet XML

Copyright (C) 2003 Andrew Keller, Xiao-jun Li

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

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/


#include "ASAPCGIDisplayParser.h"
#include "Util/RACI/RACI.h"

// lookup by name:
ASAPCGIDisplayParser::ASAPCGIDisplayParser(const char* xmlfile, const char* protein) : Parser(NULL) { 
  heavy2light_ = False;
  protein_ = new char[strlen(protein)+1];
  strcpy(protein_, protein);
  inputfiles_ = NULL;
  inputfiles_array_ = NULL;
  pro_ratio_ = NULL;
  init(xmlfile);
}
// lookup by name, alos store index
ASAPCGIDisplayParser::ASAPCGIDisplayParser(const char* xmlfile, const char* protein, unsigned int index) : Parser(NULL) { 
  heavy2light_ = False;
  protein_ = new char[strlen(protein)+1];
  strcpy(protein_, protein);
  inputfiles_ = NULL;
  inputfiles_array_ = NULL;
  pro_ratio_ = NULL;
  index_ = index;
  init(xmlfile);
}
// lookup by index:
ASAPCGIDisplayParser::ASAPCGIDisplayParser(const char* xmlfile, unsigned int index) : Parser(NULL) {
  heavy2light_ = False;
  protein_ = NULL;
  inputfiles_ = NULL;
  inputfiles_array_ = NULL;
  pro_ratio_ = NULL;
  index_ = index;
  init(xmlfile);
}

void ASAPCGIDisplayParser::parse(const char* xmlfile) {
  //open file and pass along
  //  int line_width = 10000;
  char *nextline = new char[line_width_];
  char* data = NULL;
  Tag* tag;
  Boolean heavy2light = False;

  Boolean filterParams = False;
  Boolean quantHighBG = False;
  Boolean zeroBG = False;
  double mzBound = 0.5;
  bool wavelet = false;
  norm_ = NULL;
  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cout << "ASAPCGIDisplayParser: error opening " << xmlfile << endl;
    exit(1);
  }

  Array<Tag*>* tags = NULL;

  TagFilter* summary_filter = new TagFilter("ASAP_analysis_summary", 1);
  TagFilter* ratio_filter = new TagFilter("ASAPCGIDisplay");

  Array<char*>* inputfiles = new Array<char*>;

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);
  //cerr << "writing data to " << outfile << endl;
  ofstream fout(outfile.c_str());
  if(! fout) {
    cout << "cannot write output to file " << outfile << endl;
    exit(1);
  }
  pval_link_ = NULL;
  while(fin.getline(nextline, line_width_)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);

      if (tag->isStart() && ! strcmp(tag->getName(), "analysis_summary") 
	  && ! strcmp(tag->getAttributeValue("analysis"), "asapratio")) {
	filterParams = True;
      }

      if (filterParams && tag->isEnd() && ! strcmp(tag->getName(), "analysis_summary")) {
	filterParams = False;
      }
      if (filterParams && ! strcmp(tag->getName(), "parameter")) {
	if (! strcmp(tag->getAttributeValue("name"), "quantHighBG") && 
	    ! strcmp(tag->getAttributeValue("value"), "True")) {
	  quantHighBG = True;
	}
	else if (! strcmp(tag->getAttributeValue("name"), "zeroBG") && 
		 ! strcmp(tag->getAttributeValue("value"), "True")) {
	  zeroBG = True;
	}
	else if (! strcmp(tag->getAttributeValue("name"), "wavelet") && 
		 ! strcmp(tag->getAttributeValue("value"), "True")) {
	  wavelet = true;
	}
	else if (! strcmp(tag->getAttributeValue("name"), "mzBound")) {
	  mzBound = atof(tag->getAttributeValue("value"));
	  if (mzBound <= 0 || mzBound >= 1) {
	    mzBound = 0.5;
	  }
	}
      }

      if(! strcmp(tag->getName(), "ASAP_pvalue_analysis_summary") && tag->isStart() ) {
	// get the parameters
	pValueStrct params;
	params.mean = atof(tag->getAttributeValue("background_ratio_mean"));
	params.stddev = atof(tag->getAttributeValue("background_ratio_stdev"));
	params.merr = atof(tag->getAttributeValue("background_fitting_error"));
	const char* tmp = tag->getAttributeValue("analysis_distribution_file");
	pval_link_ = new char[strlen(tmp)+1];
	sprintf(pval_link_, "%s", tmp);
	norm_ = new Normalization(params);
      }

      setFilter(tag);
      //tag->write(cout);

      // parse input file names
      if(tag->isStart() && ! strcmp(tag->getName(), "protein_summary_header")) {
	const char* files = tag->getAttributeValue("source_files_alt");
	// parse through
	int i = 0;
	int last_i = 0;
	char* nextfile;
	while(files[i]) {
	  i++;
	  if(files[i] == '+' || !files[i]) {
	    nextfile = new char[i - last_i + 1];
	    for(int z = last_i; z < i; z++)
	      nextfile[z-last_i] = files[z];
	    nextfile[i - last_i] = 0;
	    inputfiles->insertAtEnd(nextfile);
	    last_i = i+1;
	  }
	} // while

	// now call to set
	inputfiles_ = new char* [inputfiles->length()+1];
	for(int k = 0; k < inputfiles->length(); k++)
	  inputfiles_[k] = (*inputfiles)[k];
	inputfiles_[inputfiles->length()] = strdup("");

	inputfiles_array_ = inputfiles;
      } // if protein summary header

      if(filter_) { 
	//tag->write(cout);
	if(!protein_ && !strcmp(tag->getName(), "protein")) {
	  protein_ = strdup(tag->getAttributeValue("protein_name"));
	}
	if(! strcmp(tag->getName(), "ASAPRatio")) {
	  if(tag->isStart()) {
	    tags = new Array<Tag*>;
	    tags->insertAtEnd(tag);
	  }
	  else { // end tag
	    if(tags == NULL) {
	      cerr << "null tags" << endl;
	      exit(1);
	    }

	    pro_ratio_ = getProDataStruct(tags);
	    if(pro_ratio_ == NULL)
	      cout << "NULL PRORATIO" << endl;

	    for(int k = 0; k < tags->length(); k++)
	      if((*tags)[k] != NULL)
		delete (*tags)[k];

	    if(tags != NULL) {
	      delete tags;
	      tags = NULL;
	    }

	    fin.close();
	    return;
	  }
	} // if start of ratio
	else {
	  if(tags != NULL && tag->isStart())
	    tags->insertAtEnd(tag);
	  else delete tag;
	}
      } // if filter and ready....
      else {
	delete tag;
      }
      data = strstr(data+1, "<");
    }
  }

  fin.close();
  delete[] nextline;
}

proDataStrct* ASAPCGIDisplayParser::getProDataStrct() {
  return pro_ratio_;
}

Normalization* ASAPCGIDisplayParser::getNormalized() {
  return norm_;
}

char* ASAPCGIDisplayParser::getPvalLink() {
  return pval_link_;
}

Array<char*>* ASAPCGIDisplayParser::getInputFilesArray() {
  return inputfiles_array_;
}

char** ASAPCGIDisplayParser::getInputFiles() {
  return inputfiles_;
}

// from beginning to end of each protein
void ASAPCGIDisplayParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  const char* stro;

  if(protein_) {  //lookup by protein name
    if(tag->isStart() && ! strcmp(tag->getName(), "protein") && ! strcmp(tag->getAttributeValue("protein_name"), protein_)) {
      filter_ = True;
    }
    else if(filter_ && ! strcmp(tag->getName(), "protein") && ! tag->isStart())
      filter_memory_ = True;
  }
  else {  //lookup by index only if not found by name, since looking by index doesn't work for multi protein groups
    if(tag->isStart() && ! strcmp(tag->getName(), "protein_group") && ( stro = tag->getAttributeValue("group_number"))) {
      if((unsigned)atoi(stro) == index_) {
	filter_ = True;
      }
      else if(filter_)
	filter_memory_ = True;
    }
  }
}

// This function frees a proDataStrct.
void ASAPCGIDisplayParser::freeProDataStrct(proDataStrct* data) {
  int i, j;

  // sequences
  for (i = 0; i < data->dataNum; ++i) {
    // peaks
    for (j = 0; j < data->sequences[i].dataNum; ++j) {
      free(data->sequences[i].peaks[j].dataIndx);
      free(data->sequences[i].peaks[j].dataCnts);
    }
    free(data->sequences[i].peaks);

    // dataCnts
    free(data->sequences[i].dataCnts);
  } // for (i = 0; i < data.dataNum; ++i) {
  free(data->sequences);

  // dataCnts
  free(data->dataCnts);

  return;
}

proDataStrct* ASAPCGIDisplayParser::getProDataStruct(Array<Tag*>* tags) {
  if(tags == NULL || tags->length() < 1) {
    return NULL;
  }
  proDataStrct* data = (proDataStrct *)calloc(1, sizeof(proDataStrct));

  // first one should be ASAPCGIDisplay
  int index = 0;

  if(! (*tags)[index]->isStart() || strcmp((*tags)[index]->getName(), "ASAPRatio")) {
    cout << "first tag is not ASAPRatio" << endl;
    exit(1);
  }
  //(*tags)[index]->write(cout);
  data->ratio[0] = (double)atof((*tags)[index]->getAttributeValue("ratio_mean"));
  data->inv_ratio[0] = (double)atof((*tags)[index]->getAttributeValue("heavy2light_ratio_mean"));

  if(data->ratio[0] == -1.0)
    data->ratio[0] = -2.0;
  else if(data->ratio[0] >= 9999.0)
    data->ratio[0] = -1.0;

  if(data->inv_ratio[0] == -1.0)
    data->inv_ratio[0] = -2.0;
  else if(data->ratio[0] >= 9999.0)
    data->inv_ratio[0] = -1.0;

  data->ratio[1] = (double)atof((*tags)[index]->getAttributeValue("ratio_standard_dev"));
  data->inv_ratio[1] = (double)atof((*tags)[index]->getAttributeValue("heavy2light_ratio_standard_dev"));

  data->dataNum = atoi((*tags)[index]->getAttributeValue("ratio_number_peptides"));
  data->indx = atoi((*tags)[index]->getAttributeValue("status"));
  // sequences
  data->sequences = (seqDataStrct *)
    calloc(data->dataNum, sizeof(seqDataStrct));
  // dataCnts
  data->dataCnts = (int *) calloc(data->dataNum, sizeof(int));
  index++;

  // now go through each seq
  for(int k = 0; k < data->dataNum; k++) {
    //    cout << "working on seq " << (k+1) << " of " << data->dataNum << endl;
    index = setProSeqStruct(data, k, tags, index);
  } //  next seq
  if(index != (int)tags->length()) {
    cerr << "error, " << index << " != " << tags->length() << endl;
    exit(1);
  }
  return data;
}

int ASAPCGIDisplayParser::setProSeqStruct(proDataStrct* data, int seq_ind, Array<Tag*>* tags, int tag_ind) {
  // first one should be ASAPCGIDisplay
  int index = tag_ind;

  if(strcmp((*tags)[index]->getName(), "ASAP_Seq")) {
    cout << "first tag is not ASAP_Seq" << endl;
    (*tags)[tag_ind]->write(cerr);
    exit(1);
  }
  seqDataStrct* seq = data->sequences + seq_ind;
  seq->ratio[0] = (double)atof((*tags)[index]->getAttributeValue("ratio_mean"));
  seq->inv_ratio[0] = (double)atof((*tags)[index]->getAttributeValue("heavy2light_ratio_mean"));
  // revert changes
  seq->ratio[1] = (double)atof((*tags)[index]->getAttributeValue("ratio_standard_dev"));
  seq->inv_ratio[1] = (double)atof((*tags)[index]->getAttributeValue("heavy2light_ratio_standard_dev"));
  seq->dataNum = atoi((*tags)[index]->getAttributeValue("datanum"));
  seq->indx = atoi((*tags)[index]->getAttributeValue("status"));
  seq->weight = atof((*tags)[index]->getAttributeValue("weight"));
  if((*tags)[index]->getAttributeValue("light_sequence") != NULL) {
    strcpy(seq->lightSeq, (*tags)[index]->getAttributeValue("light_sequence"));
  }

  data->dataCnts[seq_ind] = atoi((*tags)[index]->getAttributeValue("include"));
  // peaks
  seq->peaks = (dataStrct *)
    calloc(seq->dataNum, sizeof(dataStrct));
  // dataCnts
  seq->dataCnts = (int *) calloc(seq->dataNum, sizeof(int));
  index++;
  // now go through each peak
  for(int k = 0; k < seq->dataNum; k++) {
    //cout << "working on peak " << (k+1) << " of " << seq->dataNum << endl;
    index = setProPeakStruct(seq, k, tags, index);
  } //  next seq
  //cout << "index: " << index << endl;
  return index;
}

int ASAPCGIDisplayParser::setProPeakStruct(seqDataStrct* seq, int peak_ind, Array<Tag*>* tags, int tag_ind) {
  // first one should be ASAPCGIDisplay
  int index = tag_ind;

  if(strcmp((*tags)[index]->getName(), "ASAP_Peak")) {
    cout << "first tag is not ASAP_Seq for index " << index << endl;
    (*tags)[tag_ind]->write(cerr);
    exit(1);
  }
  //cout << "peak name: " << (*tags)[index]->getName() << "<p>" << endl;
  //cout << "Peak num: " << peak_ind << " of tot: " << seq->dataNum << "<p>" << endl;

  dataStrct* peak = seq->peaks + peak_ind;
  peak->ratio[0] = (double)atof((*tags)[index]->getAttributeValue("ratio_mean"));
  peak->inv_ratio[0] = (double)atof((*tags)[index]->getAttributeValue("heavy2light_ratio_mean"));
  //cout << "done1" << endl;
  peak->ratio[1] = (double)atof((*tags)[index]->getAttributeValue("ratio_standard_dev"));
  peak->inv_ratio[1] = (double)atof((*tags)[index]->getAttributeValue("heavy2light_ratio_standard_dev"));
  //cout << "done2" << endl;
  peak->dataNum = atoi((*tags)[index]->getAttributeValue("datanum"));
  peak->bofIndx = atoi((*tags)[index]->getAttributeValue("peptide_binary_ind"));

  //DDS:
  //peak->msms_run_idx = atoi((*tags)[index]->getAttributeValue("peptide_binary_ind"));

  peak->indx = atoi((*tags)[index]->getAttributeValue("status"));
  peak->weight = (double)atof((*tags)[index]->getAttributeValue("weight"));
  seq->dataCnts[peak_ind] = atoi((*tags)[index]->getAttributeValue("include"));

  // peaks
  peak->dataIndx = (int *)calloc(peak->dataNum, sizeof(int));

  // dataCnts
  peak->dataCnts = (int *) calloc(peak->dataNum, sizeof(int));

  index++;
  // now go through each dta
  for(int k = 0; k < peak->dataNum; k++) {
    //cout << "working on data " << (k+1) << " out of " << peak->dataNum << endl;
    index = setProDtaStruct(peak, k, tags, index);
  } //  next seq
  //cout << "index: " << index << endl;
  return index;
}

int ASAPCGIDisplayParser::setProDtaStruct(dataStrct* peak, int dta_ind, Array<Tag*>* tags, int tag_ind) {
  // first one should be ASAPCGIDisplay
  int index = tag_ind;

  if(strcmp((*tags)[index]->getName(), "ASAP_Dta")) {
    (*tags)[tag_ind]->write(cerr);
    exit(1);
  }
  if(index >= (int)tags->length()) {
    cout << "error: tags length " << tags->length() << " less than index " << index << endl;
    exit(1);
  }

  peak->dataCnts[dta_ind] = atoi((*tags)[index]->getAttributeValue("include"));
  peak->dataIndx[dta_ind] = atoi((*tags)[index]->getAttributeValue("peptide_index"));
  index++; // done

  return index;
}
