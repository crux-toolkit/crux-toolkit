/*
Program       : ASAPRatioPvalueParser
Author        : Andrew Keller <akeller@systemsbiology.org>
                Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN rev       : $Id: ASAPRatioPvalueParser.cpp 8030 2020-02-18 19:25:48Z mhoopmann $

Overwrites ASAPRatio protein pvalue information into protXML

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


#include "ASAPRatioPvalueParser.h"
#include "Common/util.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "ASAPRatioPvalueParserTagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "Util/RACI/RACI.h"

ASAPRatioPvalueParser::ASAPRatioPvalueParser() : Parser("asapratio_pvalue") {
  pngfile_ = NULL;
  testMode_ = NULL;
}

ASAPRatioPvalueParser::ASAPRatioPvalueParser(const char* xmlfile, const char *testMode) : Parser("asapratio_pvalue") {
  use_seq_ratios_ = True;
  pngfile_ = NULL;
  testMode_ = testMode?strdup(testMode):NULL; // regression test stuff - bpratt Insilicos LLC, Nov 2005
  init(xmlfile);
}

ASAPRatioPvalueParser::ASAPRatioPvalueParser(const char* xmlfile, const char* pngfile, const char *testMode) : Parser("asapratio_pvalue") { 
  use_seq_ratios_ = True;
  pngfile_ = new char[strlen(pngfile)+1];
  strcpy(pngfile_, pngfile);
  testMode_ = testMode?strdup(testMode):NULL; // regression test stuff - bpratt Insilicos LLC, Nov 2005
  init(xmlfile);
}

ASAPRatioPvalueParser::~ASAPRatioPvalueParser() {
  if(pngfile_ != NULL)
    delete[] pngfile_;
  if(testMode_ != NULL)
    free(testMode_);
}

void ASAPRatioPvalueParser::stripWebRoot(char* filename) {
  char szWebserverRoot[1000];
  char* tmpFileName = new char[strlen(filename)+1];
  char tStr[1000] = "/";
  const char *pStr=getWebserverRoot();
  if (pStr==NULL) {
    printf("WARNING: Environment variable WEBSERVER_ROOT does not exist.\n\n");
  }
  else {
    strcpy(tStr, pStr);
  }
  pStr = tStr;
  strcpy(szWebserverRoot,pStr);

  int tmpLen=(int)strlen(szWebserverRoot);
  if (szWebserverRoot[tmpLen-1] == '/') {
    szWebserverRoot[tmpLen-1] = '\0';
  }
  unCygwinify(szWebserverRoot); // no effect in cygwin builds
  unCygwinify(filename); // no effect in cygwin builds
  pStr = strstri(filename, szWebserverRoot);

  if (pStr == NULL) {
    cout << "warning: outfile " << filename << " must begin with the WEBSERVER_ROOT path for use with a web browser." << endl;
  }
  else {
    pStr += sizeof(char)*strlen(szWebserverRoot);
    sprintf(tmpFileName, "%s", pStr);
    strcpy(filename, tmpFileName);
    delete[] tmpFileName;
  }
}

// in this case, must go through the data 2x
void ASAPRatioPvalueParser::parse(const char* xmlfile) {
  //open file and pass along
  //  int line_width = 10000;
  char *nextline = new char[line_width_];
  char* data = NULL;
  Tag* tag;

  //cout << "xmlfile: " << xmlfile << endl;
  char* modelfile = getModelFile(xmlfile);
  if(modelfile == NULL) {
    cout << "error: cannot derive modelfile from xmlfile: " << xmlfile << endl;
    exit(1);
  }

  RACI fin1(xmlfile); // can read gzipped xml
  if(! fin1) {
    cerr << "1: error opening " << xmlfile << endl;
    exit(1);
  }

  TagFilter* summary_filter = new TagFilter("analysis_summary", 1);
  summary_filter->enterRequiredAttributeVal("analysis", getName());
  TagFilter* result_filter = new TagFilter("analysis_result");
  result_filter->enterRequiredAttributeVal("analysis", getName());

  Tag* summary = NULL;
  Tag* summary_model_start = new Tag("ASAP_pvalue_analysis_model", True, False);
  Tag* summary_model_stop  = new Tag("ASAP_pvalue_analysis_model", False, True);
  Tag* result_start = new Tag("analysis_result", True, False);
  result_start->setAttributeValue("analysis", getName());
  //  result_start->setAttributeValue("time", time_);
  Tag* result_stop = new Tag("analysis_result", False, True);
  Tag* summary_start = new Tag("analysis_summary", True, False);
  summary_start->setAttributeValue("analysis", getName());
  summary_start->setAttributeValue("time", time_);
  summary_start->setAttributeValue("id", "1");
  Tag* summary_stop = new Tag("analysis_summary", False, True);

  // construct a tempfile name, possibly in the tmp dir if so configured
  std::string outfile = make_tmpfile_name(xmlfile);

  ofstream fout(outfile.c_str());
  if(! fout) {
    cout << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  //
  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  //
  Array<Tag*> test_tags;
  eTagListFilePurpose testType;
  char *testFileName = NULL;
  checkRegressionTestArgs(testMode_,testType);
  if (testType!=NO_TEST) {
    testFileName = constructTagListFilename(xmlfile, // input file
					    testMode_, // program args
					    "ASAPRatioPvalueParser", // program name
					    testType); // user info output
  }

#define RECORD(tag) {(tag)->write(fout);if (testType!=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}


  Array<double>* ratios = new Array<double>;

  double* altratios = NULL;

  while(fin1.getline(nextline, line_width_)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      if(
	 (use_seq_ratios_ && tag->isStart() && ! strcmp(tag->getName(), "ASAP_Seq") && 
	  ! strcmp(tag->getAttributeValue("include"), "1")) ||
	 (! use_seq_ratios_ && tag->isStart() && ! strcmp(tag->getName(), "ASAPRatio"))) {

	ratios->insertAtEnd(atof(tag->getAttributeValue("ratio_mean")));
	// attach std dev here as well (ratio[2]) ?

      }
      delete tag;

      data = strstr(data+1, "<");
    }

  }

  fin1.close();

  // now process this information
  altratios = new double[ratios->length()];
  for(int k = 0; k < ratios->length(); k++)
    altratios[k] = (double)((*ratios)[k]);

  pValueStrct params = getNormParams(altratios, ratios->length(), modelfile); // can be written in summary

  if(altratios != NULL)
    delete[] altratios;
  if(ratios != NULL)
    delete ratios;

  // for now only continue if pvalue analysis successful

  Normalization* norm = new Normalization(params);
  //cout << "norm results: " << norm->mean_ << " " << norm->stddev_ << " " << norm->merr_ << endl;
  if(norm->merr_ < 0.0) {
    cout << "no pvalues available for this dataset" << endl;
    fout.close();
    verified_unlink(outfile); // kill tempfile
    return; // done	
  }

  RACI fin2(xmlfile); // can read gzipped xml
  if(! fin2) {
    cout << "2: error opening " << xmlfile << endl;
    fout.close();
    verified_unlink(outfile); // kill tempfile
    exit(1);
  }

  double nextratio[2];
  double next_inv_ratio[2];
  double nextpval;
  char next[500];
  Tag* pvalue_tag = NULL;
  while(fin2.getline(nextline, line_width_)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      //tag->write(cout);
      // filter out entries below min prob, and exclude all previous ASAPRatio calculations

      if(tag->isStart() && ! strcmp(tag->getName(), "ASAPRatio")) {
	if(pvalue_tag != NULL) {
	  delete pvalue_tag;
	  pvalue_tag = NULL;
	}

	if(strcmp(tag->getAttributeValue("status"), "-1")) {

	  // get ratio and stdev
	  nextratio[0] = (double)atof(tag->getAttributeValue("ratio_mean"));
	  nextratio[1] = (double)atof(tag->getAttributeValue("ratio_standard_dev"));

	  //TODO: DDS should really model both h2l and l2h ratios with own adjustment models
	  //But for now we'll first invert, then adjust, then invert again
	  double tmp[2];
	  tmp[0] = (double)atof(tag->getAttributeValue("heavy2light_ratio_mean"));
	  tmp[1] = (double)atof(tag->getAttributeValue("heavy2light_ratio_standard_dev"));

	  next_inv_ratio[0] = 1 / tmp[0];
	  next_inv_ratio[1] = tmp[1]/(tmp[0]*tmp[0]);
	  nextpval = (double)norm->normalize(next_inv_ratio);

	  tmp[0] = next_inv_ratio[0];
	  tmp[1] = next_inv_ratio[1];

	  next_inv_ratio[0] = 1/ tmp[0];
	  next_inv_ratio[1] = tmp[1]/(tmp[0]*tmp[0]);

	  // compute modifications
	  nextpval = (double)norm->normalize(nextratio);

	  if(nextratio[0] >= 0.0) {
	    pvalue_tag = new Tag("ASAPRatio_pvalue", True, True);

	    sprintf(next, "%0.2f", nextratio[0]);
	    pvalue_tag->setAttributeValue("adj_ratio_mean", next);
	    sprintf(next, "%0.2f", nextratio[1]);
	    pvalue_tag->setAttributeValue("adj_ratio_standard_dev", next);

	    // heavy2light here...
	    if(nextratio[0] == -2.0) {
	      pvalue_tag->setAttributeValue("heavy2light_adj_ratio_mean", "-1.00");
	      pvalue_tag->setAttributeValue("heavy2light_adj_ratio_standard_dev", "0.00");
	    }
	    else if(nextratio[0] == -1.0) {
	      pvalue_tag->setAttributeValue("heavy2light_adj_ratio_mean", "0.00");
	      pvalue_tag->setAttributeValue("heavy2light_adj_ratio_standard_dev", "0.00");
	    }
	    else if(nextratio[0] == 0.0) {
	      pvalue_tag->setAttributeValue("heavy2light_adj_ratio_mean", "999.");
	      pvalue_tag->setAttributeValue("heavy2light_adj_ratio_standard_dev", "0.00");
	    }
	    else {
	      sprintf(next, "%0.2f", next_inv_ratio[0]);
	      pvalue_tag->setAttributeValue("heavy2light_adj_ratio_mean", next);
	      if (isnan(next_inv_ratio[1])) {
		pvalue_tag->setAttributeValue("heavy2light_adj_ratio_standard_dev", "NaN");
	      }
	      else {
		sprintf(next, "%0.2f",  next_inv_ratio[1]);
		pvalue_tag->setAttributeValue("heavy2light_adj_ratio_standard_dev", next);
	      }
	    }

	  }
	  if(pvalue_tag != NULL && nextpval >= 0.0) {
	    sprintf(next, "%0.2e", nextpval);
	    pvalue_tag->setAttributeValue("pvalue", next);
	    sprintf(next, "%0.8f", nextpval);
	    pvalue_tag->setAttributeValue("decimal_pvalue", next);
	  }
	} // if status is not -1
      }
      /*
      if(tag->isEnd() && ! strcmp(tag->getName(), "analysis_result") && pvalue_tag != NULL) {

	//      if(tag->isStart() && ! strcmp(tag->getName(), "peptide") && pvalue_tag != NULL) {
	result_start->write(fout);
	pvalue_tag->write(fout);
	result_stop->write(fout);
	delete pvalue_tag;
	pvalue_tag = NULL;
      }
      */
      if(! summary_filter->filter(tag) && ! result_filter->filter(tag))
	RECORD(tag); // write it

      if(tag->isEnd() && ! strcmp(tag->getName(), "analysis_result") && pvalue_tag != NULL) {

	//      if(tag->isStart() && ! strcmp(tag->getName(), "peptide") && pvalue_tag != NULL) {
	RECORD(result_start);
	RECORD(pvalue_tag);
	RECORD(result_stop);
	delete pvalue_tag;
	pvalue_tag = NULL;
      }

      /*
      if(summary != NULL && ! strcmp(tag->getName(), "ASAP_prot_analysis_summary")) {
	summary_start->write(fout);
	summary->write(fout); // do this now after xpress
	summary_stop->write(fout);
	delete summary;
	summary = NULL;
      }
      */
      if(tag->isEnd() && ! strcmp(tag->getName(), "protein_summary_header")) {
	summary = new Tag("ASAP_pvalue_analysis_summary", True, True);
	// get time info
	char next[20];
	// here set the pvalue specific summary info
	sprintf(next, "%0.3f", norm->mean_);
	summary->setAttributeValue("background_ratio_mean", next);
	sprintf(next, "%0.3f", norm->stddev_);
	summary->setAttributeValue("background_ratio_stdev", next);
	sprintf(next, "%0.3f", norm->merr_);
	summary->setAttributeValue("background_fitting_error", next);

	// do this before we change modelfile in the else clause below
	char tmpPvalFile[strlen(modelfile)+11];
	sprintf(tmpPvalFile, "%s.model.data", modelfile);

	if(pngfile_ != NULL) {
	  stripWebRoot(pngfile_);
	  summary->setAttributeValue("analysis_distribution_file", pngfile_);
	  summary->setAttributeValue("full_analysis_distribution_file", modelfile); // only nec if differs from above
	}
	else {
	  summary->setAttributeValue("full_analysis_distribution_file", modelfile);
	  stripWebRoot(modelfile);
	  summary->setAttributeValue("analysis_distribution_file", modelfile);
	}

	RECORD(summary_start);
	RECORD(summary); // do this now after xpress

	FILE *fin;
	char* pvalline = new char[line_width_];

	//	cout << "debug: attempting to read " << tmpPvalFile  << " for Pvalue model..." << endl;

	// Write model data points?
	if ((fin = fopen(tmpPvalFile,"r"))!=NULL) {
	  RECORD(summary_model_start);
	  while(fgets(pvalline, line_width_, fin)) {
	    fout << pvalline;
	  }
	  RECORD(summary_model_stop);
	  fclose(fin);
	  verified_unlink(tmpPvalFile);
	}
	else {
	  cout << "warning: could not open pvalue model file " << tmpPvalFile  << ". Pvalue model will not be included in prot.xml file." << endl;
	}

	RECORD(summary_stop);
	delete[] pvalline;
	delete summary;
	summary = NULL;
      }

      delete tag;
      data = strstr(data+1, "<");
    } // while

  }
  fin2.close();
  fout.close();
  if(modelfile != NULL)
    delete[] modelfile;
  if(norm != NULL)
    delete norm;

  //return;
  if(! overwrite(xmlfile, outfile.c_str(), "</protein_summary>")) {
    cout << "error: no ASAPRatioPvalue data written to file " << xmlfile << endl;
  }

  if (testType!=NO_TEST) {
    //
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    //
    ASAPRatioPvalueParserTagListComparator("ASAPRatioPvalueParser",testType,test_tags,testFileName);
    delete[] testFileName;
    for(int k = test_tags.length(); k--;) {
      delete test_tags[k];
    }
  }

  delete summary_filter;
  delete result_filter;
  delete summary;
  delete summary_model_start;
  delete summary_model_stop;
  delete result_start;
  delete result_stop;
  delete summary_start;
  delete summary_stop;
  delete data;
  delete[] nextline;
}


void ASAPRatioPvalueParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "protein")){
    if(tag->isStart()) {
      //tag->print();
      filter_ = True;
    }else{
      filter_memory_ = True;
    }
  }

}

char* ASAPRatioPvalueParser::getModelFile(const char* xmlfile) {
  char* output = NULL;
  char suff[] = "-pval.png";
  const char* result = hasValidProtXMLFilenameExt(xmlfile);
  int dir_size = 5000;
  char directory[5000];

  if(result != NULL) {
    if(isAbsolutePath(xmlfile)) { // full path
      output = new char[strlen(xmlfile)-strlen(result)+strlen(suff)+1];
      strncpy(output, xmlfile, strlen(xmlfile)-strlen(result));
      output[strlen(xmlfile)-strlen(result)] = 0;
      strcat(output, suff);
    }
    else {
      safepath_getcwd(directory, dir_size);
      output = new char[strlen(directory)+1+strlen(xmlfile)-strlen(result)+strlen(suff)+1];
      strcpy(output, directory);
      strcat(output, "/");
      strncat(output, xmlfile, strlen(xmlfile)-strlen(result));
      output[strlen(directory)+1+strlen(xmlfile)-strlen(result)] = 0;
      strcat(output, suff);
    }
  }

  return output;
}


void ASAPRatioPvalueParser::writeProteinRatio(ostream& os, proDataStrct* pro_ratio) {
  Tag* protag = new Tag("ASAPRatio", True, False);
  char next[50];

  double ratiomean = (double)pro_ratio->ratio[0];
  if(ratiomean == -2.0) 
    sprintf(next, "%0.2f", -1.0);
  else if(ratiomean == -1.0) 
    strcpy(next, "9999.");
  else 
    sprintf(next, "%0.2f", ratiomean);

  protag->setAttributeValue("ratio_mean", next);
  sprintf(next, "%0.2f", (double)pro_ratio->ratio[1]);
  protag->setAttributeValue("ratio_standard_dev", next);
  sprintf(next, "%d", pro_ratio->dataNum);
  protag->setAttributeValue("ratio_number_peptides", next);
  sprintf(next, "%d", pro_ratio->indx);
  protag->setAttributeValue("status", next);
  protag->write(os);
  // now the protein stuff
  delete protag;

  for(int k = 0; k < pro_ratio->dataNum; k++) {
    Tag* seqtag = new Tag("ASAP_Seq", True, False);
    sprintf(next, "%d", pro_ratio->sequences[k].indx);
    seqtag->setAttributeValue("status", next);
    sprintf(next, "%d", pro_ratio->dataCnts[k]);
    seqtag->setAttributeValue("include", next);
    sprintf(next, "%d", pro_ratio->sequences[k].dataNum);
    seqtag->setAttributeValue("datanum", next);
    sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].ratio[0]);
    seqtag->setAttributeValue("ratio_mean", next);
    sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].ratio[1]);
    seqtag->setAttributeValue("ratio_standard_dev", next);
    sprintf(next, "%0.2f", pro_ratio->sequences[k].weight);
    seqtag->setAttributeValue("weight", next);
    seqtag->write(os);
    delete seqtag;

    pro_ratio->sequences[k].sort_for_output(); // place in a logical output order
    for(int j = 0; j < pro_ratio->sequences[k].dataNum; j++) {
      Tag* peaktag = new Tag("ASAP_Peak", True, False);
      sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].indx);
      peaktag->setAttributeValue("status", next);
      sprintf(next, "%d", pro_ratio->sequences[k].dataCnts[j]);
      peaktag->setAttributeValue("include", next);
      sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].dataNum);
      peaktag->setAttributeValue("datanum", next);
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].ratio[0]);
      peaktag->setAttributeValue("ratio_mean", next);
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].ratio[1]);
      peaktag->setAttributeValue("ratio_standard_dev", next);
      //cout << "wt: " << pro_ratio->sequences[k].peaks[j].weight << endl;
      sprintf(next, "%0.1f", (double)pro_ratio->sequences[k].peaks[j].weight);
      peaktag->setAttributeValue("weight", next);

      //DDS:
      //sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].bofIndx);
      sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].msms_run_idx);
      peaktag->setAttributeValue("peptide_binary_ind", next);

      peaktag->write(os);
      delete peaktag;

      for(int i = 0; i < pro_ratio->sequences[k].peaks[j].dataNum; i++) {
	Tag* dtatag = new Tag("ASAP_Dta", True, True);
	sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].dataIndx[i]);
	dtatag->setAttributeValue("peptide_index", next);
	sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].dataCnts[i]);
	dtatag->setAttributeValue("include", next);

	dtatag->write(os);
	delete dtatag;
      } // next dta

      // close it
      peaktag = new Tag("ASAP_Peak", False, True);
      peaktag->write(os);
      delete peaktag;
    } // next peak

    // close it
    seqtag = new Tag("ASAP_Seq", False, True);
    seqtag->write(os);
    delete seqtag;

  } // next seq

  protag = new Tag("ASAPRatio", False, True);
  protag->write(os);
  delete protag;
}

// This function evaluates the normalization parameters.
pValueStrct ASAPRatioPvalueParser::getNormParams(double *ratios, int dataNum, const char *pngFile)
{
  PvalueModel* pvalue = new PvalueModel(ratios, dataNum, pngFile, NULL!=testMode_);
  pValueStrct result = pvalue->getParams();
  delete pvalue;
  return result;
}
