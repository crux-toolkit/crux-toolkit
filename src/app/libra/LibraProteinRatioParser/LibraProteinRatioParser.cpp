/*
  Program       : LibraProteinRatioParser
  Author        : Andrew Keller <akeller@systemsbiology.org>
                  *Jimmy Eng (jeng@systemsbiology.org>
  Date          : 11.27.02
  SVN info      : $Id: LibraProteinRatioParser.cpp 8873 2023-02-27 08:14:30Z real_procopio $

  Computes LIBRA ratios and errors for proteins, then overwrites
  that information onto ProteinProphet XML

  Copyright (C) 2003 Andrew Keller, P.Pedrioli

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

#include "LibraProteinRatioParser.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "Util/RACI/RACI.h" // for reading gzipped files with efficient seeks

LibraProteinRatioParser::LibraProteinRatioParser(const char * xmlfile, const char * conditionFileName, char *testArg) : Parser("libra") { 

  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  testMode_ = testArg?strdup(testArg):NULL;

  // read condition file to get arguments for LibraGroupPeptideParsers
  pLibraConditionHandler2 = new LibraConditionHandler2();
  pLibraConditionHandler2->setConditionFileName(conditionFileName);
  pLibraConditionHandler2->readFile();

  int tmp_norm_channel = pLibraConditionHandler2->getNormalizationChannel();

  if (tmp_norm_channel > 0)
    norm_channel_ = tmp_norm_channel;
  else
    norm_channel_ = 0;

  norm_channel_ = pLibraConditionHandler2->getNormalizationChannel();

  quantitationFileName = pLibraConditionHandler2->getQuantitationFileName();
  if (quantitationFileName.size() < 1)
    quantitationFileName = "quantitation.tsv";

  initializeQuantitationFile(pLibraConditionHandler2->getNumberOfReagentLines());

  // if value returned is -99., condition file did not specify a value, so
  // will carry through no use of minimum threshhold to LibraGroupPeptideParser
  minimumThreshholdIntensity = pLibraConditionHandler2->getMinimumThreshholdIntensity();

  libraGroupPeptideParser_ = NULL;
  summary_tags_ = NULL;
  result_tags_ = NULL;
  input_xmlfiles_ = new Array<const char *>;
  init(xmlfile);
}


LibraProteinRatioParser::~LibraProteinRatioParser() {
  if (result_tags_ != NULL)
    for(int k = 0; k < result_tags_->length(); k++)
      if((*result_tags_)[k] != NULL)
	delete (*result_tags_)[k];
  delete result_tags_;

  if(input_xmlfiles_ != NULL)
    for(int k = 0; k < input_xmlfiles_->length(); k++)
      delete[] (*input_xmlfiles_)[k];
  delete input_xmlfiles_;

  if(libraGroupPeptideParser_ != NULL)
    delete libraGroupPeptideParser_;

  if (pLibraConditionHandler2 != NULL)
    delete pLibraConditionHandler2;

  free(testMode_);
}


void LibraProteinRatioParser::initializeQuantitationFile(int nChannels) {
  // header:
  // protein peptide\tnr1\tnr2\tnr3t\nr4\tin1\in2\tin3\tin4\tis_rejected\tprot...
  std::string header = "#protein\tpeptide";

  int i;
  for (i =0; i < nChannels; i++)
    header = header + "\tpepratio" + itos(i);
  for (i =0; i < nChannels; i++)
    header = header + "\tinten" + itos(i);

  header = header + std::string("\tkept?");

  // protein quantitation and errors
  for (i =0; i < nChannels; i++)
    header = header + "\tprotratio" + itos(i);
  for (i =0; i < nChannels; i++)
    header = header + "\tproterr" + itos(i);

  int successful = 1;

  std::ofstream out (quantitationFileName.c_str());
  try {
    if (!out) throw 100;
    out << header << endl;
    successful = 0;

  } catch (int code) {
    if (code == 100 )
      std::cerr << "Error: could not open " << quantitationFileName << std::endl;
  }

  out.close();
}


/**
 * Parse interact-prot.xml file and add Libra protein quantitation.  
 * ProteinProphet just wrote the interact-prot.xml file.
 * @param xmlfile protXML file such as interact-prot.xml
 */
void LibraProteinRatioParser::parse(const char * xmlfile) {
  //open file and pass along
  char *nextline = new char[line_width_];
  char* data = NULL;
  Tag* tag;
  Boolean heavy2light = False;

  RACI fin(xmlfile); // can read gzipped XML
  if(! fin) {
    cerr << "LibraProteinRatioParser: error opening " << xmlfile << endl;
    exit(1);
  } else {
    cout << "Calculating and updating protein ratios in " << xmlfile << "..." << endl;
  }

  //
  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  //
  Array<Tag*> test_tags;
  eTagListFilePurpose testType;
  char *testFileName=NULL;
  checkRegressionTestArgs(testMode_,testType);
  if (testType!=NO_TEST) {
    testFileName = constructTagListFilename(xmlfile, // input file
					    testMode_, // program args
					    "libraProteinRatioParser", // program name
					    testType); // user info output
  }

#define RECORD(tag) {(tag)->write(fout);if (testType!=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}

  double MIN_PEP_WT = 0.5;
  double MIN_PEP_PROB = 0.5;
  Array<Tag*>* tags = NULL;
  Array<const char*>* peps = NULL;
  const char* prot_name = NULL;
  double next_prot_prob = 2.0;
  double MIN_PROB = 0.2;
  Boolean done = False;

  TagFilter* summary_filter = new TagFilter("analysis_summary", 1);
  summary_filter->enterRequiredAttributeVal("analysis", getName());

  TagFilter* ratio_filter = new TagFilter("analysis_result");
  ratio_filter->enterRequiredAttributeVal("analysis", getName());

  Tag* result_start = new Tag("analysis_result", True, False);
  result_start->setAttributeValue("analysis", getName());

  Tag* result_stop = new Tag("analysis_result", False, True);

  Tag* summary_start = new Tag("analysis_summary", True, False);
  summary_start->setAttributeValue("analysis", getName());
  summary_start->setAttributeValue("time", time_);
  summary_start->setAttributeValue("id", "1");

  Tag* summary_stop = new Tag("analysis_summary", False, True);

  char** altpeps = NULL;

  // construct a tempfile name, possibly in the tmp dir if so configured
  std::string outfile = make_tmpfile_name(xmlfile);

  ofstream fout(outfile.c_str());
  if(! fout) {
    cout << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  while(fin.getline(nextline, line_width_)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      //tag->write(cout);

      setFilter(tag);
      Boolean stored = False;
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
	    input_xmlfiles_->insertAtEnd(nextfile);
	    last_i = i+1;
	  }
	} // while

	// now call to set
      } // if protein summary header

      if(input_xmlfiles_ != NULL && input_xmlfiles_->length() > 0 && summary_tags_ == NULL) {

	// grab LIBRA channel settings from pepXML
	libraGroupPeptideParser_ = new LibraGroupPeptideParser(input_xmlfiles_, quantitationFileName.c_str(), MIN_PEP_PROB, norm_channel_, minimumThreshholdIntensity);
	if(libraGroupPeptideParser_ != NULL) {
	  summary_tags_ = libraGroupPeptideParser_->getProtXMLSummaryTags(MIN_PEP_PROB, MIN_PEP_WT, MIN_PROB);
	  num_channels_ = libraGroupPeptideParser_->getNumChannels();
	} else {
	  cout << "error: could not get summary tags" << endl;
	  exit(1);
	}

      }

      if(! done && tag->isStart() && ! strcmp(tag->getName(), "protein_group")) {
	const char* next = tag->getAttributeValue("probability");
	done = (next == NULL || atof(next) < MIN_PROB);

	next = tag->getAttributeValue("group_number");
	//	cout << "(" << next << ") ";
	if (atoi(next) % 10 == 0)
	  cout << "..." << next;

	if (atoi(next) % 100 == 0)
	  cout << endl;

	if (done)
	  cout << "...done! " << endl;

	flush(cout);
      }

      if(! done && tag->isStart() && ! strcmp(tag->getName(), "protein")) {
        prot_name = tag->getAttributeValue("protein_name");
	const char* next = tag->getAttributeValue("probability");
	next_prot_prob = next == NULL ? 0.0 : atof(next);
      }

      // filter out entries below min prob, and exclude all previous Libra calculations
      if(! done && filter_ && next_prot_prob >= MIN_PROB) {

	if(tags == NULL)
	  tags = new Array<Tag*>;

	if(! ratio_filter->filter(tag)) {
	  tags->insertAtEnd(tag);
	  stored = True;
	}

	if(peps == NULL)
	  peps = new Array<const char*>;

	if(tag->isStart() && ! strcmp(tag->getName(), "peptide")) {

	  // check that weight and prob are above minimum.....
	  if(atof(tag->getAttributeValue("weight")) > MIN_PEP_WT && atof(tag->getAttributeValue("nsp_adjusted_probability")) >= MIN_PEP_PROB)
	    enterUnique(peps, tag->getAttributeValue("peptide_sequence")); // add peptide to list for current protein
	}
	else if(filter_memory_ && tags != NULL && peps != NULL) {

	  // compute protein ratio only if have enough peps
	  if(peps->length() > 0)
	    getRatio(prot_name, peps, MIN_PEP_PROB);

	  // add LIBRA result to protein info
	  int k;
	  for(k = 0; k < tags->length(); k++) {
	    RECORD((*tags)[k]); //print();
	    if(k == 0 && result_tags_ != NULL) {
	      RECORD(result_start);

	      for(int j = 0; j < result_tags_->length(); j++) {
		if((*result_tags_)[j] != NULL) {
		  RECORD((*result_tags_)[j]);
		  delete (*result_tags_)[j];
		}
	      }
	      RECORD(result_stop);

	      delete result_tags_;
	      result_tags_ = NULL;
	    }
	  }

	  for(k = 0; k < tags->length(); k++)
	    if((*tags)[k] != NULL)
	      delete (*tags)[k];

	  if(tags != NULL) {
	    delete tags;
	    tags = NULL;
	  }

	  if(peps != NULL) {
	    delete peps;
	    peps = NULL;
	  }
	}

      } else {
	// write LIBRA summary info here
	if(! summary_filter->filter(tag)) {
	  RECORD(tag); //print();
	}
	if(tag->isEnd() && ! strcmp(tag->getName(), "protein_summary_header")) {
	  if(summary_tags_ != NULL) {
	    RECORD(summary_start);
	    for(int k = 0; k < summary_tags_->length(); k++) {
	      RECORD((*summary_tags_)[k]);
	      delete (*summary_tags_)[k];
	    }
	    delete summary_tags_;
	    RECORD(summary_stop);
	  }
	}
      }
      if(! stored)
	delete tag;

      data = strstr(data+1, "<");
    }
  }

  fin.close();
  fout.close();

  if(! overwrite(xmlfile, outfile.c_str(), "</protein_summary>")) {
    cout << "error: no " << getName() << " data written to file " << xmlfile << endl;
  }

  if (testType!=NO_TEST) {
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    TagListComparator("LibraProteinRatioParser",testType,test_tags,testFileName);
    delete[] testFileName;
    for(int k = test_tags.length(); k--;) {
      delete test_tags[k];
    }
  }

  delete summary_filter;
  delete ratio_filter;
  delete result_start;
  delete result_stop;
  delete summary_start;
  delete summary_stop;
  delete[] nextline;
}

void LibraProteinRatioParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "protein")) {
    if(tag->isStart()) {
      //tag->print();
      filter_ = True;
    }else{
      filter_memory_ = True;
    }
  }
}


void LibraProteinRatioParser::enterUnique(Array<const char*>* uniques, const char* next) {
  for(int k = 0; k < uniques->length(); k++)
    if(! strcmp((*uniques)[k], next))
      return;
  uniques->insertAtEnd(next);
}


// must substitute # -> ~
char* LibraProteinRatioParser::getPeptideString(Array<const char*>* peps, const char* link) {

  const int max_pep_length = 500;
  char encoded[max_pep_length];

  char* output = NULL;
  if(peps == NULL || peps->length() == 0)
    return output;

  int k, totlen= 0;
  for(k = 0; k < peps->length(); k++) {
    totlen += (int)strlen((*peps)[k]);
    if(k > 0)
      totlen += (int)strlen(link);
  }
  output = new char[totlen+1];

  if(strstr((*peps)[0], "#") == NULL)
    strcpy(output, (*peps)[0]);
  else {
    for(int j = 0; ((*peps)[0])[j]; j++)
      if(((*peps)[0])[j] == '#')
	encoded[j] = '~';
      else encoded[j] = ((*peps)[0])[j];
    encoded[strlen((*peps)[0])] = 0;
    strcpy(output, encoded);
  }

  for(k = 1; k < peps->length(); k++) {
    strcat(output, link);
    if(strstr((*peps)[k], "#") == NULL)
      strcat(output, (*peps)[k]);
    else {
      for(int j = 0; j < ((*peps)[k])[j]; j++)
	if(((*peps)[k])[j] == '#')
	  encoded[j] = '~';
	else encoded[j] = ((*peps)[k])[j];
      encoded[strlen((*peps)[k])] = 0;
      strcat(output, encoded);
    }
  }
  output[totlen] = 0;

  return output;
}

void LibraProteinRatioParser::getRatio(const char* protein_name, Array<const char*>* peptides, double dProbability) {
  double mean = 0.0;
  double meansq = 0.0;
  int num = 0;
  /*
    cout << "inputfiles..." << endl;
    for(int k = 0; k < input_xmlfiles_->length(); k++)
    cout << (*input_xmlfiles_)[k] << endl;

    cout << "peptides: " << endl;
    for(int k = 0; k < peptides->length(); k++)
    cout << (*peptides)[k] << endl;
  */

  if(libraGroupPeptideParser_ == NULL) {
    // get the LIBRA result for this set of peptides
    libraGroupPeptideParser_ = new LibraGroupPeptideParser(input_xmlfiles_, peptides, dProbability, norm_channel_, minimumThreshholdIntensity);

    if(libraGroupPeptideParser_ != NULL) {
      result_tags_ = libraGroupPeptideParser_->getProtXMLTags();
      delete libraGroupPeptideParser_;
      libraGroupPeptideParser_ = NULL;

    } else {
      cout << "Error: null parser for inputfiles: ";
      int k;
      for(k = 0; k < input_xmlfiles_->length(); k++)
	cout << (*input_xmlfiles_)[k] << " ";
      cout << " and peptides: ";
      for(k = 0; k < peptides->length(); k++)
	cout << (*peptides)[k] << " ";
      cout << endl;
      exit(1);
    }
  }
  else {
    if (dbug) {
      cout << "peps:";
      for(int k = 0; k < peptides->length(); k++)
	cout << (*peptides)[k] << "+";
      cout << endl;
    }
    result_tags_ = libraGroupPeptideParser_->getProtXMLTags(protein_name,peptides);
    if (dbug)
      cout << "<--"<<endl;
  }

}
