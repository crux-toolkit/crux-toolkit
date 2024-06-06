/*
Program       : ASAPCGIParser
Authors       : Andrew Keller <akeller@systemsbiology.org>
                Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPCGIParser.cpp 8445 2021-04-20 01:01:01Z real_procopio $

Overwrites modified ASAPRatio protein information to ProteinProphet XML

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


#include "ASAPCGIParser.h"
#include "Util/RACI/RACI.h"

ASAPCGIParser::ASAPCGIParser(const char* xmlfile, const char* protein, proDataStrct* pro_ratio) : Parser(NULL) { 
  pro_ratio_ = pro_ratio;
  heavy2light_ = False;
  protein_ = new char[strlen(protein)+1];
  strcpy(protein_, protein);

  init(xmlfile);
}

void ASAPCGIParser::parse(const char* xmlfile) {
  //open file and pass along
  char *nextline = new char[line_width_];
  char* data = NULL;
  Tag* tag;
  Boolean heavy2light = False;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cout << "ASAPCGIParser: error opening " << xmlfile << endl;
    exit(1);
  }

  norm_ = NULL;

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);

  ofstream fout(outfile.c_str());
  if(! fout) {
    cout << "cannot write output to file " << outfile << endl;
    exit(1);
  }
  Boolean boycott = False;

  while(fin.getline(nextline, line_width_)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      setFilter(tag);
      //tag->write(cout);
      if(! strcmp(tag->getName(), "ASAP_pvalue_analysis_summary")) {
	tag->write(fout); 
	// get the parameters
	pValueStrct params;
	params.mean = atof(tag->getAttributeValue("background_ratio_mean"));
	params.stddev = atof(tag->getAttributeValue("background_ratio_stdev"));
	params.merr = atof(tag->getAttributeValue("background_fitting_error"));
	norm_ = new Normalization(params);
      }

      else if(filter_) { 
	if(! strcmp(tag->getName(), "ASAPRatio")) {
	  if(tag->isStart()) {
	    // display the new one and remove all the old ones
	    writeProteinRatio(fout, pro_ratio_, tag->getName());
	    boycott = True;
	  }
	  if(tag->isEnd()) {
	    boycott = False;
	  }
	} // if start of ratio
	else if(! strcmp(tag->getName(), "ASAPRatio_pvalue")) {
	  if(tag->isStart()) {
	    // display the new one and remove all the old ones
	    writeProteinRatio(fout, pro_ratio_, tag->getName());
	    boycott = True;
	  }
	  if(tag->isEnd()) {
	    boycott = False;
	  }
	} // if start of ratio
	else if(! boycott)
	  tag->write(fout); 
      } // if filter and ready....
      else
	tag->write(fout); 

      delete tag;
      data = strstr(data+1, "<");
    }
  }
  // overwrite here
  fout.close();
  fin.close();
  //return;

  if(overwrite(xmlfile, outfile.c_str(), "</protein_summary>")) {
    //cout << "changes written to file " << xmlfile << "<br>" << endl;
    //cout << "<font color=\"red\">refresh ProteinProphet XML Viewer to display</font>" << endl << endl;
    cout << "<script language='JavaScript'>"
	 << "messages.push(\"Changes written to file " << xmlfile << "\");"
	 << "messages.push(\"Refresh ProtXMLViewer to display updated analysis.\");"
	 << "</script>" << endl;
  }
  else {
    //cout << "error: no changes written to file " << xmlfile << endl;
    cout << "<script language='JavaScript'>"
	 << "messages.push(\"Error: No changes written to file " << xmlfile << "\");"
	 << "</script>" << endl;
  }
  delete[] nextline;  
}


// from beginning to end of each protein
void ASAPCGIParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(tag->isStart() && ! strcmp(tag->getName(), "protein") && ! strcmp(tag->getAttributeValue("protein_name"), protein_)) 
    filter_ = True;

  else if(filter_ && ! strcmp(tag->getName(), "protein") && ! tag->isStart())
    filter_memory_ = True;
}

void ASAPCGIParser::writeProteinRatio(ostream& os, proDataStrct* pro_ratio, const char* name) {
  Tag* protag = NULL;
  char next[50];
  double adjratio[2];
  double adj_inv_ratio[2];
  double pvalue;
  
  double ratiomean = (double)pro_ratio->ratio[0];
  //  cout << "ratiomean: " << ratiomean << endl;
  if(ratiomean == -2.0) 
    sprintf(next, "%0.2f", -1.0);
  else if(ratiomean == -1.0) 
    strcpy(next, "999.");
  else {
    sprintf(next, "%0.2f", ratiomean);
  }
  if(! strcmp(name, "ASAPRatio_pvalue") && norm_ != NULL) {
    protag = new Tag("ASAPRatio_pvalue", True, True);
    adjratio[0] = (double)ratiomean; // decode
    adjratio[1] = pro_ratio->ratio[1];

    //TODO: DDS should really model both h2l and l2h ratios with own adjustment models
    //But for now we'll first invert, then adjust, then invert again
    adj_inv_ratio[0] = 1/pro_ratio->inv_ratio[0];
    adj_inv_ratio[1] = pro_ratio->inv_ratio[1]/(pro_ratio->inv_ratio[0]*pro_ratio->inv_ratio[0]);
    pvalue = (double)norm_->normalize(adj_inv_ratio);

    double tmp[2];
    tmp[0] = adj_inv_ratio[0];
    tmp[1] = adj_inv_ratio[1];

    adj_inv_ratio[0] = 1/ tmp[0];
    adj_inv_ratio[1] = tmp[1]/(tmp[0]*tmp[0]);

    pvalue = (double)norm_->normalize(adjratio);

    sprintf(next, "%0.2f", adjratio[0]);
    protag->setAttributeValue("adj_ratio_mean", next);
    sprintf(next, "%0.2f", adjratio[1]);
    protag->setAttributeValue("adj_ratio_standard_dev", next);

    sprintf(next, "%0.2f", adj_inv_ratio[0]);
    protag->setAttributeValue("heavy2light_adj_ratio_mean", next);
    if (isnan(adj_inv_ratio[1])) {
      protag->setAttributeValue("heavy2light_adj_ratio_standard_dev", "NaN");
    }
    else {
      sprintf(next, "%0.2f", adj_inv_ratio[1]);
      protag->setAttributeValue("heavy2light_adj_ratio_standard_dev", next);
    }

    if(pvalue >= 0.0) {
      sprintf(next, "%0.2e", pvalue);
      protag->setAttributeValue("pvalue", next);
      sprintf(next, "%0.8f", pvalue);
      protag->setAttributeValue("decimal_pvalue", next);
    }
    protag->write(os);
    // now the protein stuff
    delete protag;
    return;
  } // asap pvalue

  if(strcmp(name, "ASAPRatio")) 
    return; // go away

  protag = new Tag("ASAPRatio", True, False);
  protag->setAttributeValue("ratio_mean", next);
  sprintf(next, "%0.2f", (double)pro_ratio->ratio[1]);
  protag->setAttributeValue("ratio_standard_dev", next);
  sprintf(next, "%d", pro_ratio->dataNum);
  protag->setAttributeValue("ratio_number_peptides", next);

  sprintf(next, "%0.2f", (double)pro_ratio->inv_ratio[0]);
  protag->setAttributeValue("heavy2light_ratio_mean", next);
  if (isnan(pro_ratio->inv_ratio[1])) {
    protag->setAttributeValue("heavy2light_ratio_standard_dev", "NaN");
  }
  else {
    sprintf(next, "%0.2f", (double)pro_ratio->inv_ratio[1]);
    protag->setAttributeValue("heavy2light_ratio_standard_dev", next);
  }

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
    sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].inv_ratio[0]);
    seqtag->setAttributeValue("heavy2light_ratio_mean", next);
    if (isnan((double)pro_ratio->sequences[k].inv_ratio[1])) {
      seqtag->setAttributeValue("heavy2light_ratio_standard_dev", "NaN");
    }
    else {
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].inv_ratio[1]);
      seqtag->setAttributeValue("heavy2light_ratio_standard_dev", next);
    }
    sprintf(next, "%0.4f", pro_ratio->sequences[k].weight);
    seqtag->setAttributeValue("weight", next);
    seqtag->setAttributeValue("light_sequence", pro_ratio->sequences[k].lightSeq);
    seqtag->write(os);
    delete seqtag;

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
      sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].inv_ratio[0]);
      peaktag->setAttributeValue("heavy2light_ratio_mean", next);
      if (isnan((double)pro_ratio->sequences[k].peaks[j].inv_ratio[1])) {
	peaktag->setAttributeValue("heavy2light_ratio_standard_dev", "NaN");
      }
      else {
	sprintf(next, "%0.4f", (double)pro_ratio->sequences[k].peaks[j].inv_ratio[1]);
	peaktag->setAttributeValue("heavy2light_ratio_standard_dev", next);
      }
      sprintf(next, "%0.1f", (double)pro_ratio->sequences[k].peaks[j].weight);
      peaktag->setAttributeValue("weight", next);

      //DDS:
      //sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].msms_run_idx);
      sprintf(next, "%d", pro_ratio->sequences[k].peaks[j].bofIndx);
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
