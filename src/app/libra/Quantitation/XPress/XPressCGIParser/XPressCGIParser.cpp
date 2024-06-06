/*

Program       : XPressCGIParser                                                 
Author        : Andrew Keller <akeller@systemsbiology.org> 
                Jimmy Eng (jeng@systemsbiology.org>                                                      
Date          : 11.27.02 

Overwrites specified modified XPRESS protein ratio onto
ProteinProphet XML

Copyright (C) 2003 Andrew Keller, Jimmy Eng

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


#include "XPressCGIParser.h"
#include "Util/RACI/RACI.h"

XPressCGIParser::XPressCGIParser(const char* xmlfile, const char* protein, double ratio, double error, double h2l_ratio, double h2l_error, int numpeps) { 
  protein_ = new char[strlen(protein)+1];
  strcpy(protein_, protein);
  sprintf(ratio_, "%0.2f", ratio);
  sprintf(error_, "%0.2f", error);
  sprintf(h2l_ratio_, "%0.2f", h2l_ratio);
  sprintf(h2l_error_, "%0.2f", h2l_error);
  sprintf(numpeps_, "%d", numpeps);
  init(xmlfile);
}

void XPressCGIParser::parse(const char* xmlfile) {
  //open file and pass along
  const int line_width = 10000;
  char nextline[line_width];
  char* data = NULL;
  Tag* tag;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cout << "XPressCGIParser: error opening " << xmlfile << endl;
    exit(1);
  }

  Boolean done = False;

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);
  //cerr << "writing data to " << outfile << endl;
  ofstream fout(outfile.c_str());
  if(! fout) {
    cout << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  while(fin.getline(nextline, line_width)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      setFilter(tag);

      if(! done && filter_) {

	if(tag->isStart() && ! strcmp(tag->getName(), "XPressRatio")) {

	  tag->setAttributeValue("ratio_mean", ratio_);
	  tag->setAttributeValue("ratio_standard_dev", error_);
	  tag->setAttributeValue("heavy2light_ratio_mean", h2l_ratio_);
	  tag->setAttributeValue("heavy2light_ratio_standard_dev", h2l_error_);
	  tag->setAttributeValue("ratio_number_peptides", numpeps_);
	  
	  done = True;
	      
	}

      }

      tag->write(fout);

      delete tag;
      data = strstr(data+1, "<");
    }

  }
  fin.close();
  fout.close();

  if(overwrite(xmlfile, outfile.c_str(), "</protein_summary>")) {
    cout << "changes written to file " << xmlfile << endl;
    cout << "refresh ProteinProphet XML Viewer to display" << endl << endl;
  }
  else {
    cout << "error: no changes written to file " << xmlfile << endl;
  }
}

void XPressCGIParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(tag->isStart() && ! strcmp(tag->getName(), "protein") && 
     ! strcmp(tag->getAttributeValue("protein_name"), protein_)) 
      filter_ = True;

}


