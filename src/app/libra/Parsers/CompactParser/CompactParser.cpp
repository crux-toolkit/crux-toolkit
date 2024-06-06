/*

Program       : CompactParser                                                   
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

Primary data object holding all mixture distributions for each precursor ion charge

Copyright (C) 2003 Andrew Keller

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

#include "CompactParser.h"
#include "Util/RACI/RACI.h" // for reading gzipped files with efficient seeks

CompactParser::CompactParser(const char* xmlfile) : Parser(NULL) {

  init(xmlfile);
}

void CompactParser::parse(const char* xmlfile) {

  Tag* tag = NULL;

  //  int line_width = 10000;
  char* data = NULL;

  Tag* lasttag = NULL;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "CompactParser: error opening " << xmlfile << endl;
    exit(1);
  }

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);

  //cerr << "writing data to " << outfile << endl;
  ofstream fout(outfile.c_str());
  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }




  char *nextline = new char[line_width_];
  while(fin.getline(nextline, line_width_)) {
    //cout << "next: " << nextline << endl;

    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);

      //tag->write(cout);

      // check if same as last tag
      if(lasttag != NULL && tag != NULL) {

	if(tag->isEnd() && ! strcmp(lasttag->getName(), tag->getName())) { // combine first
	  lasttag->setEnd(); // make an end
	  lasttag->write(fout);
	  delete lasttag;
	  lasttag = NULL;
	  // now set this tag to null (do not write it)
	  //delete tag;
	    
	}
	else {

	  lasttag->write(fout);
	  delete lasttag;
	  lasttag = tag;
	}
      }
      else {
	lasttag = tag;
      }
      data = strstr(data+1, "<");
    } // next tag

  } // next line

  fin.close();
  if(lasttag != NULL) {
    lasttag->write(fout);
    delete lasttag;
  }

  fout.close();
  if(!(tag_is_at_tail(outfile.c_str(), "</msms_pipeline_analysis>") && overwrite(xmlfile, outfile.c_str(), "</protein_summary>"))) {
    cerr << "error: no data written to file " << xmlfile << endl;
  }

  delete [] nextline;
}


void CompactParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "search_result")){ 
    if(tag->isStart()) {
      //tag->print();
      filter_ = True;
    }else{
      filter_memory_ = True;
    }
  }

}

