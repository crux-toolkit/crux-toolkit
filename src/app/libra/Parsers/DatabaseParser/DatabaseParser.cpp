/*

Program       : DatabaseParser                                                   
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 
SVN Info      : $Id: DatabaseParser.cpp 8800 2023-01-06 09:26:53Z real_procopio $

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

#include "DatabaseParser.h"
#include "Util/RACI/RACI.h" // for reading gzipped files with efficient seeks

DatabaseParser::DatabaseParser(const char* xmlfile) : Parser(NULL) {
  databases_ = new Array<char*>;
  init(xmlfile);
}

void DatabaseParser::parse(const char* xmlfile) {
  char current_database[500];
  Tag* tag = NULL;
  char* data = NULL;

  current_database[0] = 0;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "error opening " << xmlfile << endl;
    exit(1);
  }
  char *nextline = new char[line_width_];
  while(fin.getline(nextline, line_width_)) {

    if(strstr(nextline, "msms_run_summary") != NULL ||
       strstr(nextline, "refresh_timestamp") != NULL ||
       strstr(nextline, "search_database") != NULL) {
    
      data = strstr(nextline, "<");
      while(data != NULL) {
	tag = new Tag(data);

	if(tag != NULL) {
	  if(tag->isStart() && ! strcmp(tag->getName(), "search_database")) {
	    memcpy(current_database, 
		   tag->getAttributeValue("local_path"), 
		   strlen(tag->getAttributeValue("local_path"))+1);
	  }
	  else if(tag->isStart() && ! strcmp(tag->getName(), "database_refresh_timestamp")) {
	    strcpy(current_database, tag->getAttributeValue("database"));
	  }
	  if(current_database[0] && tag->isEnd() && ! strcmp(tag->getName(), "msms_run_summary")) { // process
	    enter(current_database);
	    current_database[0] = 0;
	  }
	  delete tag;
	} // if not null

	data = strstr(data+1, "<");
      } // next tag
    
    } // if have reson to parse tags
  } // next line
  fin.close();

  delete [] nextline;
}


Boolean DatabaseParser::enter(char* db) {
  for(int k = 0; k < databases_->length(); k++) {
    if(! strcmp(db, (*databases_)[k])) {
      return False; // already seen
    }
    // uh oh - possible conflict, examine closer
    std::string dbk((*databases_)[k]);
    std::string dbnew(db);
    unCygwinify(dbk); // remove any old cygwin stuff
    unCygwinify(dbnew); // remove any old cygwin stuff
    std::string dbk0 = resolve_root(dbk.c_str());
    std::string dbnew0 = resolve_root(dbnew.c_str());
    int diff = strcmp(dbk0.c_str(),dbnew0.c_str());
    if (!diff) {
      return False; // already seen
    }
  }
  char* next = new char[strlen(db)+1];
  strcpy(next, db);
  databases_->insertAtEnd(next);
  return True;
} 

std::string DatabaseParser::getDatabases() {
  std::string output;
  for(int k = 0; k < getNumDatabases(); k++) {
    if (k) {
      output += ",";
    }
    output += (*databases_)[k];
  }
  return output;
}

int DatabaseParser::getNumDatabases() { 
  if(databases_ == NULL)
    return 0;
  return databases_->length();
}
