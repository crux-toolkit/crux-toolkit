/*
Program       : XML Parser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: Parser.cpp 8017 2020-02-12 18:34:06Z mhoopmann $

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

#include "Parser.h"
#include "Util/RACI/RACI.h" // for reading gzipped files with efficient seeks
#include <time.h>
#include <fcntl.h>
#include <errno.h>

Parser::Parser() {
  name_ = NULL;
  time_ = NULL;
  xmlfile_ = NULL;
  success_ = false;
}

Parser::Parser(const char* name) {
  if(name != NULL) {
    name_ = new char[strlen(name)+1];
    strcpy(name_, name);
  }
  else
    name_ = NULL;
  time_ = NULL;
  xmlfile_ = NULL;
  success_ = false;
}
/*
Parser::Parser(char* xmlfile) {
  init(xmlfile);
}
*/

Parser::~Parser() {
  if(time_ != NULL)
    delete[] time_;

  if(name_ != NULL)
    delete[] name_;

  if (xmlfile_ != NULL)
    delete[] xmlfile_;
}

void Parser::init(const char* xmlfile) {
  time_ = getDateTime(USE_LOCAL_TIME);
  success_ = false;
  if(time_ == NULL) {
    cout << "error: null analysis time" << endl;
    exit(1);
  }  
  filter_ = False;
  filter_memory_ = False;
  if(xmlfile != NULL) {
    xmlfile_ = new char[strlen(xmlfile)+1];
    strcpy(xmlfile_, xmlfile);
    parse(xmlfile_);
  }
  else
    parse(xmlfile);
}

void Parser::parse_and_print(const char* xmlfile) {
  //open file and pass along
  char *data = NULL;
  Tag* tag;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "Parser: error opening " << xmlfile << endl;
    exit(1);
  }

  char *nextline = new char[line_width_];
  while(fin.getline(nextline, line_width_)) {
    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      setFilter(tag);
      if(filter_) {
	tag->print();
      }
      delete tag;
      data = strstr(data+1, "<");
    }

  }
  fin.close();
  delete[] nextline;
}

Boolean Parser::setTagValue(Tag* tag, char* tagname, char* attr_name, int* index) {
  char text[10];
  if(! tag->isStart() || strcmp(tag->getName(), tagname))
    return False; // nothing to do
  //cout << "still here with " << *index << endl;
  sprintf(text, "%d", (*index)++);
  tag->setAttributeValue(attr_name, text);
  return True;
}

void Parser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  filter_ = True;
}


char* Parser::getUniqueNameSpace(Tag* tag, char* orig) {
  char* output = new char[strlen(orig)+3];
  strcpy(output, orig);
  char prefix[] = "xlmn:";
  char name[200];
  int index = 0;
  sprintf(name, "%s%s", prefix, output);
  while(index < 100 && tag->getAttributeValue(name) != NULL) {
    if(index == 0 || index == 10) {
      output[strlen(output)+1] = 0;
      orig[strlen(output)] = '0'; // enlarge
    }
    else 
      output[strlen(output)-1]++;
    sprintf(name, "%s%s", prefix, output);
    index++;
  } // while
  if(index == 100) {
    if(output != NULL)
      delete[] output;
    return NULL;
  }
  return output;
}


char* Parser::getOrigFile(char* file) {
  char suffix[] = ".orig";
  char* result = strstr(file, suffix);
  char* output = NULL;
  if(result != NULL && strlen(result) == strlen(suffix)) {
    output = new char[strlen(file)+1];
    strcpy(output, file);
  }
  else {
    output = new char[strlen(file)+strlen(suffix)+1];
    strcpy(output, file);
    strcat(output, suffix);
  }
  return output;
}

const char* Parser::getName() {
  return name_;
}

#ifdef _MSC_VER
#include <fcntl.h>
#endif

Boolean Parser::tag_is_at_tail(const char* xmlfile, const char* last_tag) {
  Boolean result = False;
  const int SIZE_BUF2 = 10000;
  char *szBuf = new char[SIZE_BUF2+1];
  int handle = open(xmlfile,O_RDONLY);

  if ( -1 == handle) {
    printf(" Error - cannot open input file %s\n\n", xmlfile);
    exit(1);
  }
#ifdef _MSC_VER
  _lseeki64(handle,-SIZE_BUF2,SEEK_END);
#else
  lseek(handle,-SIZE_BUF2,SEEK_END);
#endif
  int nread = read(handle,szBuf, SIZE_BUF2);
  close(handle);
  if (nread>0) {
    // look for tag in the last few lines
    szBuf[nread] = 0;
    int nlines = 0;
    for (int i=nread;i-- && (nlines < 4);) {
      if (szBuf[i]=='\n') {
	nlines++;
      }
      if (strstr(szBuf+i, last_tag) != NULL) {
	result = true;
	break;
      }
    }
  }
  delete[] szBuf;
  return result;
}

Boolean Parser::overwrite(const char* xmlfile, const char* outfile, const char* last_tag) {
  Boolean result = False;
  for (int retry = 3;retry-- && !result;) {
    if (tag_is_at_tail(outfile,last_tag)) {
      result= True; // file was written properly
      if (strcmp(xmlfile, outfile) != 0) { // not the same filename
	if (safe_rename(outfile, xmlfile)) {
	  printf("failed to update %s from temp file %s, error %d \"%s\"",
		 xmlfile,outfile,errno,strerror(errno));
	  result = False; // rename failed
	}
	if (result && isDotGZ(xmlfile)) { // perform gzip if needed
	  std::string x(xmlfile);
	  do_gzip(x);
	}
	chmod(xmlfile, 00664);  // why is this necessary?
      }
    } // end if tag is at tail
    else if (retry) {
      // slow filesystem?
      sleep(2);
    }
    else {
      printf("Error: file update failed (did not find closing tag \"%s\"), check file %s for completeness\n", last_tag, outfile);
    }
  }
  success_ = result; // set flag for access by Parser::success()
  return result;
}

char* Parser::getDateTime(Boolean local) {
  // 2002-07-23T23:04:44
  time_t now;
  time(&now);
  struct tm* tmstruct = NULL;
  if(local)
    tmstruct = localtime(&now);
  else 
    tmstruct = gmtime(&now);
  char* output = new char[100];
  sprintf(output, "%04d-%02d-%02dT%02d:%02d:%02d", tmstruct->tm_year+1900, tmstruct->tm_mon+1, tmstruct->tm_mday, tmstruct->tm_hour, tmstruct->tm_min, tmstruct->tm_sec);
  return output;
}
