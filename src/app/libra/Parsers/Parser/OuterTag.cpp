/*
Program       : OuterTag
Author        : ?
Date          : ?
SVN info      : $Id: OuterTag.cpp 7714 2017-12-20 00:16:03Z real_procopio $


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

#include "OuterTag.h"
#include <stdlib.h>
#include <assert.h>

OuterTag::OuterTag(const Tag* tag) : Tag(*tag) {

  // now determine the namespaces....
  for(int k = 0; k < getAttributeCount(); k++)
    if(! strncmp(getAttribute(k), "xmlns", strlen("xmlns"))) {
      //cout << attributes_[k] << ": " << values_[k] << endl;
      namespace_prefixes_.insertAtEnd(strCopy(getAttribute(k)));
      namespaces_.insertAtEnd(strCopy(values_[k]));
    }
    else if(strstr(getAttribute(k), "schemaLocation") != NULL) {
      //cout << attributes_[k] << ": " << values_[k] << endl;

      char text[500];
      text[0] = 0; // reset
      char text2[500];
      int start = 0;
      int stop = 0;
      while(values_[k][stop]) {
	if(values_[k][stop] == ' ') {
	  if(!text[0]) {
	    strncpy(text, values_[k]+start, stop-start);
	    text[stop-start] = 0;
	    start = stop+1;
	  }
	  else {
	    strncpy(text2, values_[k]+start, stop-start);
	    text2[stop-start] = 0;
	    start = stop+1;
	    referenced_namespaces_.insertAtEnd(strCopy(text));
	    referenced_schemas_.insertAtEnd(strCopy(text2));
	    text[0] = 0; // reset
	  }
	}
	stop++;
      } // while
      // last
      strcpy(text2, values_[k]+start);
      referenced_namespaces_.insertAtEnd(strCopy(text));
      referenced_schemas_.insertAtEnd(strCopy(text2));

      //for(int k = 0; k < referenced_schemas_->size(); k++)
      //	cout << (*referenced_schemas_)[k] << " <> " << (*referenced_namespaces_)[k] << endl;

      /*
      while(nextpair != NULL && nextpair[0]) {
	sscanf(nextpair, "%s %s", text, text2);
	char* next = new char[strlen(text)+1];
	strcpy(next, text);
	referenced_namespaces_->insertAtEnd(next);
	next = new char[strlen(text2)+1];
	strcpy(next, text2);
	referenced_schemas_->insertAtEnd(next);
	nextpair = strstr(nextpair, text);
	if(nextpair != NULL)
	  nextpair += strlen(text);
      } // while
      */
    } // if schemaLocation
#ifdef DEBUG
  assert(namespaces_.size() == namespace_prefixes_.size());
  assert(referenced_namespaces_.size() == referenced_schemas_.size());
#endif
}

OuterTag::~OuterTag() {
}


void OuterTag::enterRefs(const OuterTag* outer) {
  for(int k = 0; k < outer->getNumNameSpaces(); k++) {
    for(int j = 0; j < outer->getNumRefs(); j++) {
      if(! strcmp(outer->getNameSpace(k), outer->getRefNameSpace(j))) {
	enter(outer->getNamespace_prefix(k), outer->getNameSpace(k), outer->getRefSchema(j));
	j = outer->getNumRefs();
      }
    } // next ref
  } // next ns
}

int OuterTag::getNumRefs() const {
  return referenced_schemas_.size();
}
int OuterTag::getNumNameSpaces() const {
  return namespace_prefixes_.size();
}
const char* OuterTag::getNamespace_prefix(int k) const {
  return namespace_prefixes_[k];
}
const char* OuterTag::getRefNameSpace(int k) const {
  return referenced_namespaces_[k];
}
const char* OuterTag::getNameSpace(int k) const {
  return namespaces_[k];
}
const char* OuterTag::getRefSchema(int k) const {
  return referenced_schemas_[k];
}
bool verbose = false;
void OuterTag::enter(const char* namespace_prefix, const char* namesp, const char* ref_schema) {
  if (verbose)
    cout << "entering ns: " << namespace_prefix << " namesp: " << namesp << " ref " << ref_schema << endl;
  // check if already exists under different namespace
  int k;
  for(k = 0; k < namespaces_.size(); k++) {
    if(! strcmp(namespaces_[k], namespace_prefix)) {
      if(strcmp(namespace_prefixes_[k], namesp)) {
	cout << "error: " << namespace_prefixes_[k] << " associated with multiple namespaces: " << namespaces_[k] << " and " << namesp << endl;
	exit(1);
      }
    }
  }

  // make attribute for this new guy
  setAttributeValue(namespace_prefix, namesp);
  for(k = 0; k < getAttributeCount(); k++)
    if(strstr(getAttribute(k), "schemaLocation") != NULL) {
      char* text = new char[strlen(namesp)+strlen(ref_schema)+2];
      sprintf(text, "%s %s", namesp, ref_schema);
      appendAttributeValue(getAttribute(k), text, True);
      delete[] text;
    }

  namespace_prefixes_.insertAtEnd(strCopy(namespace_prefix));
  namespaces_.insertAtEnd(strCopy(namesp));
  referenced_namespaces_.insertAtEnd(strCopy(namesp));
  referenced_schemas_.insertAtEnd(strCopy(ref_schema));
#ifdef _DEBUG
  assert(namespaces_.size() == namespace_prefixes_.size());
  assert(referenced_namespaces_.size() == referenced_schemas_.size());
#endif
}
