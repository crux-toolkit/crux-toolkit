/*
Program       : Tag
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: Tag.cpp 7869 2018-11-28 19:47:54Z dshteyn $

Copyright (C) 2003 Andrew Keller

Performance work Copyright (C) 2007 Insilicos and Labkey

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

#include "Tag.h"
#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <ctype.h>


#ifndef NO_USETABLES
// reduce memory alloc thrash by exploiting the fact that
// most XML tags are heavily repeated
// don't make copies of every tag and attribute name
// WARNING not threadsafe (but not much is in TPP-land)


#include <vector>

static tagname_set tagnames;
static tagname_set attrnames;

int compare_charptr_int_pairs(const void *a,const void *b) {
   return strcmp(((charptr_int_pair *)a)->charptr_,((charptr_int_pair *)b)->charptr_);
}

#define CLEANUP_ON_EXIT
#ifdef CLEANUP_ON_EXIT
// clean up on exit
bool firstTag = true;
static void clearTagNames() {
  std::vector<const char *>ptrs;
  for (tagname_set::iterator i = tagnames.begin();i!=tagnames.end();i++) {
    ptrs.push_back(*i);
  }
  for (tagname_set::iterator j = attrnames.begin();j!=attrnames.end();j++) {
    ptrs.push_back(*j);
  }
  tagnames.clear();
  attrnames.clear();
  for (size_t n=ptrs.size();n--;) {
    delete[] (char *)(ptrs[n]); // cast away const
  }
}
#endif

static const char *makeTagName(const char *name, int length=0) {
  if (!name || !*name) {
    return NULL;
  }
  tagname_set::iterator iter;
  if (length) { // cut it out of the buffer
    char hold = name[length];
    ((char *)name)[length] = 0; // casting away const, but we'll restore
    iter = tagnames.find( name );
    ((char *)name)[length] = hold; // restore
  } else {
    iter = tagnames.find( name );
  }
  if (iter == tagnames.end()) {
#ifdef CLEANUP_ON_EXIT
    if (firstTag) {
      firstTag = false;
      atexit(clearTagNames);
    }
#endif
    if (!length) {
      length = (int)strlen(name);
    }
    char *result = new char[length+1];
    memmove(result,name,length);
    result[length] = 0;
    tagnames.insert(result);
    return result;
  } else {
    return *iter;
  }
}

static const char *makeAttributeName(const char *name, int length=0) {
  if (!name || !*name) {
    return NULL;
  }
  tagname_set::iterator iter;
  if (length) { // cut it out of the buffer
    char hold = name[length];
    ((char *)name)[length] = 0; // casting away const, but we'll restore
    iter = attrnames.find( name );
    ((char *)name)[length] = hold; // restore
  } else {
    iter = attrnames.find( name );
  }
  if (iter == attrnames.end()) {
    if (!length) {
      length = (int)strlen(name);
    }
    char *result = new char[length+1];
    memmove(result,name,length);
    result[length] = 0;
    attrnames.insert(result);
    return result;
  } else {
    return *iter;
  }
}

#define copyAttributeName(s) s
#define copyTagName(s) s



#else
// literal copy of tags and attrs, no reuse
static char *makeTagName(const char *buf) {
  return strCopy(buf);
}
static char *makeTagName(const char *buf, int len) {
  char *result = new char[len+1];
  memmove(result,buf,len);
  result[len] = 0;
  return result;
}
#define copyTagName strCopy
#define makeAttributeName makeTagName
#define copyAttributeName strCopy
#endif



Tag::Tag(const char* data,const TagInclusionTest* inclusionTest) {
  start_ = False;
  end_ = False;
  namespace_ = NULL;
  line_number_ = 0;

  // must parse through
  if(data[0] == '<' && data[1]) {
    size_t start;
    size_t stop;

    if(data[1] == '?') {
      start_ = False;
      end_ = False; // special case
      start=2;
      stop=2;
    }
    else if(data[1] == '/') {
      start=2;
      stop=2;
      end_ = True;
    }
    else {
      start = 1;
      stop = 1;
      start_ = True;
    }

    int k = (int)start;
    while(data[k] && 
	  (data[k] != ' ' && data[k] != '>' && 
	   (!data[k+1] || data[k] != '/' || data[k+1] != '>'))) 
      k++;

    stop = k;
    bool isBang; // watch out for directives which may contain colons in path
    if((isBang=(data[1] == '!'))) { // special case
      k = (int)strlen(data)-1;
      while(k >= 0 && data[k] != '>')
	k--;
      stop = k;
    } else {
      k = (int)start;
      while(data[k] && 
	    (data[k] != ' ' && data[k] != '>' && 
	     (data[k] != '/' || data[k+1] != '>'))) 
	k++;
      stop = k;
    }

    // check for namespace here
    size_t i,divisor = 0;
    for(i = start; i < stop; i++) {
      if(data[i] == ':') {
	divisor = i;
      }
    }

    if((divisor > 0)&&!isBang) { // namespace
      if (inclusionTest && !inclusionTest->includesName(data+divisor+1,(int)(stop - divisor -1))) {
	name_ = NULL;
	return;
      }
      namespace_ = makeTagName(data+start,(int)(divisor - start));
      name_ = makeTagName(data+divisor+1,(int)(stop - divisor -1));
    }
    else {
      if (inclusionTest && !inclusionTest->includesName(data+start,(int)(stop - start))) {
	name_ = NULL;
	return;
      }
      name_ = makeTagName(data+start,(int)(stop - start));
    }

    // now read the attributes
    while(readAttribute(data, k, inclusionTest));

    if(! end_) {
      int final = 0;
      Boolean open_sn_quote = False;
      Boolean open_db_quote = False;
      while(data[final] && 
	    (data[final] != '>' || open_db_quote || open_sn_quote)) {
	if (open_sn_quote) {
	  if (data[final] == '\'') { 
	    open_sn_quote = False;
	  }
	}
	else if (open_db_quote) {
	  if (data[final] == '\"') {
	    open_db_quote = False;
	  }
	}
	else {
	  if (data[final] == '\"') {
	    open_db_quote = True;
	  }
	  else if (data[final] == '\'') {
	    open_sn_quote = True;
	  }
	}
	final++;
      }
      end_ = data[final] == '>' && final > 0 && data[final-1] == '/';

    }
  }
}

Tag::Tag(const char* name, Boolean start, Boolean end) {
  start_ = start;
  end_ = end;
  name_ = makeTagName(name);
  namespace_ = NULL;
  line_number_ = 0;
}

Tag::Tag(const char* name, Boolean start, Boolean end, const char* nspace) {
  start_ = start;
  end_ = end;
  name_ = makeTagName(name);
  namespace_ = makeTagName(nspace);
  line_number_ = 0;
}

Tag::Tag() {
  start_ = end_ = False;
  name_ = NULL;
  namespace_ = NULL;
  line_number_ = 0;
}

Tag::Tag(const Tag &rhs) {
  start_ = rhs.start_;
  end_ =rhs.end_; 
  name_ = copyTagName(rhs.name_);
  namespace_ = copyTagName(rhs.namespace_);
  attributes_.clear();
  for (int i=0;i<rhs.attributes_.size();i++) {
    attributes_.add(copyAttributeName(rhs.attributes_.getPair(i)->charptr_),rhs.attributes_.getPair(i)->index_);
  }
  values_ = rhs.values_;
  line_number_ = rhs.line_number_;
}

Tag::~Tag() {
#ifdef NO_USETABLES // are we doing the re-use thing?
  delete[] name_;
  delete[] namespace_;
  charptr_int_map::const_iterator iter(attributes_);
  while (iter!=attributes_.end()) {
    delete[] iter->first;
  }
#endif
}

void Tag::setName(const char *newname) {
  name_ = makeTagName(newname);
}

void Tag::setNamespace(const char *newnamespace) {
  namespace_ = makeTagName(newnamespace);
}

Tag* Tag::copy() const {
  Tag* output = new Tag(*this);
  return output;
}

void Tag::copyAttributesTo(Tag* output) const {
  if(output == NULL)
    return;
  for(int k = 0; k < attributes_.size(); k++)
    output->setAttributeValue(getAttribute(k), values_[k]);
}

const char* Tag::getAttributeValue(int n) const {
  if (n < attributes_.size()) {
    return values_[n];
  }
  // still here
  return NULL;
}


Boolean Tag::appendAttributeValue(const char* attribute, const char* value, Boolean unique) {
  int k = getAttributeIndex(attribute);
  if(k>=0) {
    if(values_[k] != NULL) {
      // now check if value is already present
      if(unique) {
	char* result = strstr(values_[k], value);
	if(result != NULL) { // check that exact match on flanks
	  if(
	     (strlen(result) == strlen(values_[k]) || (result-1)[0] == ' ') &&
	     (strlen(result) == strlen(value) || result[strlen(value)] == ' '))  // then already present
	    return False;
	  // cout << "still here...." << endl;
	  char* newvalue = new char[strlen(values_[k]) + strlen(value) + 2];
	  strcpy(newvalue, values_[k]);
	  strcat(newvalue, " ");
	  strcat(newvalue, value);
	  delete[] values_[k];
	  values_.replace(k, newvalue);
	  return True;
	} // if potential match
      } // if unique
      // still here?
      char* newvalue = new char[strlen(values_[k]) + strlen(value) + 2];
      strcpy(newvalue, values_[k]);
      strcat(newvalue, " ");
      strcat(newvalue, value);
      delete[] values_[k];
      values_.replace(k, newvalue);
      return True;
    } // not null
    else {
      values_.replace(k,strCopy(value));
      return True;
    }
  } // if atttribute match

  // must make new entry
  // still here
  insert_attribute(makeAttributeName(attribute),strCopy(value));
  return True;
}

void Tag::insert_attribute(const char *attr,char *value) {
  attributes_.add(attr); // set name,index pair - index set as attributes_.size() before add
  values_.insertAtEnd(value);
}

void Tag::setAttributeValue(const char* attribute, const char* value) {
  char* newvalue = strCopy(value);
  for (char *cp=newvalue;*cp;cp++) { // enforce legal XML
    if ((!isprint(*cp))&&(!isspace(*cp))) {
      *cp = ' ';
    }
  }
  // are we adding a new one, or replacing value for an old one?
  int k = getAttributeIndex(attribute);
  if (k>=0) {
    // replacing an old one
    delete[] values_[k];
    values_.replace(k, newvalue); 
    return;
  }
  // adding a new one
  const char* newattr = makeAttributeName(attribute);
  insert_attribute(newattr,newvalue);
}

void Tag::addAttributeValue(const char* attribute, const char* value) {
  char* newvalue = strCopy(value);
  for (char *cp=newvalue;*cp;cp++) { // enforce legal XML
    if ((!isprint(*cp))&&(!isspace(*cp))) {
      *cp = ' ';
    }
  }
  const char* newattr = makeAttributeName(attribute);
  insert_attribute(newattr,newvalue);
}

void Tag::setEnd() {
  end_ = True;
}

Boolean Tag::readAttribute(const char* data, int& index, const TagInclusionTest *test) {
  while(data[index] == ' ')
    index++;

  int att_start = index;
  int att_stop = index;
  int val_start = index;
  int val_stop = index;
  char* next = NULL;

  // go ahead till hit end of tag or second quote
  while(data[index] && data[index] != '=' && data[index] != '>')
    index++;

  if(!data[index] || data[index] == '>')
    return False;

  if(data[index] == '=') {
    att_stop = index;
    index++;
  }
  if((data[index] == '"' || data[index] == '\'')) {
    char del = data[index++]; // are we looking for " or ' as closing quote?
    val_start = index;
    while(data[index]) {
      if( (data[index++] == del ) ) {
	val_stop = index - 1;
	break;
      }
    }
  }
  if(val_stop > att_start) { // have one
    if (!test || test->includesAttribute(data+att_start,att_stop-att_start)) {
      // include if no test, or if listed
      next = new char[val_stop - val_start + 1];
      memcpy(next,data+val_start,val_stop - val_start);
      next[val_stop - val_start] = 0;
      //cout << next << endl;
      insert_attribute(makeAttributeName(data+att_start,att_stop-att_start),next);
    }
    return True;
  }
  return False;
}

void Tag::print() const {
  if(name_ == NULL || !name_[0])
    return;

  if(attributes_.size() > 0) {
    cout << "<TABLE frame=\"border\" rules=\"all\" cellpadding=\"2\" bgcolor=\"white\" style=\"font-family: \'Courier New\', Courier, mono; font-size: 10pt;\">" << endl;
    cout << "<tr><td><b><font color=\"brown\">parameter</font></b></td><td><b><font color=\"brown\">value</font></b></td></tr>" << endl;
    for(int k = 0; k < (int)attributes_.size(); k++)
      cout << "<tr><td>" << getAttribute(k) << "</td><td>" << values_[k] << "</td></tr>" << endl;

    cout << "</TABLE>" << endl;
  }
}


void Tag::write(ostream& os) const {
  if(name_ && name_[0]) {

    os << "<";
    if(! start_ && end_) {
      os << "/";
    } else if(! start_ && ! end_) {
      os << "?";
    }
  
    if(namespace_ != NULL) {
      os << namespace_ << ":";
    }
    os << name_;
    for(int k = 0; k < (int)attributes_.size(); k++) {
      os << " " << getAttribute(k) << "=\"" << values_[k] << "\"";
    }
    if(start_ && end_) {
      os << "/";
    } else if(! start_ && ! end_) {
      os << "?";
    }
    os << ">" << endl;
  }
  os.flush();
}

Boolean Tag::isIdentical(const Tag& rhs) const {
  if (!name_ && !rhs.name_) {
    return true; // both empty, so same
  }
  if (!name_ || !rhs.name_) {
    return false; // one is empty
  }
  if (strcmp(name_,rhs.name_)) {
    return false; // different name
  }
  if ((start_ != rhs.start_) || (end_ != rhs.end_)) {
    return false;
  }
  if (namespace_==NULL) {
    if (rhs.namespace_) {
      return false;
    }
  } else if (rhs.namespace_==NULL) {
    return false;
  } else if (strcmp(namespace_,rhs.namespace_)) {
    return false;
  }
  if (attributes_.size() != rhs.attributes_.size()) {
    return false;
  }
  for(int k = 0; k < (int)attributes_.size(); k++) {
    if (strcmp(getAttribute(k),rhs.getAttribute(k)) ||
	strcmp(values_[k],rhs.values_[k])) {
      return false;
    }
  }
  return true;
}

// now instantiate it for FILE and ogzstream
//template void Tag::writeTraditional(FILE * fp) const;
//template void Tag::writeTraditional(ogzstream * fp) const;

// for use with FILE or ogztream
//template <class T> void Tag::writeTraditional(T * fp) const {
void Tag::writeTraditional(ogzstream * fp) const {
  if(name_ == NULL || !name_[0])
    return;

  fprintf(fp, "<");
  if(! start_ && end_)
    fprintf(fp, "/");
  else if(! start_ && ! end_)
    fprintf(fp, "?");

  if(namespace_ != NULL)
    fprintf(fp, "%s:", namespace_);
 
  fprintf(fp, "%s", name_);

  for(int k = 0; k < (int)attributes_.size(); k++)
    fprintf(fp, " %s=\"%s\"",  getAttribute(k), values_[k]);

  if(start_ && end_)
    fprintf(fp, "/");
  else if(! start_ && ! end_)
    fprintf(fp, "?");

  fprintf(fp, ">\n");
}

void Tag::writeTraditional(FILE * fp) const {
  if(name_ == NULL || !name_[0])
    return;

  fprintf(fp, "<");
  if(! start_ && end_)
    fprintf(fp, "/");
  else if(! start_ && ! end_)
    fprintf(fp, "?");

  if(namespace_ != NULL)
    fprintf(fp, "%s:", namespace_);
 
  fprintf(fp, "%s", name_);

  for(int k = 0; k < (int)attributes_.size(); k++)
    fprintf(fp, " %s=\"%s\"",  getAttribute(k), values_[k]);

  if(start_ && end_)
    fprintf(fp, "/");
  else if(! start_ && ! end_)
    fprintf(fp, "?");

  fprintf(fp, ">\n");
}

//
// a helper class for optionally skipping loading of tags - 
// if tag name isn't listed, an empty tag is created and getName() returns NULL
// if tag attribute isn't listed, that attribute isn't saved
// this can save on heap thrash and fragmentation
//
TagInclusionTest::TagInclusionTest(const char *interesting_tags[], // only take tags with these names
                                   const char *interesting_attributes[]) // only take attributes with these names 
{
  for (int i=0;interesting_tags[i];i++) {
    tagNamesToInclude_.insert(interesting_tags[i]);
  }
  for (int j=0;interesting_attributes[j];j++) {
    tagAttributesToInclude_.insert(interesting_attributes[j]);
  }
}

TagInclusionTest::~TagInclusionTest() {
}

bool TagInclusionTest::includesName(const char *name, int namelen) const {
  char *end = (char *)name+namelen; // ugly const cast, but we'll restore before return
  char endval = *end;
  *end = 0;
  bool result = (tagNamesToInclude_.find(name)!=tagNamesToInclude_.end());
  *end = endval;
  return result;
}

bool TagInclusionTest::includesAttribute(const char *attr, int attrlen) const {
  char *end = (char *)attr+attrlen; // ugly const cast, but we'll restore before return
  char endval = *end;
  *end = 0;
  bool result = tagAttributesToInclude_.find(attr)!=tagAttributesToInclude_.end();
  *end = endval;
  return result;
}
