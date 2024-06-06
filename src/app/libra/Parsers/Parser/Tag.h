#ifndef TAG_H
#define TAG_H

/*
Program       : Tag
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: Tag.h 7816 2018-09-14 19:53:17Z dshteyn $

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

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <math.h>
#include <fstream>
#include <ostream>
#include "gzstream.h"
#include "Common/constants.h"
#include "Common/util.h"
#include "Common/Array.h"
using namespace std;

//
// a helper class for optionally skipping loading of tags - 
// if tag name isn't listed, an empty tag is created and getName() returns NULL
// if tag attribute isn't listed, that attribute isn't saved
// this can save on heap thrash and fragmentation
//

struct tagltstr {
  bool operator()(const char* s1, const char* s2) const {
    return strcmp(s1, s2)<0;
  }
};
struct tageqstr {
  bool operator()(const char* s1, const char* s2) const {
    return !strcmp(s1, s2);
  }
};
#ifdef _MSC_VER
#include <hash_set>
#ifdef _STLP_HASH_SET
typedef std::hash_set<const char *,std::hash<const char *>,tageqstr> tagname_set;
#else // assuming VC8's STL
typedef stdext::hash_set<const char *,stdext::hash_compare<const char *,tagltstr>> tagname_set;
#endif
#else
#include <ext/hash_set>
using namespace __gnu_cxx;
typedef hash_set<const char *,__gnu_cxx::hash<const char *>,tageqstr> tagname_set;
#endif

class TagInclusionTest {
 public:
  TagInclusionTest(const char *interesting_tags[], // only take tags with these names
		   const char *interesting_attributes[]); // only take attributes with these names
  ~TagInclusionTest();
  bool includesName(const char *name, int namelen) const;
  bool includesAttribute(const char *attr, int attrlen) const;
 private:
  tagname_set tagNamesToInclude_;
  tagname_set tagAttributesToInclude_;
};

//
// now a helper class for quickly finding the index of a value for an attribute
//

#include <vector>
#include <algorithm>
class charptr_int_pair {
 public:
  charptr_int_pair() {};
 charptr_int_pair(const char *charptr) : charptr_(charptr) {};
 charptr_int_pair(const char *charptr,int index) : charptr_(charptr),index_(index) {};
  bool operator == (const charptr_int_pair & rhs) const {
    return (index_==rhs.index_)&&(charptr_==rhs.charptr_);
  }
  const char *charptr_;
  int index_;
}; 
extern int compare_charptr_int_pairs(const void *a,const void *b); // in Tag.cxx, returns strcmp(a->charptr_,b->charptr_);

class charptr_int_map {
 public:
  charptr_int_map() {
    dirty_ = false;
  };
  void clear() {
    charptrs_.clear();
  }
  void trim() { // remove excess capacity
    charptrs_.trim();
  }
  int size() const {
    return charptrs_.size();
  }
  void add(const char *charptr,int index=-1) {
    dirty_ = true; // will need sorting
    charptrs_.insertAtEnd(charptr_int_pair(charptr,(index<0)?charptrs_.size():index));
  }
  const charptr_int_pair *find(const char *charptr) const {
    if (dirty_) { // needs sorting for efficient search
      charptrs_.sort(compare_charptr_int_pairs);
      dirty_ = false;
    }
    charptr_int_pair tmp(charptr);
    return charptrs_.find_sorted(&tmp,compare_charptr_int_pairs);
  }

  const char* getNthCharptr(int n) const {
    int i=charptrs_.size();
    if (i) {
      const charptr_int_pair* pair = &(charptrs_[0]);
      for (;i--;pair++) {
	if (n==pair->index_) {
	  return pair->charptr_;
	}
      }
    }
    return NULL;
  }

  const charptr_int_pair *getPair(int n) const {
    if (n < charptrs_.size()) {
      return &charptrs_[n];
    }
    return NULL;
  }

 private:
  mutable Array<charptr_int_pair> charptrs_;
  mutable bool dirty_;
};

//
// and now the Tag class itself
//
class Tag {

 public:

  Tag();
  Tag(const char* data,const TagInclusionTest* inclusionTest=NULL);
  Tag(const char* name, Boolean start, Boolean end);
  Tag(const char* name, Boolean start, Boolean end, const char* nspace);
  Tag(const Tag &rhs); // copy constructor
  ~Tag();
  Boolean isEnd() const {
    return end_;
  }
  Boolean isStart() const {
    return start_;
  }
  int getAttributeCount() const {
    return attributes_.size();
  }
  const char* getAttribute(int n) const {
    return attributes_.getNthCharptr(n);
  }
  const char* getAttributeValue(int n) const; 
  const char* getAttributeValue(const char* attr) const {
    const charptr_int_pair *pair = attributes_.find(attr);
    return pair?values_[pair->index_]:NULL;
  }

  const char* getName() const {
    return name_;
  }
  const char* getNameSpace() const {
    return namespace_;
  }
  void print() const;
  void write(ostream& os) const;
  //  template <class T> void writeTraditional(T * fp) const;
  void writeTraditional(ogzstream * fp) const;
  //  writeTraditional(fp);
  //}
  void writeTraditional(FILE * fp) const;
  //  writeTraditional(fp);
  //}
  void setAttributeValue(const char* attribute, const char* value); 
  // similar to setAttributeValue, but does not check for attribute existence
  // this can be more efficient when creating a new tag, since the attribute map 
  // doesn't have to be sorted for search
  void addAttributeValue(const char* attribute, const char* value);
  Boolean appendAttributeValue(const char* attribute, const char* value, Boolean unique);
  Tag* copy() const;
  Boolean isIdentical(const Tag &rhs) const; // compare
  void setEnd();
  void copyAttributesTo(Tag* output) const;
  void setLineNumber(int num) {
    line_number_ = num;
  }
  int getLineNumber() const {
    return line_number_;
  }

  void trim() { // minimize storage use by trimming excess vector capacity
    attributes_.trim();
    values_.trim();
  }

 protected:

  void setName(const char *newname); // for use by child classes
  void setNamespace(const char *newnamespace); // for use by child classes
  Boolean readAttribute(const char* data, int& index, const TagInclusionTest *test); // populate from text string
  Boolean start_; // whether start tag
  Boolean end_; // whether end tag

 protected:
  // storing attributes for fast lookup while preserving read order
  charptr_int_map attributes_;  // paired attribute name,read_order
  StringArray values_; // attribute values as values_[attributes_[n].second]

 private:
  const char* name_;
  const char* namespace_; // optional prefix
  int line_number_; // optional, useful when reading taglist files
  void insert_attribute(const char *attr, char *value); // add values with bookkeeping
  int getAttributeIndex(const char *attr) const {
    const charptr_int_pair *pair = attributes_.find(attr);
    return pair?pair->index_:-1;
  }
};


//
// convenience class for destroying tags on exit
//
class TagArray : public Array<Tag *> {
 public:
  ~TagArray() {
    for (int i=this->size();i--;) {
      delete this->arrPtr_[i];
    }
  }
  void clear() {
    for (int i=this->size();i--;) {
      delete this->arrPtr_[i];
    }
    Array<Tag *>::clear();
  }
  TagArray & operator = (const TagArray &rhs) {
    clear();
    for (int j=0; j<rhs.size();j++) {
      this->insertAtEnd(rhs.arrPtr_[j]->copy());
    }
    return *this;
  }

  // const element access
  const Tag & operator [] (int iIndex) const { 
    ARRAY_H_DEBUG_ASSERT((iIndex >= 0) && (iIndex < iUsed_));
    ARRAY_H_DEBUG_ASSERT(this->arrPtr_[iIndex]!=NULL);
    return *(this->arrPtr_[iIndex]);
  }

  // non-const element access
  Tag & operator [] (int iIndex) { 
    ARRAY_H_DEBUG_ASSERT((iIndex >= 0) && (iIndex < iUsed_));
    ARRAY_H_DEBUG_ASSERT(this->arrPtr_[iIndex]!=NULL);
    return *(this->arrPtr_[iIndex]);
  }

  void write(ostream& os) const {
    for(int k = 0; k < size(); k++) {
      arrPtr_[k]->write(os);
    }
  }
  void writeTraditional(FILE* fp) const {
    for(int k = 0; k < size(); k++) {
      arrPtr_[k]->writeTraditional(fp);
    }
  }
    
};

class TagArrayArray : public Array<TagArray*> {
 public:
  ~TagArrayArray() {
    for (int i=this->size();i--;) {
      delete this->arrPtr_[i];
    }
  }
};


#endif
