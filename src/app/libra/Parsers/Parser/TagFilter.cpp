/*

Program       : TagFilter                                                       
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

#include "TagFilter.h"

TagFilter::TagFilter(const char* name) {
  name_ = new char[strlen(name)+1];
  strcpy(name_, name);
  filter_ = False;
  counter_ = 0;
  num_ = -5; // any number
  done_ = False;
  attribute_names_ = new Array<char*>;
  attribute_vals_ = new Array<char*>;
}

TagFilter::TagFilter(const char* name, int num) {
  name_ = new char[strlen(name)+1];
  strcpy(name_, name);
  filter_ = False;
  counter_ = 0;
  num_ = num;
  done_ = False;
  attribute_names_ = new Array<char*>;
  attribute_vals_ = new Array<char*>;
}


TagFilter::~TagFilter() {
  if(attribute_names_ != NULL) {
    for(int k = 0; k < attribute_names_->length(); k++)
      if((*attribute_names_)[k] != NULL)
	delete [] (*attribute_names_)[k];
    delete attribute_names_;
  }
  if(attribute_vals_ != NULL) {
    for(int k = 0; k < attribute_vals_->length(); k++)
      if((*attribute_vals_)[k] != NULL)
	delete [] (*attribute_vals_)[k];
    delete attribute_vals_;
  }
  delete[] name_;
}


Boolean TagFilter::filter(Tag* tag) {
  if(num_ >= 0 && counter_ > num_)
    return False; // already done

  Boolean match = strstr(tag->getName(), name_) != NULL;
  if(tag->isStart()) {
    for(int k = 0; k < attribute_names_->length(); k++)
      if(tag->getAttributeValue((*attribute_names_)[k]) == NULL ||
       strcmp(tag->getAttributeValue((*attribute_names_)[k]), (*attribute_vals_)[k])) {
	match = False;
	k = attribute_names_->length();
      }
  } // only look for required attribute values if start tag

  if(done_) {
    filter_ = False; // done
    counter_++;
    done_ = False;
  }
  if(num_ == -1) { // special case, look for substring
    //cout << "here for " << name_ << endl;
    
    if(match) {
    
      //cout << "match of " << name_ << " within: " << tag->getName() << endl;
      if(tag->isStart()) {
	filter_ = True;
	  //done_ = False;
      }
      else if(! tag->isStart() && ! tag->isEnd()) { // special case tags
	filter_ = True;
	  //done_ = True;
      }
      //if(tag->isEnd())
      done_ = True;
    }
  } // num -1
  else if(match && ! strcmp(name_, tag->getName())) {
    if(tag->isStart()) {
      filter_ = True;
      done_ = False;
    }
    else if(! tag->isStart() && ! tag->isEnd()) { // special case tags
      filter_ = True;
      done_ = True;
    }
    if(tag->isEnd())
      done_ = True;
  }

  //   if(filter_) {
  //      cout << "true for " << tag->getName() << " with filter " << name_ << endl;
  //   }
  return filter_;
}

    
void TagFilter::enterRequiredAttributeVal(const char* attr_name, const char* attr_val) {
  char* next = new char[strlen(attr_name)+1];
  strcpy(next, attr_name);
  attribute_names_->insertAtEnd(next);
  next = new char[strlen(attr_val)+1];
  strcpy(next, attr_val);
  attribute_vals_->insertAtEnd(next);

}

    
