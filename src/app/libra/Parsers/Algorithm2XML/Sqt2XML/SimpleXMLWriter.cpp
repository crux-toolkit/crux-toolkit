// -*- mode: c++ -*-


/*
    File: SimpleXMLWriter.cpp
    Description: Basic stack-based XML writer.
    Date: July 25, 2007

    Copyright (C) 2007 Natalie Tasman, ISB Seattle


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
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/


#include "SimpleXMLWriter.h"

#include <sstream>
#include <iostream>

using namespace std;

SimpleXMLWriter::SimpleXMLWriter() : 
condenseAttr_(false),
pOut_(NULL),
indent_(0),
tagOpen_(false),
hasAttr_(false),
hasData_(false),
indentStr_(""), 
spaceStr_(" ")
{
}


SimpleXMLWriter::~SimpleXMLWriter() {
}


void 
SimpleXMLWriter::startDocument(void) {
	(*pOut_) << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
}


void 
SimpleXMLWriter::open(const string& tagname) {
	if (tagOpen_) {
		if (hasAttr_) {
			(*pOut_) << " ";
		}
		(*pOut_) << ">";
		tagOpen_ = false;
	} 
	if (!tags_.empty()) {
		(*pOut_) << '\n';
	}
	tags_.push(tagname);

	tagOpen_ = true;
	indent_ += 1;
	setIndentStr();
	(*pOut_) << indentStr_ << "<" << tagname;
	hasAttr_ = false;
	hasData_ = false;
}

#ifndef _MSC_VER
template<>
inline void SimpleXMLWriter::attr(const std::string& attrname, const std::string& val) {
	if (!hasAttr_) {
	  // whitespace formatting to assist parsing:
	  // first attribute should appear after element name,
	  // separated by one space character, on same line
	  // as the element start.
	  (*pOut_) << " ";
	  hasAttr_ = true;
	}
	else if (!condenseAttr_) { 
	(*pOut_) << '\n' << indentStr_ << spaceStr_;
	}
	else {
	  (*pOut_) << spaceStr_;
	}
	(*pOut_) << attrname << "=\"";
	for( size_t i=0; i < val.length(); ++i )
		switch( val[i] )
		{
			case '"':	(*pOut_) << "&quot;"; break;
			case '\'':	(*pOut_) << "&apos;"; break;
			case '<':	(*pOut_) << "&lt;"; break;
			case '>':	(*pOut_) << "&gt;"; break;
			case '&':	(*pOut_) << "&amp;"; break;
			default:	(*pOut_) << val[i]; break;
		}
	(*pOut_) << "\"";
}
#endif

void 
SimpleXMLWriter::attr(const vector< pair< string, string> > & attrlist) {
	vector< pair< string, string> >::size_type numAttr  = attrlist.size();
	if (numAttr > 0) {
		hasAttr_ = true;
		for (vector< pair< string,string> >::size_type i=0; i<attrlist.size(); ++i) {
			(*pOut_) << " " << attrlist[i].first << "=" << "\"" << attrlist[i].second << "\"";
		}
	}
}


void
SimpleXMLWriter::noattr(void) {
	if (tagOpen_) {
		(*pOut_) << ">";
		tagOpen_ = false;
	} 
}

void 
SimpleXMLWriter::data(const string& data) {
	if (tagOpen_) {
		if (hasAttr_) {
			(*pOut_) << " ";
		}
		(*pOut_) << ">";
		tagOpen_ = false;
	} 
	(*pOut_) << data;
	hasData_ = true;
}



void 
SimpleXMLWriter::close() {
	if (tagOpen_) {
		(*pOut_) << " />";
	} else {
		if (!hasData_) {
			(*pOut_) << '\n' << indentStr_;
		}
		(*pOut_) << "</" << tags_.top() << ">";
	}

	tagOpen_ = false;
	hasData_ = false;

	--indent_;
	setIndentStr();

	tags_.pop();
	if (tags_.empty()) {
		(*pOut_) << '\n';
	}
}



void 
SimpleXMLWriter::closeAll() {
	while (!tags_.empty()) {
		this->close();
	}
}





void 
SimpleXMLWriter::setIndentStr() {
	indentStr_ = "";
	for (int i=1; i<indent_; i++) {
		indentStr_ += spaceStr_;
	}
}
