// -*- mode: c++ -*-


/*
    File: SimpleXMLWriter.h
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



#ifndef _INCLUDED_SIMPLEXMLWRITER_H_
#define _INCLUDED_SIMPLEXMLWRITER_H_

#include <string>
#include <stack>
#include <vector>
#include <iostream>

class SimpleXMLWriter {
public:
	SimpleXMLWriter();
	virtual ~SimpleXMLWriter();

	void setOutputStream( std::ostream& o ) { pOut_ = &o; }

	void startDocument(void);

	void open(const std::string& tagname);

	void open(const std::string& tagname, const std::vector< std::pair<std::string, std::string> > & attrlist);

	template<typename StreamableType>
	void attr(const std::string& attrname, const StreamableType& val) {
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
		(*pOut_) << attrname << "=\"" << val << "\"";
	}

#ifdef _MSC_VER // MSVC won't use this specialization unless it's in the header
	template<>
	inline void attr(const std::string& attrname, const std::string& val) {
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

	void attr(const std::vector< std::pair< std::string, std::string> > & attrlist);

	void noattr(void);

	void data(const std::string& data);

	void close();

	void closeAll();

	bool condenseAttr_;

protected:

	std::ostream* pOut_; // must be set before use.  Add error checking for unset case.
	void setIndentStr();

	int indent_;
	bool tagOpen_;
	bool hasAttr_;
	bool hasData_;

	std::stack<std::string> tags_;
	std::string indentStr_;
	std::string spaceStr_;
};


#endif // _INCLUDED_SIMPLEXMLWRITER_H_
