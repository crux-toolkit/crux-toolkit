//
// derived class for XpressPeptideParser regression test use
//
// Copyright (C) Insilicos LLC 2005 All Rights Reserverd
/*
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
*/

class XPressPeptideParserTagListComparator : public TagListComparator {
public:
	XPressPeptideParserTagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const char *rhsFilename); // read from file then compare
   // returns 0 iff it finds the differences acceptable  
   static int local_exceptions(const char *attribute, // return 1 if these values are too different for this attribute
                             const char *value1, const char *value2);
protected:
   virtual int unacceptable_difference(const char *attribute, // return 1 if these values are too different for this attribute
       const char *value1, const char *value2);
};

