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

#include "Parsers/Parser/TagListComparator.h"
#include "XPressPeptideParserTagListComparator.h"
#include "Quantitation/ASAPRatio/ASAPRatioPeptideParser/ASAPRatioPeptideParserTagListComparator.h"
#include <math.h>

XPressPeptideParserTagListComparator::XPressPeptideParserTagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const char *rhsFilename) :
  TagListComparator(progname,why,lhs,rhsFilename) { // read from file then compare
}

int XPressPeptideParserTagListComparator::unacceptable_difference(const char *attribute, // return 1 if these values are too different for this attribute
                            const char *value1, const char *value2) {
   if (TagListComparator::unacceptable_difference(attribute,value1,value2) &&
       ASAPRatioPeptideParserTagListComparator::local_exceptions(attribute,value1,value2)) // in case ASAPRatio has already been invoked
      {
      // check for local exceptions
      return local_exceptions(attribute,value1,value2);
   }
   return 0;
}

// returns 0 iff it finds the differences acceptable
int XPressPeptideParserTagListComparator::local_exceptions(const char *attribute, // return 1 if these values are too different for this attribute
                                                           const char *value1, const char *value2) {
   // check for local exceptions
   if (!strcmp(attribute,"heavy2light_ratio")) {
      const char *r1=strchr(value1,':');
      const char *r2=strchr(value2,':');
      if (r1 && r2) {
         double v1= atof(++r1);
         double v2= atof(++r2);
         if ((fabs(v1-v2)/Min(v1,v2)) <= .01) {
            return 0; // within a reasonable tolerance
         }
      }
   }
   
   return 1;
}
