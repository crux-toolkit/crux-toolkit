// derived class for ASAPRatioPvalueParser regression test use
//
// Copyright (C) Insilicos LLC 2005 All Rights Reserved

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
#include "ASAPRatioPvalueParserTagListComparator.h"
#include <math.h>

ASAPRatioPvalueParserTagListComparator::ASAPRatioPvalueParserTagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const char *rhsFilename) :
  TagListComparator(progname,why,lhs,rhsFilename) { // read from file then compare
}

static int checktoler(const char *value1, const char *value2, double toler) {
  double r1=atof(value1);
  double r2=atof(value2);
  if (r1 || r2) {
    double r= fabs((r2-r1)/Max(r1,r2));
    if (r <= toler) {
      return 0; // within a reasonable tolerance
    }
  }
  return 1;
}

int ASAPRatioPvalueParserTagListComparator::unacceptable_difference(const char *attribute, // return 1 if these values are too different for this attribute
								    const char *value1, const char *value2) {
  if (TagListComparator::unacceptable_difference(attribute,value1,value2)) {
    // check for local exceptions due to differences in smoothing functions
    if (strstr(attribute,"adj_ratio_standard_dev")) {
      return checktoler(value1,value2,.4);
    } else if ((!strcmp(attribute,"pvalue")) ||
	       (!strcmp(attribute,"decimal_pvalue"))) {
      // given that this is an iterative calculation with a 
      // random seed, be generous in tolerance especially with smaller values
      double r1=atof(value1);
      double r2=atof(value2);
      return checktoler(value1,value2,(Min(r1,r2)<.01)?.7:.5);
    } else if (strstr(attribute,"adj_ratio_mean")) {
      // heavy rand() influence here, be generous with tolerance
      return checktoler(value1,value2,.02);
    } else if ((!strcmp(attribute,"background_ratio_mean")) ||
	       (!strcmp(attribute,"background_ratio_stdev")) ||
	       (!strcmp(attribute,"background_fitting_error"))) {
      // heavy rand() influence here, be generous with tolerance
      return checktoler(value1,value2,.4);
    }
    return 1;
  }
  return 0;
}
