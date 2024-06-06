//
// derived class for ASAPRatioPeptideParser regression test use
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
#include "ASAPRatioPeptideParserTagListComparator.h"
#include "Quantitation/XPress/XPressPeptideParser/XPressPeptideParserTagListComparator.h"
#include <math.h>

ASAPRatioPeptideParserTagListComparator::ASAPRatioPeptideParserTagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const char *rhsFilename) :
  TagListComparator(progname,why,lhs,rhsFilename) { // read from file then compare
}

int ASAPRatioPeptideParserTagListComparator::unacceptable_difference(const char *attribute, // return 1 if these values are too different for this attribute
								     const char *value1, const char *value2) {
  if (TagListComparator::unacceptable_difference(attribute,value1,value2) &&
      XPressPeptideParserTagListComparator::local_exceptions(attribute,value1,value2)) // in case xpress has already been invoked
    {
      return local_exceptions(attribute,value1,value2);
    }
  return 0;
}

// returns 0 iff it finds the differences acceptable
int ASAPRatioPeptideParserTagListComparator::local_exceptions(const char *attribute, // return 1 if these values are too different for this attribute
							      const char *value1, const char *value2) {
  // check for local exceptions due to differences in smoothing functions
  if ((!strcmp(attribute,"ratio")) ||
      (!strcmp(attribute,"error")) ||
      (!strcmp(attribute,"time_width"))) {
    double r1=atof(value1);
    double r2=atof(value2);
    if (r1 || r2) {
      double r= fabs((r2-r1)/Max(r1,r2));
      if (r <= .02) {
	return 0; // within a reasonable tolerance
      }
      if ((r <= .1) && !strcmp(attribute,"time_width")) {
	return 0; // within a reasonable tolerance
      }
    }
  } else if ((!strcmp(attribute,"area")) ||
	     (!strcmp(attribute,"area_error")) ||
	     (!strcmp(attribute,"background"))) {
    double r1=atof(value1);
    double r2=atof(value2);
    if (r1 || r2) {
      double r= fabs((r2-r1)/Max(r1,r2));
      if (r <= .05) {
	return 0; // within a reasonable tolerance
      }
    }
  } else if ((!strcmp(attribute,"is_heavy")) ||
	     (!strcmp(attribute,"is_light")) ||
	     (!strcmp(attribute,"left_valley")) ||
	     (!strcmp(attribute,"right_valley"))) {
    int scan1=atoi(value1);
    int scan2=atoi(value2);
    if (abs(scan2-scan1)<=2) {
      return 0; // peak center determination is nearly same
    }
  }

  return 1;
}
