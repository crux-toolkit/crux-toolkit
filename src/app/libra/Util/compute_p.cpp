/*

Program       : compute_ttest_p
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>
Date          : 03.22.2018
SVN Info      : $Id$

Copyright (C) 2018 David Shteynberg

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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA3

*/
#include "Common/util.h"
#include "Common/TPPVersion.h"
#include <iostream>
#ifndef _LGPL_
#include <gsl/gsl_cdf.h>
using namespace std;
int main(int argc, char** argv) {
  //TODO: Add input error handling
  if(argc < 3) {
    //    cerr <<  argv[0] << " (" << szTPPVersionInfo << ")" << endl;
    cerr << "USAGE: " << argv[0] << " <TTEST|CHISQ|FTEST> <value> <degrees_of_freedom1> (<degrees_of_freedom2>)  " << endl;
    return 1;
  }

    double val = atof(argv[2]);
    double df = atof(argv[3]);
    double df2 = df;
    if (argc > 4) {
      df2  = atof(argv[4]);
    }
    double p = -1;
    if (!strcmp(argv[1], "FTEST")) {
      p = gsl_cdf_fdist_P(val, df, df2);
    } 
    else if (!strcmp(argv[1], "TTEST")) {
      p = gsl_cdf_tdist_P(val, df);
    } 
    else if (!strcmp(argv[1], "CHISQ")) {
      p = gsl_cdf_chisq_P(val, df);
    }

    cout << p << endl;

    return 0;

}
#endif
