#include "EnzymeSpecificity.h"

/*

Program       : EnzymeSpecificity for PeptideProphet                                                       
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

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

EnzymeSpecificity::EnzymeSpecificity() { }


// factory for producing enzyme digestions
// register all new enzyme digestions here
EnzymeDigestion* EnzymeSpecificity::getEnzymeDigestion(char* enz) {
  if(enz == NULL) // default
    return new TrypticEnzymeDigestion();
  if(strcmp(enz, "tryptic") == 0)
    return new TrypticEnzymeDigestion();
  if(strcmp(enz, "gluC") == 0)
    return new GluCEnzymeDigestion();
  if(strcmp(enz, "gluC_bicarb") == 0)
    return new GluC_bicarbEnzymeDigestion();
  if(strcmp(enz, "chymotryptic") == 0)
    return new ChymotrypticEnzymeDigestion();
  if(strcmp(enz, "elastase") == 0)
    return new Elastase();
  if(strcmp(enz, "nonspecific") == 0) 
    return new NonspecificEnzymeDigestion();
  if(strcmp(enz, "tca") == 0)
    return new TrypChymAspnEnzymeDigestion();
  if(strcmp(enz, "CNBr") == 0)
    return new CNBrEnzymeDigestion();
  if(strcmp(enz, "AspN") == 0)
    return new AspNEnzymeDigestion();
  if(strcmp(enz, "tryptic/CNBr") == 0)
    return new TrypticCNBrEnzymeDigestion();


  // enzyme not registered
  return NULL;
}
