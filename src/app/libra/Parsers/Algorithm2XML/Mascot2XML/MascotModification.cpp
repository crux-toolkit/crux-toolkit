#include "MascotModification.h"

/*

Program       : Mascot2XML converter to pepXML
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

using namespace std;

MascotModification::MascotModification(char* modname, Boolean variable, int index) { 
  cerr << "error: stub function called" << endl;
  exit(1);
}





MascotModification::MascotModification(char residue, double mass1, double massdiff1, Boolean variable1, int index1) { 
  if(residue == '1') { // protein n term
    aa = 'n';
    prot_term = True;
  }
  else if(residue == '2') { // protein n term
    aa = 'c';
    prot_term = True;
  }
  else {
    aa = residue;
    prot_term = False;
  }
  //  cout << "adding mod with aa: " << aa << endl;
  mass = mass1;
  massdiff = massdiff1;
  variable = variable1;
  n_terminal = false;
  c_terminal = false;
  index = index1;
}

MascotModification::MascotModification(char residue, double mass1, double massdiff1, Boolean variable1, int index1, Boolean n_term, Boolean c_term) { 
  if(residue == '1') { // protein n term
    aa = 'n';
    prot_term = True;
  }
  else if(residue == '2') { // protein n term
    aa = 'c';
    prot_term = True;
  }
  else {
    aa = residue;
    prot_term = False;
  }
  //  cout << "adding mod with aa: " << aa << endl;
  mass = mass1;
  massdiff = massdiff1;
  variable = variable1;
  n_terminal = n_term;
  c_terminal = c_term;
  index = index1;
}

void MascotModification::print() {
  cout << "aa: " << aa;
  if(prot_term)
    cout << "*";
  cout << " mass: " << mass << " index: " << index;
  if(variable)
    cout << " variable" << endl;
  else
    cout << " static" << endl;

  /*
  if(variable)
    cout << "aa: " << aa << " mass: " << mass << " index: " << index <<  " variable" << endl;
  else
    cout << "aa: " << aa << " mass: " << mass << " index: " << index << " static" << endl;
  */
}
