/*
Program       : QuantitationFile
Author        : Nichole King
Date          : 10.16.05
SVN info      : $Id: QuantitationFile.cpp 8873 2023-02-27 08:14:30Z real_procopio $


Copyright (C) 2005 Nichole King

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

Nichole King
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
*/

#include "QuantitationFile.h"


/**
* constructor sets output file name 
*/
QuantitationFile::QuantitationFile(string fname) {
  setFileName(fname);
  cout<<"QuantitationFile constructed" << endl;
};


/**
*/
QuantitationFile::~QuantitationFile() {
  deleteVectors();
  cout<<"QuantitationFile destroyed" << endl;
};


/**
* get putput file name
* @return fname output file name
*/
string QuantitationFile::getFileName() {
  return fileName;
}

/**
* set putput file name
* @param fname output file name
*/
void QuantitationFile::setFileName(string fname) {
  fileName = fname;
};


/**
* create empty file and close it
*/
void QuantitationFile::init() {
  allocateVectors();
  initFile();
}


void QuantitationFile::allocateVectors() {
  try {
    peptide_quantitation_lines = new vector<string>;
    protein_quantitation_lines = new vector<string>;

  } catch (std::bad_alloc xa) {
    cout << "Couldn't allocate enough memory for vectors" << endl;
    exit(0);
  };

} // end allocateVectors



void QuantitationFile::deleteVectors() {
  delete peptide_quantitation_lines;
  delete protein_quantitation_lines;
} // end deleteVectors


/**
* open new output file, write a header line, and close file
* @return 0 for successful
*/
int QuantitationFile::initFile() {
  int successful = 1;
  std::ofstream out ( fileName.c_str() );

  cout << "opening file: " << fileName.c_str() << endl;

  try {
    if (!out ) throw 100; 
    successful = 0;

  } catch (int code) {
    if (code == 100 )
      cerr << "Error: could not open " << fileName << endl;
  }

  out.close();
  cout << "closed file: " << fileName.c_str() << endl;

  return successful;
}


/**
* open existing output file
* @return 0 for successful
*/
/*
int QuantitationFile::openExistingFile()
{

    int successful = 1;

    try
    {

        if (!fileName.size() > 0) throw 99;

        m_fout.open( fileName.c_str(), ios::app );

        if (!m_fout ) throw 100;

        successful = 0;

    } catch (int code)
    {

        if (code == 99 )
            cerr << "Error: quantitation file name not set "<< endl;

        if (code == 100 )
            cerr << "Error: could not open " << fileName << endl;

    }

    return successful;


}
*/


/**
* close output file
* @return 0 for successful
*/
/*
int QuantitationFile::closeFile()
{

    int successful = 1;

    if ( !m_fout.is_open() )
    {

        cerr << "Error: cannot close the already closed file " 
        << fileName << endl;

    } else
    {

        m_fout.close();

        successful = 0;

    }

    return successful;

}
*/

/**
* append string to existing outfile (adds end-of-line too)
* @param string to write 
* @return 0 for successful
*/
int QuantitationFile::writeStringToOutfile(string line) {
  int successful = 1;

  try {
    if (!(fileName.size() > 0)) throw 99;

    m_fout.open( fileName.c_str(), ios::app );

    if (!m_fout ) throw 100;

    m_fout << line << endl;
    successful = 0;

  } catch (int code) {
    if (code == 99 )
      cerr << "Error: quantitation file name not set "<< endl;

    if (code == 100 )
      cerr << "Error: could not open " << fileName << endl;

  }

  m_fout.close();
  return successful;
}


/**
* append strings to existing outfile (adds end-of-line too)
* @param string vector to write 
* @return 0 for successful
*/
int QuantitationFile::writeStringVectorToOutfile(vector<string>* stringVector) {
  int successful = 1;

  try {
    if (!(fileName.size() > 0)) throw 99;

    m_fout.open( fileName.c_str(), ios::app );

    if (!m_fout ) throw 100;

    for (int i=0; i < (int) (*stringVector).size(); i++) {
      m_fout << (*stringVector)[i] << endl;
    }

    successful = 0;

  } catch (int code) {
    if (code == 99 )
      cerr << "Error: quantitation file name not set "<< endl;

    if (code == 100 )
      cerr << "Error: could not open " << fileName << endl;
  }

  m_fout.close();
  return successful;
}


/**
* set protein name
* @param protein name
*/
void QuantitationFile::setProteinName( string protein) {
  cout << "set protein:" <<endl;
  cout << protein_name <<endl;
  protein_name = protein;
}


/**
* enter a string into  peptide_quantitation_lines.
* expecting entries to be of format: 
* peptide<tab>nr1<tab>nr2<tab>nr3<tab>nr4<tab>in1<tab>in2<tab>in3<tab>in4<tab>is_rejected
* this puts a prefix of protein<tab> before storing in structure
* @param line of peptide quantition information
*/
void QuantitationFile::setPeptideLine( string line) {
  line = protein_name + string("\t") + line;
  (*peptide_quantitation_lines).push_back(line);
}


/**
* set number of reagent lines
*/
void QuantitationFile::setNumReagentLines( int numR) {
  num_reagent_lines = numR;
}


/**
* enter a string into  protein_quantitation_lines.
* expecting entries to be of format: 
* r1<tab>r2<tab>r3<tab>r4<tab>
* or
* e1<tab>e2<tab>e3<tab>e4<tab>
* This puts a prefix of protein<tab> and 2X n_channels of <tab>'s and another <tab>
* @param line of protein quantition information
*/
void QuantitationFile::setProteinLine( string line) {
  string prefix = protein_name;

  for (int i =1; i <= num_reagent_lines; i++) {
    prefix = prefix + string("\t");
  }

  //tab placeholder for is_rejected column:
  prefix = prefix + string("\t");

  line = prefix + string("\t") + line;

  (*protein_quantitation_lines).push_back(line);
}


/**
*clear protein name, peptide_quantitation_lines, and protein_quantitation_lines
*/
void QuantitationFile::reset() {
  cout << "IN RESET " << endl;

  writeQuantitationToOutfile();

  protein_name.erase();
  (*peptide_quantitation_lines).clear();
  (*protein_quantitation_lines).clear();
}


/**
* writes protein quantitation info to outfile.
* writing format: (this is the most up to date note, so overrides other specs in this file)
*  protein peptide\tnr1\tnr2\tnr3t\nr4\tin1\in2\tin3\tin4\kept?
*  protein peptide\tnr1\tnr2\tnr3t\nr4\tin1\in2\tin3\tin4\tkept?
*  protein peptide\tnr1\tnr2\tnr3t\nr4\tin1\in2\tin3\tin4\tkept?
*  protein 2Xn_channels+1 of tabs r1\tr2\tr3\tr4\te1\te2\te3\te4
*  
*/
void QuantitationFile::writeQuantitationToOutfile() {
  cout << "in write quantitation" << endl;
  cout << "protein is : " <<  endl;
  cout << protein_name <<  endl;

  writeStringVectorToOutfile(peptide_quantitation_lines);
  writeStringVectorToOutfile(protein_quantitation_lines);
}
