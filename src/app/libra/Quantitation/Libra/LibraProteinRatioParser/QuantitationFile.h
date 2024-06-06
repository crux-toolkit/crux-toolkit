#ifndef QFILE_H
#define QFILE_H

#include <cstdlib> 
#include<iostream>
#include<fstream>
#include<string>
#include <vector>

using std::ios;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;
using std::fstream;


class QuantitationFile
{

public:

    std::ofstream m_fout;

    /** 
    * name of output file to write quantitation too
    */
    string fileName;

    /**
    * open existing file
    * @return 0 for successful
    */
    //int openExistingFile();

    /**
    * open file
    * @return 0 for successful
    */
    int initFile();

    /**
    * close file
    * @return 0 for successful
    */
    //int closeFile();

    /**
    * create new empty file
    */
    void init();

    // xxxxxxx the following should be stored in separate object one day,
    // but need to split parsers into read write classes, and use a protein 
    // quantitation object to store the protein and it;s peptide info.
    // that's where the peptide removal should be done too.

    // Expecting the following format to try to make it usable input for user in Excel:
    //protein peptide nr1 nr2 nr3 nr4 in1 in2 in3 in4 is_rejected
    //protein peptide nr1 nr2 nr3 nr4 in1 in2 in3 in4
    //protein peptide nr1 nr2 nr3 nr4 in1 in2 in3 in4
    //protiein \t     \t  \t  \t  \t  \t  \t  \t  \t  \t r1 r2 r3 r4
    //protiein \t     \t  \t  \t  \t  \t  \t  \t  \t  \t e1 e2 e3 e4
    string protein_name;

    /**
    * expecting entries to be of format: 
    * peptide<tab>nr1<tab>nr2<tab>nr3<tab>nr4<tab>in1<tab>in2<tab>in3<tab>in4<tab>is_rejected
    */
    vector<string>* peptide_quantitation_lines;

    /**
    * protein 2Xn_channels of tabs r1\tr2\tr3\tr4
    * protein 2Xn_channels of tabs e1\te2\te3\te4
    */
    vector<string>* protein_quantitation_lines;

    /** 
    * append new lines to file
    * @param string line to append to file
    * @return 0 for successful
    */
    int writeStringVectorToOutfile( vector<string>* );

    /**
    * number of reagent lines
    */
    int num_reagent_lines;


    QuantitationFile(string fname);

    ~QuantitationFile();

    /** 
    * append new line to file
    * @param string line to append to file
    * @return 0 for successful
    */
    int writeStringToOutfile( string );


    /** 
    * get file name
    * @return file name
    */
    string getFileName( );

    /** 
    * set file name
    * @param file name
    */
    void setFileName( string );

    /** 
    * set latest protein name
    * @param latest protein name
    */
    void setProteinName( string );

    /**
    * set line of peptide information to be written to quantitation file,
    * expecting entries to be of format: 
    * peptide<tab>nr1<tab>nr2<tab>nr3<tab>nr4<tab>in1<tab>in2<tab>in3<tab>in4<tab>is_rejected
    */
    void setPeptideLine( string );

    /**
    * set line of protein information to be written to quantitation file
    * protein 2Xn_channels of tabs r1\tr2\tr3\tr4
    * protein 2Xn_channels of tabs e1\te2\te3\te4
    */
    void setProteinLine( string );


    /**
    * write quantitation information to outfile
    */
    void writeQuantitationToOutfile();


    // clear protein name, peptide_quantitation_lines, and protein_quantitation_lines
    void reset();

    /**
    * set number of reagent lines
    * @param number of reagent lines
    */
    void setNumReagentLines( int );

    void allocateVectors();

    void deleteVectors();
};

#endif
