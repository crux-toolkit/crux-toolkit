#ifndef CONDITION2_H
#define CONDITION2_H

#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include "Parsers/Parser/Tag.h"
#include "Common/Array.h"

#define FALSE 0
#define TRUE 1
#define ddbug FALSE


class LibraConditionHandler2
{
private:

    std::ofstream m_fout;

public:

    const char * conditionFileName;

    /**
    * minimum integrated intensity of a reagent line to keep during outlier
    * removal stages in LibraProteinRatioParser package (default value is 0).
    */
    double minThresInten;   

    /**
    * name for output quantitation file (default is null)
    */
    std::string quantitationFileName;

    /**
    * number of reagent lines
    */
    int  numberOfReagentLines;

    /**
    * channel to use as reference channel in normalization
    */
    int normalizationChannel;


    /**
    *  number of pseudo-counts to use
    */
    unsigned int pseudoCount;


  
    LibraConditionHandler2();

    ~LibraConditionHandler2();

    void setConditionFileName( const char* );

    const char* getConditionFileName();

    void setMinimumThreshholdIntensity( double);

    double getMinimumThreshholdIntensity();

    void setQuantitationFileName( std::string );

    std::string getQuantitationFileName();

    void setNormalizationChannel( int );

    int  getNormalizationChannel();

    unsigned int  getPseudoCount();

    int readFile();

    int getNumberOfReagentLines();

};

#endif
