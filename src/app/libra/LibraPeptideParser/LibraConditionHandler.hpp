#ifndef CONDITION_HPP
#define CONDITION_HPP

#include "typedefs.hpp"

#include<iostream>
#include<fstream>
#include<vector>

#include "Parsers/Parser/Tag.h"
#include "Common/Array.h"


class LibraConditionHandler
{
private:

    const char* conditionFileName;

    std::ofstream m_fout;


    /**
    * Normalize intensities of the reagent mass lines using
    * value set in condition file.
    *   -2: for normalize against TIC (not recommended)
    *   -1: for normalize against the sum of the reagent lines
    *    0: for normalize against most intense peak (not recommended)
    *    1: for normalize against 1st listed reagent mass
    * ...n: for normalize against nth listed reagent mass
    * Note that the normalization in LibraProtein code, in other words,
    * the final protein quantitation, does not use the normalization here.
    * This normalization is just for user viewing in the interact.xml file.
    */
    int m_normalPosition;


public:

    /**
     * gets condition parameters into xml format
     * needed by TPP
     * @return xml tags in PepXML format
     */
    Array<Tag*>* getPepXMLTags();


    /**
    * reagent masses
    */
    std::vector<double> m_mass;


    /**
    * Matrix of isotopic contributions where indices are
    * contributingMz channel index and 
    * affected_channel index, and value is
    * affected_channel correction.
    * m_massIsotopes[contributingMz channel index - 1]
    * [affected_channel index - 1] = correction
    */
    vvf m_massIsotopes;


    bool m_useIsotopicCorrection;


    /**
    * mass tolerance around the m/z for each for each isotope
    */
    double m_tolerance;   


    /**
    * 1. math average 2.Weighted
    */
    int m_centroidingPref;


    /**
    *  number of times the centroiding process is repeated
    */
    int m_centroidingIteration;


    /**
    *  number of pseudo-counts to use and max channels to apply it to
    */
    unsigned int m_pseudoCount;
    unsigned int m_maxPseudoChannels;



    /**
    *
    */
    bool m_isToNormalize;

    /**
    * The MS level we want to look into
    */
    int m_requestedMsLevel;

    /**
    * Use scanNum (1) or RT (2)
    */
    int   m_outputPrefs;

    /**
    * number of quantitative reagents used
    */
    int m_nrReagents;

    /**
    * for Thermo SPS data, read reporter ions from corresponding MS3 scan given the MS/MS scan that gave ID
    */
    bool m_reporterFromMS3 = false;

    
    /**
    * starting from the position given by isoMassPos to the end
    */
    int zeroFill(int);

    LibraConditionHandler();
    ~LibraConditionHandler();

    void setFileName( const char* );

    const char* getFileName() const;

    int readFile();
    int getFromStdIn();
    int getRequestedMsLevel() const;
    int getCentroidingPref() const;
    int getIsToNormalize() const;

    const vvf &getMassIsotopes() const;
    const vf &getReagentMasses() const;

    int getNormalPosition() const;

    double getTolerance() const;

    bool getUseIsotopicCorrection() const;
    bool getReporterFromMS3() const;

    int writeLibraConditionHandlerToOutFile( const char* );
    int openOutFile(const char* );
    int closeOutFile();
    int getOutputPrefs() const;

    int getNumCentroidingIterations() const;

    int toString();

};

#endif
