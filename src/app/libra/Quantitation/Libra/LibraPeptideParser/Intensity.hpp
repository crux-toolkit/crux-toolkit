#ifndef INTENSITY_HPP
#define INTENSITY_HPP

#include "LibraConditionHandler.hpp"
#include "Quantitation.hpp"
#include "typedefs.hpp"

#include<vector>
#include<utility>

using std::vector;
using std::pair;
using std::make_pair;

class Intensity
{
private:

    Quantitation* pQuantitation;

    LibraConditionHandler* pLibraConditionHandler;

    int scanNum;

    char* mzXMLFileName;

public:

    Intensity();
  
    Intensity(LibraConditionHandler*, char*, mzParser::RAMPFILE*, int scanNum);

    ~Intensity();

    vector<double> getTargetMass();

    vector<double> getAbsoluteIntensities();

    vector<double> getNormalizedIntensities();

    int getNumberOfChannels();

    int getScanNumber();

};

#endif
