#ifndef LIBRASUMMARY_HPP
#define LIBRASUMMARY_HPP

#include "LibraConditionHandler.hpp"

#include<vector>
#include<utility>

using std::vector;
using std::pair;
using std::make_pair;


class LibraSummary
{
private:
    LibraConditionHandler* pLibraConditionHandler;

public:

    LibraSummary( LibraConditionHandler* );

    ~LibraSummary();

    double getMassTolerance();

    bool getCentroiding();

    int getCentroidingPref();

    int getNumCentroidingIterations();

    int getNormalization();

    int getOutputType();

    vector<double> getFragmentMassValues( int targetID);

};

#endif
