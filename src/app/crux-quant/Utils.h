#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include "io/SpectrumCollectionFactory.h"
#include "IndexedMassSpectralPeak.h"
#include "io/MatchFileReader.h"

using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

namespace CruxQuant {


string calcFormula(string seq);

const int NUMISOTOPES_REQUIRED = 2;  // May need to make this a user input
const int BINS_PER_DALTON = 100;

struct Identification {
    string Sequence;
    int Charge;
    double PeptideMass = 0.0;
    double MonoisotopicMass;
    double PeakfindingMass;
    double PrecursorCharge;
    
};
Crux::SpectrumCollection* loadSpectra(const string& file, int ms_level);

unordered_map<int, vector<CruxQuant::IndexedMassSpectralPeak>> indexedMassSpectralPeaks(Crux::SpectrumCollection* spectrum_collection);

vector<Identification> createIdentifications(MatchFileReader* matchFileReader);

unordered_map<string, vector<pair<double, double>>> calculateTheoreticalIsotopeDistributions(const vector<Identification>& allIdentifications);

void SetPeakFindingMass(vector<Identification>& allIdentifications, unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution);

vector<double> createChargeStates(const vector<Identification>& allIdentifications);



}  // namespace CruxQuant
