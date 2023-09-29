#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "io/SpectrumCollectionFactory.h"
#include "IndexedMassSpectralPeak.h"
#include "PpmTolerance.h"
#include "io/MatchFileReader.h"

using std::pair;
using std::string;
using std::unordered_map;
using std::map;
using std::vector;

namespace CruxQuant {

const int BINS_PER_DALTON = 100;
const double PROTONMASS = 1.007276466879;

const int NUMISOTOPES_REQUIRED = 2;  // May need to make this a user input
const double PEAK_FINDING_PPM_TOLERANCE = 20.0; // May need to make this a user input
const double PPM_TOLERANCE = 10.0; // May need to make this a user input


string calcFormula(string seq);

struct Identification {
    string Sequence;
    int Charge;
    double PeptideMass = 0.0;
    double MonoisotopicMass;
    double PeakfindingMass;
    double PrecursorCharge;
    string spectralFile;
};

struct Ms1ScanInfo{
    int OneBasedScanNumber;
    int ZeroBasedMs1ScanIndex;
    double RetentionTime;
};

struct IndexedSpectralResults{
    
    map<int, map<int, IndexedMassSpectralPeak>> _indexedPeaks;
    
    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;
};

Crux::SpectrumCollection* loadSpectra(const string& file, int ms_level);

IndexedSpectralResults indexedMassSpectralPeaks(Crux::SpectrumCollection* spectrum_collection, const string &spectra_file);

vector<Identification> createIdentifications(MatchFileReader* matchFileReader, const string &spectra_file);

unordered_map<string, vector<pair<double, double>>> calculateTheoreticalIsotopeDistributions(const vector<Identification>& allIdentifications);

void SetPeakFindingMass(vector<Identification>& allIdentifications, unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution);

vector<double> createChargeStates(const vector<Identification>& allIdentifications);

void quantifyMs2IdentifiedPeptides(string spectraFile, const vector<Identification>& allIdentifications);

double toMz(double mass, int charge);

double toMass(double massToChargeRatio, int charge);

IndexedMassSpectralPeak* getIndexedPeak(double theorMass, int zeroBasedScanIndex, PpmTolerance tolerance, int chargeState,  map<int, map<int, IndexedMassSpectralPeak>> indexedPeaks);

}  // namespace CruxQuant
