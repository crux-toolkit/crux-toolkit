#pragma once

#include <unordered_map>
#include <vector>

#include "CruxApplication.h"


#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumListWrapper.hpp"

#include "crux-lfq/Utils.h"

using std::pair;
using std::string;
using std::unordered_map;
using std::vector;
using CruxLFQ::PSM;
using CruxLFQ::IndexedSpectralResults;
using CruxLFQ::Identification;
using CruxLFQ::IndexedMassSpectralPeak;
using CruxLFQ::Ms1ScanInfo;
using CruxLFQ::BINS_PER_DALTON;
using CruxLFQ::getScanID;

typedef pwiz::msdata::SpectrumListPtr SpectrumListPtr;

/**
 * \class CruxLFQApplication
 * \brief Application for quantifying peptides/proteins from MS/MS data
 */
class CruxLFQApplication : public CruxApplication {
 public:
    /**
     * Constructor
     */
    CruxLFQApplication();

    /**
     * Destructor
     */
    ~CruxLFQApplication();

    /**
     * Main method
     */
    virtual int main(int argc, char** argv);

    int main(const string& psm_file, const vector<string>& spec_files);

    /**
     * \returns the name of the subclassed application
     */
    virtual string getName() const;

    /**
     * \returns the description of the subclassed application
     */
    virtual string getDescription() const;

    /**
     * \returns the command arguments
     */
    virtual vector<string> getArgs() const;

    /**
     * \returns the command options
     */
    virtual vector<string> getOptions() const;

    /**
     * \returns the command outputs
     */
    virtual vector<pair<string, string> > getOutputs() const;

    /**
     * \returns the enum of the application, default MISC_COMMAND
     */
    virtual COMMAND_T getCommand() const;

    /**
     * \returns whether the application needs the output directory or not. (default false)
     */
    virtual bool needsOutputDirectory() const;

    virtual void processParams();

    map<int, PSM> create_psm_map(const string& psm_file);

    pwiz::msdata::SpectrumListPtr loadSpectra(const string& file, int ms_level);

    IndexedSpectralResults indexedMassSpectralPeaks(SpectrumListPtr spectrum_collection, const string& spectra_file);

    vector<Identification> createIdentifications(const map<int, PSM>& psm_datum, const string& spectra_file, SpectrumListPtr spectrum_collection);
};