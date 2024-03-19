#pragma once

#include <unordered_map>
#include <vector>

#include "CruxApplication.h"
#include "crux-lfq/Utils.h"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/SpectrumListWrapper.hpp"

using CruxLFQ::BINS_PER_DALTON;
using CruxLFQ::getScanID;
using CruxLFQ::Identification;
using CruxLFQ::IndexedMassSpectralPeak;
using CruxLFQ::IndexedSpectralResults;
using CruxLFQ::Ms1ScanInfo;
using CruxLFQ::PSM;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

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

    static pwiz::msdata::SpectrumListPtr loadSpectra(const string& file, int ms_level);

    static IndexedSpectralResults indexedMassSpectralPeaks(SpectrumListPtr spectrum_collection, const string& spectra_file);

    static vector<Identification> createIdentifications(const vector<PSM>& psm_data, const string& spectra_file);
};