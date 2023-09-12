#pragma once

#include <vector>
#include <unordered_map>
#include <list>
#include "CruxApplication.h"
#include "io/SpectrumCollectionFactory.h"
#include "io/MatchFileReader.h"
#include "crux-quant/IndexedMassSpectralPeak.h"


using std::string;
using std::vector;
using std::pair;
using std::list;
using std::unordered_map;

/**
 * \class CruxQuantApplication
 * \brief Application for quantifying peptides/proteins from MS/MS data
 */
class CruxQuantApplication: public CruxApplication{
    public:

        static const int BinsPerDalton = 100;
        /**
         * Constructor
         */
        CruxQuantApplication();

        /**
         * Destructor
         */
        ~CruxQuantApplication();

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
        virtual vector< pair<string, string> > getOutputs() const;

        /**
         * \returns the enum of the application, default MISC_COMMAND
         */
        virtual COMMAND_T getCommand() const;

        /**
         * \returns whether the application needs the output directory or not. (default false)
         */
        virtual bool needsOutputDirectory() const;

        virtual void processParams();

        static Crux::SpectrumCollection* loadSpectra(const string& file, int ms_level);

        static unordered_map<int, list<CruxQuant::IndexedMassSpectralPeak>> indexedMassSpectralPeaks(Crux::SpectrumCollection* spectrum_collection);

        static unordered_map<string, list<pair<double, double>>> calculateTheoreticalIsotopeDistributions(MatchFileReader* matchFileReader);

        static double getMinChargeState(MatchFileReader* matchFileReader);

        static double getMaxChargeState(MatchFileReader* matchFileReader);

        static vector<double> CreateChargeStates(double minChargeState, double maxChargeState);
};