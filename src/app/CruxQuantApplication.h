#pragma once

#include "CruxApplication.h"
#include "crux-quant/IndexedMassSpectralPeak.h"
#include "crux-quant/PeakIndexingEngine.h"

using std::string;
using std::vector;
using std::pair;

/**
 * \class CruxQuantApplication
 * \brief Application for quantifying peptides/proteins from MS/MS data
 */
class CruxQuantApplication: public CruxApplication{
    public:
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
        virtual vector<std::string> getArgs() const;

        /**
         * \returns the command options
         */
        virtual vector<std::string> getOptions() const;

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
};