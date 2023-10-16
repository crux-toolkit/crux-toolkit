#pragma once

#include <unordered_map>
#include <vector>

#include "CruxApplication.h"
#include "crux-quant/Utils.h"

using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

/**
 * \class CruxQuantApplication
 * \brief Application for quantifying peptides/proteins from MS/MS data
 */
class CruxQuantApplication : public CruxApplication {
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
};