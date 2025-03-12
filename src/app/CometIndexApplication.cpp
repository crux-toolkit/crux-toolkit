/**
 * \file CometIndexApplication.cpp
 * \brief Runs comet -i
 *****************************************************************************/
#include "CometIndexApplication.h"
#include "io/DelimitedFile.h"
#include "io/DelimitedFileWriter.h"
#include "util/AminoAcidUtil.h"
#include "util/CarpStreamBuf.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"

using namespace std;

/**
 * \returns a blank CometApplication object
 */
CometIndexApplication::CometIndexApplication() {
}

/**
 * Destructor
 */
CometIndexApplication::~CometIndexApplication() {
}

/**
 * main method for CometIndexApplication
 */
int CometIndexApplication::main(int argc, char **argv) {
    return main();
}

int CometIndexApplication::main() {

    /* Re-route stderr to log file */
    CarpStreamBuf buffer;
    streambuf *old = std::cerr.rdbuf();
    std::cerr.rdbuf(&buffer);

    /* set parameters */
    vector<string> input_files;
    vector<InputFileInfo *> pv_input_files;
    setCometParameters(input_files, pv_input_files);
    searchManager_.AddInputFiles(pv_input_files);

    bool success = false;
    // Create fragment index or run search
    success = searchManager_.CreateFragmentIndex();

    /* Recover stderr */
    std::cerr.rdbuf(old);

    return success ? 0 : 1;
}





/**
 * \returns the command name for CometApplication
 */
string CometIndexApplication::getName() const { return "comet-index"; }

/**
 * \returns the description for CometIndexApplication
 */
string CometIndexApplication::getDescription() const {
    return "[[nohtml:Create an index of fragment ions from a FASTA sequence database.]]"
           "[[html:<p>This command creates an index of fragment ions from a FASTA sequence database."
           "This indexing engine was developed by Jimmy Eng at the University of Washington "
           "Proteomics Resource.</p>]]";
}

/**
 * \returns the command arguments
 */
vector<string> CometIndexApplication::getArgs() const {
    string arr[] = {"database_name"};
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}


/**
 * \returns the command outputs
 */
vector<pair<string, string>> CometIndexApplication::getOutputs() const {
    vector<pair<string, string>> outputs;
    outputs.push_back(make_pair("comet.params.txt",
                                "a file containing the name and value of all "
                                "parameters/options for the "
                                "current operation. Not all parameters in the "
                                "file may have been used in "
                                "the operation. The resulting file can be used "
                                "with the --parameter-file "
                                "option for other crux programs."));
    outputs.push_back(make_pair(
        "comet.log.txt",
        "a log file containing a copy of all messages that were printed to "
        "standard error."));
    return outputs;
}

COMMAND_T CometIndexApplication::getCommand() const {
    return COMET_INDEX_COMMAND;
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
