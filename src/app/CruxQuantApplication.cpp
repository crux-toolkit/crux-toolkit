#include "io/carp.h"
#include "CruxQuantApplication.h"

using std::make_pair;

CruxQuantApplication::CruxQuantApplication() {}

CruxQuantApplication::~CruxQuantApplication() {}

int CruxQuantApplication::main(int argc, char** argv) {
    carp(CARP_INFO, "Running crux-quant...");
    return 0;
}

string CruxQuantApplication::getName() const {
    return "crux-quant";
}

// TODO: Add better description
string CruxQuantApplication::getDescription() const {
    return "Quantify peptides/proteins from MS/MS data";
}

vector<string> CruxQuantApplication::getArgs() const {
    string arr[] = {"spectrum-files", "output-dir"};
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> CruxQuantApplication::getOptions() const {
    string arr[] = {"overwrite", "parameter-file", "help"};
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > CruxQuantApplication::getOutputs() const {
    vector< pair<string, string> > outputs;
    outputs.push_back(make_pair("crux-quant-output-dir", "Directory containing the output files"));
    return outputs;
}

COMMAND_T CruxQuantApplication::getCommand() const {
    return CRUX_QUANT_COMMAND;
}

bool CruxQuantApplication::needsOutputDirectory() const {
    return true;
}

// TODO: Add parameter processing
void CruxQuantApplication::processParams() {

}

