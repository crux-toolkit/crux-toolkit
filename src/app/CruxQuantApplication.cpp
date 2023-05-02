#include "io/carp.h"
#include "CruxQuantApplication.h"

using std::make_pair;

CruxQuantApplication::CruxQuantApplication() {}

CruxQuantApplication::~CruxQuantApplication() {}

int CruxQuantApplication::main(int argc, char** argv) {
    carp(CARP_INFO, "Running crux-lfq...");
    return 0;
}

string CruxQuantApplication::getName() const {
    return "crux-lfq";
}

// TODO: Change the description - this description is copied from tideindex.
string CruxQuantApplication::getDescription() const {
    return
        "[[nohtml:This command reads a set of PSMs and a corresponding set of spectrum files"
        "and carries out label-free quantification (LFQ) for each detected peptide.]]"
        "[[html:<p>This command reads a set of PSMs and a corresponding set of spectrum files "
        "and carries out label-free quantification (LFQ) for each detected peptide."
        "The algorithm follows that of FlashLFQ: "
        "Millikin RJ, Solntsev SK, Shortreed MR, Smith LM. &quot;<a href=\""
        "https://pubmed.ncbi.nlm.nih.gov/29083185/\">Ultrafast Peptide Label-Free Quantification with FlashLFQ.</a>&quot;" 
        "<em>Journal of Proteome Research</em>. 17(1):386-391, 2018.</blockquote><p>]]";
}

vector<string> CruxQuantApplication::getArgs() const {
    string arr[] = {
        "peptide-spectrum matches", 
        "spectrum file"
    };
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> CruxQuantApplication::getOptions() const {
    string arr[] = {
        "score", 
        "threshold", 
        "smaller-is-better",
        "fileroot",
        "output-dir",
        "overwrite",
        "parameter-file",
        "spectrum-parser",
        "verbosity"
    };
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > CruxQuantApplication::getOutputs() const {
    vector< pair<string, string> > outputs;
    outputs.push_back(make_pair("crux-lfq.txt", 
        "A tab-delimited text file in which rows are peptides, "
        "columns correspond to the different spectrum files, "
        "and values are peptide quantifications.  "
        "If a peptide is not detected in a given run, "
        "then its corresponding quantification value is NaN."));
    outputs.push_back(make_pair("crux-lfq.params.txt",
        "A file containing the name and value of all parameters/options"
        " for the current operation. Not all parameters in the file may have"
        " been used in the operation. The resulting file can be used with the "
        "--parameter-file option for other Crux programs."));
    outputs.push_back(make_pair("crux-lfq.log.txt", 
        "A log file containing a copy of all messages that were printed to the screen during execution."
    ));
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

