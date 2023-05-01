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

// TODO: Change the description - this description is copied from tideindex.
string CruxQuantApplication::getDescription() const {
    return
        "[[nohtml:Create an index for all peptides in a fasta file, for use in "
        "subsequent calls to tide-search.]]"
        "[[html:<p>Tide is a tool for identifying peptides from tandem mass "
        "spectra. It is an independent reimplementation of the SEQUEST<sup>&reg;"
        "</sup> algorithm, which assigns peptides to spectra by comparing the "
        "observed spectra to a catalog of theoretical spectra derived from a "
        "database of known proteins. Tide's primary advantage is its speed. Our "
        "published paper provides more detail on how Tide works. If you use Tide "
        "in your research, please cite:</p><blockquote>Benjamin J. Diament and "
        "William Stafford Noble. &quot;<a href=\""
        "http://dx.doi.org/10.1021/pr101196n\">Faster SEQUEST Searching for "
        "Peptide Identification from Tandem Mass Spectra.</a>&quot; <em>Journal of "
        "Proteome Research</em>. 10(9):3871-9, 2011.</blockquote><p>The <code>"
        "tide-index</code> command performs an optional pre-processing step on the "
        "protein database, converting it to a binary format suitable for input to "
        "the <code>tide-search</code> command.</p><p>Tide considers only the "
        "standard set of 21 amino acids. Peptides containing non-amino acid "
        "alphanumeric characters (BJXZ) are skipped. Non-alphanumeric characters "
        "are ignored completely.</p>]]";
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

