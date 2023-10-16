#include "CruxQuantApplication.h"

#include <cmath>
#include "IndexedMassSpectralPeak.h"
#include "crux-quant/Utils.h"
#include "io/carp.h"
#include "util/FileUtils.h"
#include "util/Params.h"

using std::make_pair;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

CruxQuantApplication::CruxQuantApplication() {}

CruxQuantApplication::~CruxQuantApplication() {}

int CruxQuantApplication::main(int argc, char **argv) {
    string psm_file = Params::GetString("lfq-peptide-spectrum matches");
    vector<string> spec_files = Params::GetStrings("spectrum files");
    return main(psm_file, spec_files);
}

int CruxQuantApplication::main(const string &psm_file, const vector<string> &spec_files) {
    carp(CARP_INFO, "Running crux-lfq...");

    if (!FileUtils::Exists(psm_file)) {
        carp(CARP_FATAL, "PSM file %s not found", psm_file.c_str());
    }

    MatchFileReader *matchFileReader = new MatchFileReader(psm_file);

    for (const string &spectra_file : spec_files) {
        SpectrumCollection* spectra_ms1 = CruxQuant::loadSpectra(spectra_file, true);
        SpectrumCollection* spectra_ms2 = CruxQuant::loadSpectra(spectra_file, false);

        spectra_ms1->

        carp(CARP_INFO, "Read %d spectra. for MS1", spectra_ms1->SpectrumNum);
        carp(CARP_INFO, "Read %d spectra. for MS2", spectra_ms2->getNumSpectra());

        // CruxQuant::IndexedSpectralResults indexResults = CruxQuant::indexedMassSpectralPeaks(spectra_ms1, spectra_file);
        
        // // Replace spectra_file in createIdentifications with file column in the matchFileReader
        // vector<CruxQuant::Identification> allIdentifications = CruxQuant::createIdentifications(matchFileReader, spectra_file, spectra_ms2);
        // unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxQuant::calculateTheoreticalIsotopeDistributions(allIdentifications);
        
        // CruxQuant::setPeakFindingMass(allIdentifications, modifiedSequenceToIsotopicDistribution);
        // vector<double> chargeStates = CruxQuant::createChargeStates(allIdentifications);

        // CruxQuant::quantifyMs2IdentifiedPeptides(spectra_file, allIdentifications, chargeStates, indexResults._ms1Scans, indexResults._indexedPeaks);
      


       
    }
    delete matchFileReader;
    return 0;
}

string CruxQuantApplication::getName() const {
    return "crux-lfq";
}

string CruxQuantApplication::getDescription() const {
    return "[[nohtml:This command reads a set of PSMs and a corresponding set of spectrum files"
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
        "lfq-peptide-spectrum matches",
        "spectrum files"};
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
        "verbosity"};
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<pair<string, string>> CruxQuantApplication::getOutputs() const {
    vector<pair<string, string>> outputs;
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
                                "A log file containing a copy of all messages that were printed to the screen during execution."));
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
