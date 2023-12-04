#include "CruxQuantApplication.h"

#include <cmath>
#include <exception>
#include <sstream>

#include "IndexedMassSpectralPeak.h"
#include "crux-quant/Results.h"
#include "crux-quant/Utils.h"
#include "io/carp.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/crux-utils.h"

using std::make_pair;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

using pwiz::cv::MS_ms_level;
using pwiz::cv::MS_scan_start_time;
using pwiz::msdata::BinaryDataArrayPtr;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::SpectrumListSimple;
using pwiz::msdata::SpectrumListSimplePtr;
using pwiz::msdata::SpectrumPtr;

typedef pwiz::msdata::SpectrumListPtr SpectrumListPtr;

CruxQuantApplication::CruxQuantApplication() {}

CruxQuantApplication::~CruxQuantApplication() {}

int CruxQuantApplication::main(int argc, char** argv) {
    string psm_file = Params::GetString("lfq-peptide-spectrum matches");
    vector<string> spec_files = Params::GetStrings("spectrum files");
    return main(psm_file, spec_files);
}

int CruxQuantApplication::main(const string& psm_file, const vector<string>& spec_files) {
    carp(CARP_INFO, "Running crux-lfq...");

    if (!FileUtils::Exists(psm_file)) {
        carp(CARP_FATAL, "PSM file %s not found", psm_file.c_str());
    }
    map<int, CruxQuant::PSM> psm_datum = CruxQuant::create_psm_map(psm_file);
    CruxQuant::CruxLFQResults lfqResults(spec_files);

    for (const string& spectra_file : spec_files) {
        SpectrumListPtr spectra_ms1 = loadSpectra(spectra_file, 1);
        SpectrumListPtr spectra_ms2 = loadSpectra(spectra_file, 2);

        carp(CARP_INFO, "Read %d spectra. for MS1", spectra_ms1->size());
        carp(CARP_INFO, "Read %d spectra. for MS2", spectra_ms2->size());

        CruxQuant::IndexedSpectralResults indexResults = indexedMassSpectralPeaks(spectra_ms1, spectra_file);

        vector<CruxQuant::Identification> allIdentifications = createIdentifications(psm_datum, spectra_file, spectra_ms2);
        unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxQuant::calculateTheoreticalIsotopeDistributions(allIdentifications);

        CruxQuant::setPeakFindingMass(allIdentifications, modifiedSequenceToIsotopicDistribution);
        vector<double> chargeStates = CruxQuant::createChargeStates(allIdentifications);

        CruxQuant::quantifyMs2IdentifiedPeptides(
            spectra_file,
            allIdentifications,
            chargeStates,
            indexResults._ms1Scans,
            indexResults._indexedPeaks,
            modifiedSequenceToIsotopicDistribution,
            lfqResults
        );
       
        CruxQuant::runErrorChecking(spectra_file, lfqResults);
        // For now this happens in the forloop, but it should be moved out of the forloop based on FlashLFQ Code
        if(CruxQuant::QUANTIFY_AMBIGUOUS_PEPTIDES){
            lfqResults.setPeptideModifiedSequencesAndProteinGroups(allIdentifications);
        }
        
    }
    lfqResults.calculatePeptideResults(CruxQuant::QUANTIFY_AMBIGUOUS_PEPTIDES);
    lfqResults.calculateProteinResultsMedianPolish(CruxQuant::USE_SHARED_PEPTIDES_FOR_PROTEIN_QUANT);
    const std::string results_file = make_file_path("crux-lfq.txt");
    lfqResults.writeResults(results_file);

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

SpectrumListPtr CruxQuantApplication::loadSpectra(const string& file, int msLevel) {
    try {
        MSDataFile msd(file);
        SpectrumListPtr originalSpectrumList = msd.run.spectrumListPtr;
        if (!originalSpectrumList) {
            carp(CARP_FATAL, "Error reading spectrum file %s", file.c_str());
        }
        SpectrumListSimplePtr filteredSpectrumList(new SpectrumListSimple);

        for (size_t i = 0; i < originalSpectrumList->size(); ++i) {
            SpectrumPtr spectrum = originalSpectrumList->spectrum(i);

            int spectrumMSLevel = spectrum->cvParam(MS_ms_level).valueAs<int>();
            if (spectrumMSLevel == msLevel) {
                // Add the spectrum to the filtered list
                filteredSpectrumList->spectra.push_back(spectrum);
            }
        }

        return filteredSpectrumList;
    } catch (const std::exception& e) {
        carp(CARP_INFO, "Error:  %s", e.what());
        return nullptr;
    }
}

IndexedSpectralResults CruxQuantApplication::indexedMassSpectralPeaks(SpectrumListPtr spectrum_collection, const string& spectra_file) {
    string _spectra_file(spectra_file);

    map<int, map<int, IndexedMassSpectralPeak>> _indexedPeaks;
    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;

    _ms1Scans[_spectra_file] = vector<Ms1ScanInfo>();

    IndexedSpectralResults index_results{_indexedPeaks, _ms1Scans};

    if (!spectrum_collection) {
        return index_results;
    }

    int _scanIndex = 0;
    int _oneBasedScanNumber = 1;

    for (size_t i = 0; i < spectrum_collection->size(); ++i) {
        SpectrumPtr spectrum = spectrum_collection->spectrum(i);

        if (spectrum) {
            BinaryDataArrayPtr mzs = spectrum->getMZArray();
            BinaryDataArrayPtr intensities = spectrum->getIntensityArray();

            int scanIndex;
            int oneBasedScanNumber;
            std::string scanId = getScanID(spectrum->id);
            if (scanId.empty()) {
                scanIndex = _scanIndex;
                oneBasedScanNumber = _oneBasedScanNumber;
            } else {
                scanIndex = std::stoi(scanId);
                oneBasedScanNumber = scanIndex + 1;
            }

            double retentionTime = spectrum->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();

            if (mzs && intensities) {
                const std::vector<double>& mzArray = mzs->data;
                const std::vector<double>& intensityArray = intensities->data;
                for (size_t j = 0; j < mzArray.size(); ++j) {
                    FLOAT_T mz = mzArray[j];
                    int roundedMz = static_cast<int>(std::round(mz * BINS_PER_DALTON));
                    IndexedMassSpectralPeak spec_data(
                        mz,                 // mz value
                        intensityArray[j],  // intensity
                        scanIndex,          // zeroBasedMs1ScanIndex
                        retentionTime);

                    auto& indexedPeaks = index_results._indexedPeaks;
                    auto it = indexedPeaks.find(roundedMz);
                    if (it == indexedPeaks.end()) {
                        map<int, IndexedMassSpectralPeak> tmp;
                        tmp.insert({scanIndex, spec_data});
                        indexedPeaks[roundedMz] = tmp;
                    } else {
                        it->second.insert({scanIndex, spec_data});
                    }

                    Ms1ScanInfo scan = {oneBasedScanNumber, scanIndex, retentionTime};
                    index_results._ms1Scans[spectra_file].push_back(scan);
                }
            }

            _scanIndex++;
            _oneBasedScanNumber++;
        }
    }
    return index_results;
}

// Make this a multithreaded process
vector<Identification> CruxQuantApplication::createIdentifications(const map<int, PSM>& psm_datum, const string& spectra_file, SpectrumListPtr spectrum_collection) {
    carp(CARP_INFO, "Creating indentifications, this may take a bit of time, do not terminate the process...");

    vector<Identification> allIdentifications;
    string _spectra_file(spectra_file);

    for (size_t i = 0; i < spectrum_collection->size(); ++i) {
        SpectrumPtr spectrum = spectrum_collection->spectrum(i);
        if (spectrum) {
            int scanIndex;
            std::string scanId = getScanID(spectrum->id);
            if (scanId.empty()) {
                continue;
            } else {
                scanIndex = std::stoi(scanId);
            }

            auto it = psm_datum.find(scanIndex);

            if (it != psm_datum.end()) {
                double retentionTimeInSeconds = spectrum->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();

                FLOAT_T retentionTimeInMinutes = retentionTimeInSeconds / 60.0;

                Identification identification;

                identification.sequence = it->second.sequence_col;
                identification.monoIsotopicMass = it->second.peptide_mass_col;
                identification.charge = it->second.charge_col;
                identification.peptideMass = it->second.peptide_mass_col;
                identification.precursorCharge = it->second.spectrum_precursor_mz_col;
                identification.spectralFile = _spectra_file;
                identification.ms2RetentionTimeInMinutes = retentionTimeInMinutes;
                identification.scanId = it->second.scan_col;
                identification.modifications = it->second.modifications;
                allIdentifications.push_back(identification);
            }
        }
    }

    return allIdentifications;
}