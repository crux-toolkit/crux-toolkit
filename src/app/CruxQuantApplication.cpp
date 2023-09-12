#include <cmath>
#include "io/carp.h"
#include "util/Params.h"
#include "util/crux-utils.h"
#include "util/FileUtils.h"
#include "CruxQuantApplication.h"
#include "crux-quant/Utils.h"
#include "crux-quant/CMercury8.h"

using std::list;
using std::make_pair;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

CruxQuantApplication::CruxQuantApplication() {}

CruxQuantApplication::~CruxQuantApplication() {}

int CruxQuantApplication::main(int argc, char **argv)
{
    string psm_file = Params::GetString("lfq-peptide-spectrum matches");
    vector<string> input_files = Params::GetStrings("spectrum files");
    return main(psm_file, input_files);
}

int CruxQuantApplication::main(const string &psm_file, const vector<string> &input_files)
{
    carp(CARP_INFO, "Running crux-lfq...");

    if (!FileUtils::Exists(psm_file))
    {
        carp(CARP_FATAL, "PSM file %s not found", psm_file.c_str());
    }

    MatchFileReader *matchFileReader = new MatchFileReader(psm_file);

    for (const string &spectra_file : input_files)
    {
        Crux::SpectrumCollection *spectra_ms1 = loadSpectra(spectra_file, 1);

        unordered_map<int, list<CruxQuant::IndexedMassSpectralPeak>> indexes = indexedMassSpectralPeaks(spectra_ms1);
        unordered_map<string, list<pair<double, double>>> _modifiedSequenceToIsotopicDistribution = calculateTheoreticalIsotopeDistributions(matchFileReader);
        double minChargeState = getMinChargeState(matchFileReader);
        double maxChargeState = getMaxChargeState(matchFileReader);

        vector<double> chargeStates = CreateChargeStates(minChargeState, maxChargeState);

        // Print the charge states
        for (const auto& state : chargeStates) {
            std::cout << state << " ";
        }
    }
    delete matchFileReader;
    return 0;
}

string CruxQuantApplication::getName() const
{
    return "crux-lfq";
}

string CruxQuantApplication::getDescription() const
{
    return "[[nohtml:This command reads a set of PSMs and a corresponding set of spectrum files"
           "and carries out label-free quantification (LFQ) for each detected peptide.]]"
           "[[html:<p>This command reads a set of PSMs and a corresponding set of spectrum files "
           "and carries out label-free quantification (LFQ) for each detected peptide."
           "The algorithm follows that of FlashLFQ: "
           "Millikin RJ, Solntsev SK, Shortreed MR, Smith LM. &quot;<a href=\""
           "https://pubmed.ncbi.nlm.nih.gov/29083185/\">Ultrafast Peptide Label-Free Quantification with FlashLFQ.</a>&quot;"
           "<em>Journal of Proteome Research</em>. 17(1):386-391, 2018.</blockquote><p>]]";
}

vector<string> CruxQuantApplication::getArgs() const
{
    string arr[] = {
        "lfq-peptide-spectrum matches",
        "spectrum files"};
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> CruxQuantApplication::getOptions() const
{
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

vector<pair<string, string>> CruxQuantApplication::getOutputs() const
{
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

COMMAND_T CruxQuantApplication::getCommand() const
{
    return CRUX_QUANT_COMMAND;
}

bool CruxQuantApplication::needsOutputDirectory() const
{
    return true;
}

// TODO: Add parameter processing
void CruxQuantApplication::processParams()
{
}

Crux::SpectrumCollection *CruxQuantApplication::loadSpectra(const string &file, int ms_level)
{
    Crux::SpectrumCollection *spectra = SpectrumCollectionFactory::create(file);
    spectra->parse(ms_level = ms_level);
    return spectra;
}

std::unordered_map<int, std::list<CruxQuant::IndexedMassSpectralPeak>> CruxQuantApplication::indexedMassSpectralPeaks(Crux::SpectrumCollection *spectrum_collection)
{

    carp(CARP_INFO, "Read %d spectra. for MS1", spectrum_collection->getNumSpectra());

    // Define the hashmap with a list of IndexedMassSpectralPeak instances as the value type
    std::unordered_map<int, std::list<CruxQuant::IndexedMassSpectralPeak>> _indexedPeaks;

    if (spectrum_collection->getNumSpectra() <= 0)
    {
        return _indexedPeaks;
    }

    int scanIndex = 0;

    for (auto spectrum = spectrum_collection->begin(); spectrum != spectrum_collection->end(); ++spectrum)
    {
        if (*spectrum != nullptr)
        {
            for (auto peak = (*spectrum)->begin(); peak != (*spectrum)->end(); ++peak)
            {
                if (*peak != nullptr)
                {
                    FLOAT_T mz = (*peak)->getLocation();
                    int roundedMz = static_cast<int>(std::round(mz * CruxQuantApplication::BinsPerDalton));
                    CruxQuant::IndexedMassSpectralPeak spec_data(
                        mz,                             // mz value
                        (*peak)->getIntensity(),        // intensity
                        scanIndex,                      // zeroBasedMs1ScanIndex
                        (*spectrum)->getRetentionTime() // retentionTime
                    );

                    // Find the corresponding entry in _indexedPeaks and create a new entry if it doesn't exist
                    auto it = _indexedPeaks.find(roundedMz);
                    if (it == _indexedPeaks.end())
                    {
                        _indexedPeaks.insert({roundedMz, std::list<CruxQuant::IndexedMassSpectralPeak>()});
                    }
                    else
                    {
                        // Add a new IndexedMassSpectralPeak object to the list at the corresponding roundedMz entry
                        _indexedPeaks[roundedMz].emplace_back(spec_data);
                    }
                }
            }
            scanIndex++;
        }
    }

    return _indexedPeaks;
}

unordered_map<string, list<pair<double, double>>> CruxQuantApplication::calculateTheoreticalIsotopeDistributions(MatchFileReader *matchFileReader)
{

    unordered_map<string, list<pair<double, double>>> _modifiedSequenceToIsotopicDistribution;

    // Iterate through the data
    while (matchFileReader->hasNext())
    {
        // Move to the next row
        matchFileReader->next();

        // Access data from columns using the appropriate methods
        string peptide_sequence = matchFileReader->getString(SEQUENCE_COL);
        int charge = matchFileReader->getInteger(CHARGE_COL);
        

        if (_modifiedSequenceToIsotopicDistribution.find(peptide_sequence) != _modifiedSequenceToIsotopicDistribution.end())
        {
            continue;
        }

        std::vector<std::pair<double, double>> isotopicMassesAndNormalizedAbundances;

        string formula = CruxQuant::calcFormula(peptide_sequence);
        carp(CARP_INFO, "formula: %s", formula.c_str());
        char *char_array = new char[formula.length() + 1];
        char *fn = nullptr;
        strcpy(char_array, formula.c_str());
        CMercury8 dist(fn);
        dist.Echo(true);
        dist.GoMercury(char_array, charge);
        vector<double> masses;
        vector<double> abundances;
        for (auto i : dist.FixedData){
            masses.push_back(i.mass);
            abundances.push_back(i.data)
        }
        delete[] char_array;
    }

    return _modifiedSequenceToIsotopicDistribution;
}

double CruxQuantApplication::getMinChargeState(MatchFileReader *matchFileReader)
{

    auto minChargeState = matchFileReader->getDouble(SPECTRUM_PRECURSOR_MZ_COL); // Get the first element

    while (matchFileReader->hasNext())
    {
        matchFileReader->next();
        auto currentChargeState = matchFileReader->getDouble(SPECTRUM_PRECURSOR_MZ_COL);

        if (currentChargeState < minChargeState)
        {
            minChargeState = currentChargeState;
        }
    }

    return minChargeState;
}

double CruxQuantApplication::getMaxChargeState(MatchFileReader *matchFileReader)
{

   
    auto maxChargeState = matchFileReader->getDouble(SPECTRUM_PRECURSOR_MZ_COL);// Get the first element

    while (matchFileReader->hasNext())
    {
        matchFileReader->next();
        auto currentChargeState = matchFileReader->getDouble(SPECTRUM_PRECURSOR_MZ_COL);

        if (currentChargeState > maxChargeState)
        {
            maxChargeState = currentChargeState;
        }
    }

    return maxChargeState;
}

vector<double> CruxQuantApplication::CreateChargeStates(double minChargeState, double maxChargeState){
    vector<double> chargeStates;

    for (double i = minChargeState; i <= maxChargeState; i++) {
        chargeStates.push_back(i);
    }

    return chargeStates;

}