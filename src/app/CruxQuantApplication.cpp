#include "io/carp.h"
#include "io/SpectrumRecordWriter.h"
#include "io/PSMReader.h"
#include "io/MatchFileReader.h"
#include "io/MzIdentMLReader.h"
#include "io/SQTReader.h"
#include "io/PepXMLReader.h"
#include "util/Params.h"
#include "util/crux-utils.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include "TideSearchApplication.h"
#include "CruxQuantApplication.h"

using std::make_pair;

CruxQuantApplication::CruxQuantApplication() {}

CruxQuantApplication::~CruxQuantApplication() {}

int CruxQuantApplication::main(int argc, char** argv) {
    string psm_file = Params::GetString("lfq-peptide-spectrum matches");
    vector<string> spec_files = Params::GetStrings("spectrum files");
    return main(psm_file, spec_files);
}

int CruxQuantApplication::main(const string& psm_file, const vector<string>& spec_files){
    carp(CARP_INFO, "Running crux-lfq...");
    // MatchCollection* psm = getPSM(psm_file);
    vector<InputFile> ms1_spectra_files = getSpecFiles(spec_files, 1);
    vector<InputFile> ms2_spectra_files = getSpecFiles(spec_files, 2);
    
    return 0;
}

string CruxQuantApplication::getName() const {
    return "crux-lfq";
}

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
        "spectrum files"
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

// I may need to rewrite this code - loadSpectra can be found in TideSearchApplication.
// Currently all this code does is generates spec records from spec files.
//  Must revisit this code to have a better understanding of what it does.
vector<InputFile> CruxQuantApplication::getSpecFiles(const vector<string>& filepaths, int ms_level) const {
  vector<InputFile> input_sr;

  if (Params::GetString("spectrum-parser") != "pwiz") { 
    carp(CARP_FATAL, "spectrum-parser must be pwiz instead of %s", Params::GetString("spectrum-parser").c_str() ); 
  }

  for (vector<string>::const_iterator f = filepaths.begin(); f != filepaths.end(); f++) {
    string spectrum_input_url = *f;
    string spectrumrecords_url = make_file_path(FileUtils::BaseName(spectrum_input_url) + ".spectrumrecords.ms" + to_string(ms_level));
    carp(CARP_INFO, "Converting %s to spectrumrecords %s", spectrum_input_url.c_str(), spectrumrecords_url.c_str());
    carp(CARP_DEBUG, "New MS%d spectrumrecords filename: %s", ms_level, spectrumrecords_url.c_str());

    if (!FileUtils::Exists(spectrumrecords_url)) {
      if (!SpectrumRecordWriter::convert(spectrum_input_url, spectrumrecords_url, ms_level, true)) {
        carp(CARP_FATAL, "Error converting MS2 spectrumrecords from %s", spectrumrecords_url.c_str());
      }
    }
    input_sr.push_back(InputFile(*f, spectrumrecords_url, true));
  }

  return input_sr;
}

// MatchCollection* CruxQuantApplication::getPSM(string psm_file){
//     Database* data; // Figure out the import of this data and why it's needed
//     PSMReader* reader;
//     if (StringUtils::IEndsWith(psm_file, ".txt")) {
//       reader = new MatchFileReader(psm_file.c_str(), data);
//     } else if (StringUtils::IEndsWith(psm_file, ".html")) {
//       carp(CARP_FATAL, "HTML format has not been implemented yet");
//     } else if (StringUtils::IEndsWith(psm_file, ".sqt")) {
//       reader = new SQTReader(psm_file.c_str(), data);
//     } else if (StringUtils::IEndsWith(psm_file, ".pin")) {
//       carp(CARP_FATAL, "Pin format has not been implemented yet");
//     } else if (StringUtils::IEndsWith(psm_file, ".xml")) {
//       reader = new PepXMLReader(psm_file.c_str(), data);
//     } else if (StringUtils::IEndsWith(psm_file, ".mzid")) {
//       reader = new MzIdentMLReader(psm_file.c_str(), data);
//     } else {
//       carp(CARP_FATAL, "Could not determine input format, "
//            "Please name your files ending with .txt, .html, .sqt, .pin, "
//            ".xml, .mzid or use the --input-format option to "
//            "specify file type");
//     }
//     MatchCollection* collection = reader->parse();

//     return collection;
// }