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
    vector<string> input_files = Params::GetStrings("spectrum files");
    return main(psm_file, input_files);
}

int CruxQuantApplication::main(const string& psm_file, const vector<string>& input_files){
    carp(CARP_INFO, "Running crux-lfq...");
    std::map<std::string, SpectrumCollection*> spectra_;
    
    vector<InputFile> sr = getInputFiles(input_files);
    
    for (vector<InputFile>::const_iterator f = sr.begin(); f != sr.end(); f++) {
        string spectra_file = f->SpectrumRecords;
        SpectrumCollection* spectra = NULL;
        map<string, SpectrumCollection*>::iterator spectraIter = spectra_.find(spectra_file);
        if (spectraIter == spectra_.end()) {
            carp(CARP_INFO, "Reading spectrum file %s.", spectra_file.c_str());
            spectra = loadSpectra(spectra_file);
            carp(CARP_INFO, "Read %d spectra.", spectra->Size());
        } else {
            spectra = spectraIter->second;
        }
    }
    
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
vector<InputFile> CruxQuantApplication::getInputFiles(const vector<string>& filepaths) const {
   // Try to read all spectrum files as spectrumrecords, convert those that fail
  vector<InputFile> input_sr;
  for (vector<string>::const_iterator f = filepaths.begin(); f != filepaths.end(); f++) {
    SpectrumCollection spectra;
    pb::Header spectrum_header;
    string spectrumrecords = *f;
    bool keepSpectrumrecords = true;
    if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
      // Failed, try converting to spectrumrecords file
      carp(CARP_INFO, "Converting %s to spectrumrecords format", f->c_str());
      carp(CARP_INFO, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);

      
      spectrumrecords = make_file_path(FileUtils::BaseName(*f) + ".spectrumrecords.tmp");
     
      carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());
      if (!SpectrumRecordWriter::convert(*f, spectrumrecords)) {
        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", f->c_str());
      }
      carp(CARP_DEBUG, "Reading converted spectrum file %s", spectrumrecords.c_str());
      // Re-read converted file as spectrumrecords file
      if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
        carp(CARP_DEBUG, "Deleting %s", spectrumrecords.c_str());
        FileUtils::Remove(spectrumrecords);
        carp(CARP_FATAL, "Error reading spectra file %s", spectrumrecords.c_str());
      }
    }
    input_sr.push_back(InputFile(*f, spectrumrecords, keepSpectrumrecords));
  }
  return input_sr;
}

SpectrumCollection* CruxQuantApplication::loadSpectra(const string& file) {
  SpectrumCollection* spectra = new SpectrumCollection();
  pb::Header header;
  if (!spectra->ReadSpectrumRecords(file, &header)) {
    carp(CARP_FATAL, "Error reading spectrum file %s", file.c_str());
  }
  return spectra;
}