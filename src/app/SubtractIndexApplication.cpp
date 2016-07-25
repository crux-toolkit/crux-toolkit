/**
 * \file SubtractIndexApplication.cpp
 * \brief Iterative PSM meta-search via Cascade protocol
 ************************************************************/
#include "SubtractIndexApplication.h"
#include "io/OutputFiles.h"
#include "AssignConfidenceApplication.h"
#include "TideSearchApplication.h"
#include "TideIndexApplication.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "util/FileUtils.h"
#include "io/carp.h"
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#define CHECK(x) GOOGLE_CHECK(x)

using namespace std;

/**
 * \returns a blank SubtractIndexApplication object
 */
SubtractIndexApplication::SubtractIndexApplication() {
}

/**
 * Destructor
 */
SubtractIndexApplication::~SubtractIndexApplication() {
}

/**
 * main method for SubtractIndexApplication
 */
int SubtractIndexApplication::main(int argc, char** argv) {
  carp(CARP_INFO, "Running subtract-index...");

  bool overwrite = Params::GetBool("overwrite");
  bool has_decoys = false;

  //open tide index 1
  const string index1 = Params::GetString("tide index 1");  
  bool write_peptides = FileUtils::Exists(index1 + "/tide-index.peptides.target.txt");
  string peptides_file1 = index1 + "/pepix";
  string proteins_file1 = index1 + "/protix";
  string auxlocs_file1 = index1 + "/auxlocs";

  carp(CARP_INFO, "Reading index %s", index1.c_str());
  ProteinVec proteins1;
  pb::Header protein_header1;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins1,
    proteins_file1, &protein_header1)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file1.c_str());
  }
  carp(CARP_DEBUG, "Read %d proteins", proteins1.size());
  
  pb::Header peptides_header1;
  HeadedRecordReader peptide_reader1(peptides_file1, &peptides_header1);
  if (peptides_header1.file_type() != pb::Header::PEPTIDES ||
    !peptides_header1.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file1.c_str());
  }
  const pb::Header::PeptidesHeader& pepHeader1 = peptides_header1.peptides_header();
  DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader1.decoys();
  if (headerDecoyType != NO_DECOYS) {
    has_decoys = true;
    if (headerDecoyType == PROTEIN_REVERSE_DECOYS) {
      TideSearchApplication::PROTEIN_LEVEL_DECOYS = true;
    }
  }
  MassConstants::Init(&peptides_header1.peptides_header().mods(), 
    &peptides_header1.peptides_header().nterm_mods(),
    &peptides_header1.peptides_header().cterm_mods(), 0.0, 0.0);

  //open tide index 2
  const string index2 = Params::GetString("tide index 2");
  string peptides_file2 = index2 + "/pepix";
  string proteins_file2 = index2 + "/protix";
  carp(CARP_INFO, "Reading index %s", index2.c_str());
  pb::Header peptides_header2;
  HeadedRecordReader peptide_reader2(peptides_file2, &peptides_header2);
  ProteinVec proteins2;
  pb::Header protein_header2;
  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins2,
    proteins_file2, &protein_header2)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file2.c_str());
  }
  carp(CARP_DEBUG, "Read %d proteins", proteins2.size());

  //output files;
  const string index_out = Params::GetString("output index");
  string out_proteins = index_out + "/" + "protix";
  string out_peptides = index_out + "/" + "pepix";
  string out_aux = index_out + "/" + "auxlocs";

  if (create_output_directory(index_out.c_str(), overwrite) != 0) {
    carp(CARP_FATAL, "Error creating index directory");
  } else if (FileUtils::Exists(out_proteins) ||
             FileUtils::Exists(out_peptides) ||
             FileUtils::Exists(out_aux)) {
    if (overwrite) {
      carp(CARP_DEBUG, "Cleaning old index file(s)");
      remove(out_proteins.c_str());
      remove(out_peptides.c_str());
      remove(out_aux.c_str());
    } else {
      carp(CARP_FATAL, "Index file(s) already exist, use --overwrite T or a "
        "different index name");
    }
  }
  ofstream* out_target_list = NULL;
  ofstream* out_decoy_list = NULL;
  if (write_peptides) {
    out_target_list = create_stream_in_path(make_file_path(
      "tide-index.peptides.target.txt").c_str(), NULL, overwrite);
    if (has_decoys) {
      out_decoy_list = create_stream_in_path(make_file_path(
        "tide-index.peptides.decoy.txt").c_str(), NULL, overwrite);
    }
  }
  //copy aux and protein files;
  FileUtils::Copy(auxlocs_file1, out_aux);
  FileUtils::Copy(proteins_file1, out_proteins);

  pb::Header new_header;
  CHECK(peptides_header1.file_type() == pb::Header::PEPTIDES);
  CHECK(peptides_header1.has_peptides_header());
  new_header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
  subheader->CopyFrom(peptides_header1.peptides_header());
  subheader->set_has_peaks(true);
  pb::Header_Source* source = new_header.add_source();
  source->mutable_header()->CopyFrom(peptides_header1);
  HeadedRecordWriter writer(out_peptides, new_header);
  CHECK(peptide_reader1.OK());
  CHECK(peptide_reader2.OK());
  CHECK(writer.OK());

  pb::Peptide pb_peptide1, pb_peptide2;
  double pep_mass1, pre_pep_mass1, pep_mass2, pre_pep_mass2;
  string pep_str1, pep_str2;
    
  vector<pb::Peptide> peptide_container;
  pep_mass2 = pre_pep_mass2 = 0.0;
  pre_pep_mass1 = 0.0;
  bool found;

  while (!peptide_reader1.Done()) {
    peptide_reader1.Read(&pb_peptide1);
    pep_mass1 = pb_peptide1.mass();
    found = false;

    if (pre_pep_mass1 < pep_mass1) {
      peptide_container.clear();
      pre_pep_mass1 = pep_mass1;
    }

    while (pep_mass2 < pep_mass1) {
      if (!peptide_reader2.Done()) {
        peptide_reader2.Read(&pb_peptide2);
        pep_mass2 = pb_peptide2.mass();
      } else { 
        break; 
      }
    }
    while (pep_mass2 == pep_mass1) {
      peptide_container.push_back(pb_peptide2);
      if (!peptide_reader2.Done()) {
        peptide_reader2.Read(&pb_peptide2);
        pep_mass2 = pb_peptide2.mass();
      } else {
        break;
      }
    }
    pep_str1 = getModifiedPeptideSeq(&pb_peptide1, &proteins1);
    for (vector<pb::Peptide>::iterator it = peptide_container.begin(); it != peptide_container.end(); ++it) {
      pep_str2 = getModifiedPeptideSeq(&(*it), &proteins2);
      if (pep_str1 == pep_str2) {
        found = true;
        break;
      }
    }
    if (!found) {
      CHECK(writer.Write(&pb_peptide1));
      carp(CARP_DEBUG, "%s\t%lf\n", pep_str1.c_str(), pep_mass1);

      if (write_peptides) {
        if (!pb_peptide1.is_decoy()) {
          *out_target_list << pep_str1 << '\t' << pb_peptide1.mass() << endl;
        } else {
          *out_decoy_list << pep_str1 << '\t' << pb_peptide1.mass() << endl;
        }
      }
    }
  }

  return 0;
}

/**
 * \returns the command name for SubtractIndexApplication
 */
string SubtractIndexApplication::getName() const {
  return "subtract-index";
}

/**
 * \returns the description for SubtractIndexApplication
 */
string SubtractIndexApplication::getDescription() const {
  return "[[html:<p>This command takes two peptide indices, created by the tide-index "
    "command, and subtracts the second index from the first. The result is an output "
    "index that contains peptides that appear in the first index but not the second.</p>]]"
    "[[nohtml:This command takes two peptide indices, created by the tide-index command, "
    "and subtracts the second index from the first. The result is an output index "
    "that contains peptides that appear in the first index but not the second.]]";
}

/**
 * \returns the command arguments
 */
vector<string> SubtractIndexApplication::getArgs() const {
  string arr[] = {
    "tide index 1",
    "tide index 2",
    "output index"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}


/**
 * \returns the command options
 */
vector<string> SubtractIndexApplication::getOptions() const {
  string arr[] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "peptide-list",
    "output-dir"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > SubtractIndexApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("subtract-index.target.txt",
    "a <a href=\"../file-formats/txt-format.html\">tab-delimited text file</a> containing the "
    "target peptides."));
  outputs.push_back(make_pair("subtract-index.decoy.txt",
    "a <a href=\"../file-formats/txt-format.html\">tab-delimited text file</a> containing the "
    "decoy peptides."));
  outputs.push_back(make_pair("subtract-index.log.txt",
    "a log file containing a copy of all messages that were printed to stderr."));
  outputs.push_back(make_pair("subtract-index.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  return outputs;
}

/**
 * \returns the filestem for SubtractIndexApplication
 */
string SubtractIndexApplication::getFileStem() const {
  return "subtract-index";
}

COMMAND_T SubtractIndexApplication::getCommand() const {
  return CASCADE_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not.
 */
bool SubtractIndexApplication::needsOutputDirectory() const {
  return true;
}

void SubtractIndexApplication::processParams() {

}
