/*
 * The original tide-index has been implemented by Benjamin Diament, (I guess). and it has been 
 reimplemented (not form scratch) by Attila Kertesz-Farkas. The sorting on disk has been 
 implemented by Larry Frank Acquaye in March 2022.
 The pipe-line of the new tide-search is the following:
 1. Genertate all the target peptides (with redundancy). The peptides are either stored in 
    the memory or dumped in a text file.
 2. Sort the target peptides
 3. Filter the target peptides and keep the unique peptides, and collect the location 
    of the peptides in different proteins, 
 4. Generate modified target peptides, 
 5. Generate decoy peptides for each modified (and unmodified) peptides, so they are 
    paired and can be printed together nicely.
 6. Note that, in order to keep the set of target and decoy peptides disjunt, one does 
    not need to store all the peptides in a set. It is enough to keep a set of unique peptides
    with the very same neutral mass. This can be done becase the decoy peptide generation 
    does not change the mass of the peptides.
 */

#include <cstdio>
#include <fstream>
#include "io/carp.h"
#include "util/CarpStreamBuf.h"
#include "util/AminoAcidUtil.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include "GeneratePeptides.h"
#include "TideIndexApplication.h"
#include "TideMatchSet.h"
#include "app/tide/modifications.h"
#include "app/tide/records_to_vector-inl.h"
#include "ParamMedicApplication.h"
#include <boost/algorithm/string.hpp>

#include <regex>
#include <assert.h>

#ifdef _MSC_VER
#include <io.h>
#endif
#define CHECK(x) GOOGLE_CHECK(x)


// Larry's code
const char* sortedPeptideFile = "sortedPepTarget.txt";
const char* peptideFile = "pepTarget.txt";
DECLARE_int32(fifo_page_size);
// Larry's code ends here

extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
                                const string& input_filename,
                                const string& output_filename);
extern unsigned long AddMods(HeadedRecordReader* reader,
                    string out_file,
                    string tmpDir,                    
                    const pb::Header& header,
                    const vector<const pb::Protein*>& proteins,
                    VariableModTable* var_mod_table);
DECLARE_int32(max_mods);
DECLARE_int32(min_mods);
DECLARE_int32(modsoutputter_file_threshold);

TideIndexApplication::TideIndexApplication() {
}

TideIndexApplication::~TideIndexApplication() {
}

int TideIndexApplication::main(int argc, char** argv) {
  return main(Params::GetString("protein fasta file"),
              Params::GetString("index name"),
              StringUtils::Join(vector<string>(argv, argv + argc), ' '));
}

int TideIndexApplication::main(
  const string& fasta,
  const string& index,
  string cmd_line
) {
  carp(CARP_INFO, "Running tide-index...");

  if (cmd_line.empty()) {
    cmd_line = "crux tide-index " + fasta + " " + index;
  }

  // Reroute stderr
  CarpStreamBuf buffer;
  streambuf* old = cerr.rdbuf();
  cerr.rdbuf(&buffer);

  // Get options
  bool overwrite = Params::GetBool("overwrite");  
  double min_mass = Params::GetDouble("min-mass");
  double max_mass = Params::GetDouble("max-mass");
  int min_length = Params::GetInt("min-length");
  int max_length = Params::GetInt("max-length");
  bool monoisotopic_precursor = Params::GetString("isotopic-mass") != "average";
  FLAGS_max_mods = Params::GetInt("max-mods");
  FLAGS_min_mods = Params::GetInt("min-mods");
  FLAGS_modsoutputter_file_threshold = Params::GetInt("modsoutputter-threshold");
  bool allowDups = Params::GetBool("allow-dups");
  if (FLAGS_min_mods > FLAGS_max_mods) {
    carp(CARP_FATAL, "The value for 'min-mods' cannot be greater than the value "
                     "for 'max-mods'");
  }
  bool sort_on_disk = (Params::GetString("sort") == string("disk"));
  
  MASS_TYPE_T mass_type = (monoisotopic_precursor) ? MONO : AVERAGE;
  int missed_cleavages = Params::GetInt("missed-cleavages");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  ENZYME_T enzyme_t = get_enzyme_type_parameter("enzyme");
  const char* enzymePtr = enzyme_type_to_string(enzyme_t);
  string enzyme(enzymePtr);
  if ((enzyme != "no-enzyme") && 
      (digestion != FULL_DIGEST && digestion != PARTIAL_DIGEST)) {
    carp(CARP_FATAL, "'digestion' must be 'full-digest' or 'partial-digest'");
  }

  DECOY_TYPE_T decoy_type = get_tide_decoy_type_parameter("decoy-format");
  string decoyPrefix = Params::GetString("decoy-prefix");  

  ofstream* out_target_decoy_list = NULL;  
  if (Params::GetBool("peptide-list")) {
     out_target_decoy_list = create_stream_in_path(make_file_path(
      "tide-index.peptides.txt").c_str(), NULL, overwrite);
  }
  
  ofstream* out_decoy_fasta = GeneratePeptides::canGenerateDecoyProteins() ?
    create_stream_in_path(make_file_path(
      "tide-index.decoy.fasta").c_str(), NULL, overwrite) : NULL;
	  
  string out_proteins = FileUtils::Join(index, "protix");
  string out_peptides = FileUtils::Join(index, "pepix");
  string auxLocsPbFile = FileUtils::Join(index, "auxlocs");
  string modless_peptides = out_peptides + ".nomods.tmp";
  string peakless_peptides = out_peptides + ".nopeaks.tmp";

  if (create_output_directory(index.c_str(), overwrite) != 0) {
    carp(CARP_FATAL, "Error creating index directory");
  } else if (FileUtils::Exists(out_proteins) ||
             FileUtils::Exists(out_peptides) ||
             FileUtils::Exists(auxLocsPbFile)) {
    if (overwrite) {
      carp(CARP_DEBUG, "Cleaning old index file(s)");
      FileUtils::Remove(out_proteins);
      FileUtils::Remove(out_peptides);
      FileUtils::Remove(auxLocsPbFile);
      FileUtils::Remove(modless_peptides);
      FileUtils::Remove(peakless_peptides);
    } else {
      carp(CARP_FATAL, "Index file(s) already exist, use --overwrite T or a "
                       "different index name");
    }
  }
  int numDecoys;
  switch (decoy_type) {
    case NO_DECOYS:
      numDecoys = 0;
      break;
    case PEPTIDE_SHUFFLE_DECOYS:
      numDecoys = Params::GetInt("num-decoys-per-target");
      break;
    default:
      numDecoys = 1;
      break;
  }

  bool shuffle = decoy_type == PEPTIDE_SHUFFLE_DECOYS;  
  
  if (decoy_type != PEPTIDE_SHUFFLE_DECOYS && numDecoys > 1) {
    carp(CARP_FATAL, "Cannot generate multiple decoys per target in non-shuffle decoy-format!");
  }
  
  // Set up output paths
  if (!FileUtils::Exists(fasta)) {
    carp(CARP_FATAL, "Fasta file %s does not exist", fasta.c_str());
  }

 // Start tide-index
  carp(CARP_INFO, "Reading %s and computing unmodified peptides...",
       fasta.c_str());


  VariableModTable var_mod_table;
  var_mod_table.ClearTables();
  //parse regular amino acid modifications
  string mods_spec = Params::GetString("mods-spec");
  carp(CARP_DEBUG, "mods_spec='%s'", mods_spec.c_str());
  if (!var_mod_table.Parse(mods_spec.c_str())) {
    carp(CARP_FATAL, "Error parsing mods");
  }
  //parse terminal modifications
  mods_spec = Params::GetString("cterm-peptide-mods-spec");
  if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), CTPEP)) {
    carp(CARP_FATAL, "Error parsing c-terminal peptide mods");
  }
  mods_spec = Params::GetString("nterm-peptide-mods-spec");
  if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), NTPEP)) {
    carp(CARP_FATAL, "Error parsing n-terminal peptide mods");
  }
  mods_spec = Params::GetString("cterm-protein-mods-spec");
  if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), CTPRO)) {
    carp(CARP_FATAL, "Error parsing c-terminal protein mods");
  }
  mods_spec = Params::GetString("nterm-protein-mods-spec");
  if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), NTPRO)) {
    carp(CARP_FATAL, "Error parsing n-terminal protein mods");
  }
  var_mod_table.SerializeUniqueDeltas();
  if (!MassConstants::Init(var_mod_table.ParsedModTable(), 
    var_mod_table.ParsedNtpepModTable(), 
    var_mod_table.ParsedCtpepModTable(),
    var_mod_table.ParsedNtproModTable(),
    var_mod_table.ParsedCtproModTable(), MassConstants::bin_width_, MassConstants::bin_offset_)) {
    carp(CARP_FATAL, "Error in MassConstants::Init");
  }
  
  // Create protocol buffer for the protein sequences
  pb::Header proteinPbHeader;  
  proteinPbHeader.Clear();
  proteinPbHeader.set_file_type(pb::Header::RAW_PROTEINS);
  proteinPbHeader.set_command_line(cmd_line);
  pb::Header_Source* headerSource = proteinPbHeader.add_source();
  headerSource->set_filename(AbsPath(fasta));
  headerSource->set_filetype("fasta");
  HeadedRecordWriter proteinWriter(out_proteins, proteinPbHeader);


  // Generate peptide sequences via in silico cleavage. This was in the fastaToPb function before
    
  // Container for the protein header and protein seuqnces.
  ProteinVec vProteinHeaderSequence;  
  
  string proteinHeader;
  std::string proteinSequence;

  FixPt minMassFixPt = MassConstants::ToFixPt(min_mass);
  FixPt maxMassFixPt = MassConstants::ToFixPt(max_mass);
   
  ifstream fastaStream(fasta.c_str(), ifstream::in);

  unsigned long invalidPepCnt = 0;
  unsigned long failedDecoyCnt = 0;

  unsigned long targetsGenerated = 0;
  FILE* fp;
  if (sort_on_disk)
    fp = fopen(peptideFile, "w");  // Peptides stored in this file to be sorted on disk.
  
  int curProtein = -1;  
  pb::Header header_with_mods;
  
  vector<TideIndexPeptide> peptide_list;
  
  // Iterate over all proteins in FASTA file and generate target peptides (with redundancy)
  while (GeneratePeptides::getNextProtein(fastaStream, &proteinHeader, &proteinSequence)) {
	
    // Write pb::Protein
    const pb::Protein* pbProtein = writePbProtein(proteinWriter, ++curProtein, proteinHeader, proteinSequence);
	  // Store the pretein header and the protein sequence
	  vProteinHeaderSequence.push_back(pbProtein);
	
    vector<GeneratePeptides::PeptideReference> cleavedPeptides = GeneratePeptides::cleaveProteinTideIndex(
      &proteinSequence, enzyme_t, digestion, missed_cleavages, min_length, max_length);

    // Iterate over all generated peptides for this protein
    for (vector<GeneratePeptides::PeptideReference>::iterator i = cleavedPeptides.begin();
         i != cleavedPeptides.end(); ++i) {
			 
      FixPt pepMass = calcPepMassTide(&(*i), mass_type, proteinSequence);
      if (pepMass == 0) {
        // Sequence contained some invalid character
        carp(CARP_DEBUG, "Ignoring invalid sequence <%s>", std::string(proteinSequence.data()+i->pos_,i->length_).c_str());  
        ++invalidPepCnt;
        continue;
      } else if (pepMass < minMassFixPt || pepMass > maxMassFixPt) {
        // Skip to next peptide if not in mass range
        continue;
      }
      if (sort_on_disk) {
        fprintf(fp, "%u,%s,%d,%d,%d\n", (unsigned int)(pepMass), std::string(proteinSequence.data() + i->pos_, i->length_).c_str(), curProtein, i->pos_, i->length_);
      } else {
        peptide_list.push_back(TideIndexPeptide(pepMass, i->length_, &(pbProtein->residues()), curProtein, i->pos_, -1));
      }
      ++targetsGenerated;

    }
    if ((curProtein+1) % 10000 == 0) {
      carp(CARP_INFO, "Processed %d protein sequences", curProtein+1);
    }
  }
  if (sort_on_disk)
    fclose(fp);

  if (targetsGenerated == 0) {
    carp(CARP_FATAL, "No target sequences generated.  Is \'%s\' a FASTA file?",
         fasta.c_str());
  }
  if (invalidPepCnt > 0) {
    carp(CARP_INFO, "Ignoring %u peptide sequences containing unrecognized characters.", invalidPepCnt);
  }
  carp(CARP_INFO, "Generated %u targets, including duplicates.", targetsGenerated);

  // Sort peptides
  carp(CARP_INFO, "Sorting peptides");
  
  if (sort_on_disk) {
    
    remove(sortedPeptideFile);
    #ifdef _WIN32
      std::cout << "Windows OS\n";
	  std::string import_csv = "Import-Csv -Header Mass, Protein, A, B, C -Delimiter ',' " + std::string(peptideFile);
	  std::string actual_sorting = "Select-Object { [int]$_.Mass }, Protein, A,B,C";
	  std::string write_to_file = "ConvertTo-Csv -NoTypeInformation | Select-Object -Skip 1 | Set-Content -Path " + std::string(sortedPeptideFile);
	  std::string cmd = "Powershell -Command \"" + import_csv + "|" + actual_sorting + "|" + write_to_file + "\"";
	  int systemResult = system(cmd.c_str());
	  if (systemResult != 0)
		  carp(CARP_FATAL, "System sort failed, the call returned code: %i", systemResult);
    #elif __linux__
	  std::cout << "Linux Based OS\n";
      std::string cmd = "sort  -t, -k 1n,1 -k 2,2 -o " +  std::string(sortedPeptideFile) + " " + std::string(peptideFile);
      int systemResult = system(cmd.c_str()); 
      if (systemResult != 0) 
        carp(CARP_FATAL, "System sort failed, the call returned code: %i", systemResult);
    #elif __unix__
      std::cout << "Other UNIX OS\n";
	  std::string cmd = "sort  -t, -k 1n,1 -k 2,2 -o " + std::string(sortedPeptideFile) + " " + std::string(peptideFile);
	  int systemResult = system(cmd.c_str());
	  if (systemResult != 0)
		  carp(CARP_FATAL, "System sort failed, the call returned code: %i", systemResult);
    #elif __APPLE__
      std::cout << "Apple OS\n";
	  std::string cmd = "sort  -t, -k 1n,1 -k 2,2 -o " + std::string(sortedPeptideFile) + " " + std::string(peptideFile);
	  int systemResult = system(cmd.c_str());
	  if (systemResult != 0)
		  carp(CARP_FATAL, "System sort failed, the call returned code: %i", systemResult);
    #else
      std::cout << "OS not supported";
	  carp(CARP_FATAL, "OS not supported");
      std::string cmd = ""
    #endif
    // TODO Add code to check if sorted file exists, and if the size of the sorted file is the same as the pepTarget file 
  } else {
    sort(peptide_list.begin(), peptide_list.end(), less<TideIndexPeptide>());    
  }
 
  // Prepare the protocol buffer for the peptides.  
  carp(CARP_INFO, "Writting peptides");

  // pb::Header header_with_mods;
  pb::Header_PeptidesHeader& pep_header = *(header_with_mods.mutable_peptides_header());
  
  pep_header.Clear();
  pep_header.set_min_mass(min_mass);
  pep_header.set_max_mass(max_mass);
  pep_header.set_min_length(min_length);
  pep_header.set_max_length(max_length);
  pep_header.set_monoisotopic_precursor(monoisotopic_precursor);
  pep_header.set_enzyme(enzyme);
  if (enzyme != "no-enzyme") {
    pep_header.set_full_digestion(digestion == FULL_DIGEST);
    pep_header.set_max_missed_cleavages(missed_cleavages);
  }
  pep_header.mutable_mods()->CopyFrom(*(var_mod_table.ParsedModTable()));
  pep_header.mutable_nterm_mods()->CopyFrom(*(var_mod_table.ParsedNtpepModTable()));
  pep_header.mutable_cterm_mods()->CopyFrom(*(var_mod_table.ParsedCtpepModTable()));
  pep_header.mutable_nprotterm_mods()->CopyFrom(*(var_mod_table.ParsedNtproModTable()));
  pep_header.mutable_cprotterm_mods()->CopyFrom(*(var_mod_table.ParsedCtproModTable()));

  pep_header.set_decoys_per_target(numDecoys);

  header_with_mods.set_file_type(pb::Header::PEPTIDES);
  header_with_mods.set_command_line(cmd_line);
  pb::Header_Source* source = header_with_mods.add_source();
  source->mutable_header()->CopyFrom(proteinPbHeader);
  source->set_filename(AbsPath(out_proteins));

  pb::Header header_no_mods;
  header_no_mods.CopyFrom(header_with_mods);
  pb::ModTable* del = header_no_mods.mutable_peptides_header()->mutable_mods();
  del->mutable_variable_mod()->Clear();
  del->mutable_unique_deltas()->Clear();

  bool need_mods = var_mod_table.Unique_delta_size() > 0;

  string peptidePbFile = need_mods ? modless_peptides : peakless_peptides;  
  
  // Check header
  if (header_no_mods.source_size() != 1) {
    carp(CARP_FATAL, "header_no_mods had a number of sources other than 1");
  }
  
  headerSource = header_no_mods.mutable_source(0);
  if (!headerSource->has_filename() || headerSource->has_filetype()) {
    carp(CARP_FATAL, "pbHeader source invalid");
  }

  // Now check other desired settings
  if (!header_no_mods.has_peptides_header()) {
    carp(CARP_FATAL, "!header_no_mods->has_peptideHeapheader()");
  }
  const pb::Header_PeptidesHeader& settings = header_no_mods.peptides_header();
  //if (!Peptide::SetMinMaxMassAndLength(settings)) {
  //  carp(CARP_FATAL, "Error setting min/max mass/length");
  if (!settings.has_enzyme() || settings.enzyme().empty()) {
    carp(CARP_FATAL, "Enzyme settings error");
  }

  header_no_mods.set_file_type(pb::Header::PEPTIDES);
  header_no_mods.mutable_peptides_header()->set_has_peaks(false);
  header_no_mods.mutable_peptides_header()->set_decoys(decoy_type);

  pb::Peptide pbPeptide;
  unsigned long count = 0;
  unsigned long numTargets = 0;
  unsigned long numDuplicateTargets = 0;
  unsigned long peptide_cnt = 0;
  
  if (!sort_on_disk && peptide_list.size() == 0)
    carp(CARP_FATAL, "No peptides were generated.");

  ifstream sortedFile(sortedPeptideFile);
  unsigned long numLines = 0;
  TideIndexPeptide* currentPeptide;
  TideIndexPeptide* duplicatedPeptide;
  // Filter peptides and keep the unique target peptides and gather the 
  // location of the peptide in other protein sequences 
  if (sort_on_disk) {
    currentPeptide = getNextPeptide(sortedFile, vProteinHeaderSequence);  // get the first peptide  
  } else {
    currentPeptide = &(peptide_list[peptide_cnt++]);  // get the first peptide  
  }
  if (1==1) {  // This is needed because we need to destroy the peptideWriter and pbAuxLoc later. Ugly solution :/
    // Create the auxiliary locations header and writer
    pb::Header auxLocsHeader;
    auxLocsHeader.set_file_type(pb::Header::AUX_LOCATIONS);
    pb::Header_Source* auxLocsSource = auxLocsHeader.add_source();
    auxLocsSource->set_filename(peptidePbFile);
    auxLocsSource->mutable_header()->CopyFrom(header_no_mods);
    HeadedRecordWriter auxLocWriter(auxLocsPbFile, auxLocsHeader);  
    pb::AuxLocation pbAuxLoc;
    int auxLocIdx = -1;
    
    HeadedRecordWriter peptideWriter(peptidePbFile, header_no_mods); // put header in outfile	  
    while (currentPeptide != nullptr) {

      while (true) {

        if (sort_on_disk) {
          duplicatedPeptide = getNextPeptide(sortedFile, vProteinHeaderSequence);  // get a peptide  
          if (duplicatedPeptide == nullptr) {
            break;
          }
        } else {
          if (peptide_cnt >= peptide_list.size()){
            duplicatedPeptide = NULL;
            break;
          }
          duplicatedPeptide = &(peptide_list[peptide_cnt++]);  // get a peptide  
        }
        if( (*duplicatedPeptide) == (*currentPeptide)) {

          numDuplicateTargets++;
          carp(CARP_DEBUG, "Skipping duplicate %s.", currentPeptide->getSequence().c_str());
          pb::Location* location = pbAuxLoc.add_location();
          location->set_protein_id(duplicatedPeptide->getProteinId());
          location->set_pos(duplicatedPeptide->getProteinPos());
          if (sort_on_disk){ 
            delete duplicatedPeptide;
          }
        } else {
          break;
        }
        if (duplicatedPeptide->getMass() < currentPeptide->getMass()){  // Check if sorting worked properly.
          carp(CARP_INFO, "peptide mass: %lf, subsequent peptide mass %lf", currentPeptide->getMass(), duplicatedPeptide->getMass());
          carp(CARP_FATAL, "Peptides are not sorted correctly. Sorting seems to be failed. Try again and check the free disk space.");
        }
      }

      getPbPeptide(count, *currentPeptide, pbPeptide);
      // Not all peptides have aux locations associated with them. Check to see
      // if GetGroup added any locations to aux_location. If yes, only then
      // assign the corresponding array index to the peptide and write it out.
      if (pbAuxLoc.location_size() > 0) {
        pbPeptide.set_aux_locations_index(++auxLocIdx);
        auxLocWriter.Write(&pbAuxLoc);
        pbAuxLoc.Clear();
      }
      // Write the peptide AFTER the aux_locations check, in case we added an
      // aux_locations_index to the peptide.
      peptideWriter.Write(&pbPeptide);

      numTargets++;
      if (++count % 1000000 == 0) {
        carp(CARP_INFO, "Wrote %u unique target peptides", count);
      }
      numLines++;
      if (sort_on_disk){      
        delete currentPeptide;
      }
      currentPeptide = duplicatedPeptide;
    }
  }
  carp(CARP_DETAILED_INFO, "%u peptides in file", numLines);
  
  // Release the memory allocated.
  peptide_list.clear();
  vector<TideIndexPeptide> tmp;
  peptide_list.swap(tmp);
  // peptide_list = vector<TideIndexPeptide>();
  
  carp(CARP_INFO, "Skipped %u duplicate targets.",
       numDuplicateTargets);
  
  carp(CARP_INFO, "Generated %u unique target peptides.", numTargets);

  peptidePbFile = modless_peptides;  
  if (need_mods) {
    carp(CARP_INFO, "Computing modified peptides...");
    HeadedRecordReader reader(modless_peptides, NULL, 1024 << 10); // 1024kb buffer
    numTargets = AddMods(&reader, peakless_peptides, Params::GetString("temp-dir"), header_with_mods, vProteinHeaderSequence, &var_mod_table);
    carp(CARP_INFO, "Created %u modified and unmodified target peptides.", numTargets);
    peptidePbFile = peakless_peptides;
  }
  
  if (numDecoys > 0) {
      carp(CARP_INFO, "Generating %d decoy(s) per target peptide", numDecoys);
  } else {
      carp(CARP_INFO, "No decoy peptides will be generated");
  }
  //Reader for the peptides:
  pb::Header aaf_peptides_header;
  HeadedRecordReader aaf_peptide_reader(peptidePbFile, &aaf_peptides_header);

  if (aaf_peptides_header.file_type() != pb::Header::PEPTIDES ||
    !aaf_peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptidePbFile.c_str());
  }

  FifoAllocator fifo_alloc_peptides_(FLAGS_fifo_page_size << 20);
  RecordReader* reader_;
  reader_ = aaf_peptide_reader.Reader();
  pb::Peptide current_pb_peptide_;

  MassConstants::Init(&aaf_peptides_header.peptides_header().mods(),
                    &aaf_peptides_header.peptides_header().nterm_mods(),
                    &aaf_peptides_header.peptides_header().cterm_mods(),
                    &aaf_peptides_header.peptides_header().nprotterm_mods(),
                    &aaf_peptides_header.peptides_header().cprotterm_mods(),
                    MassConstants::bin_width_, MassConstants::bin_offset_);

  bool success;
  vector<int> decoy_peptide_idx;
  int startLoc;
  int protein_id;
  int mod_code;
  int decoy_index;
  int mod_index;
  int unique_delta;
  double delta;
  unsigned long decoy_count = 0;
  double mass;
  pb::Protein* decoy_pd_protein;
  int generateAttemptsMax = 6;

  string target_peptide_with_mods;
  string decoy_peptide_with_mods;
  int prot_id, pos, len;

  
  pb::Header new_header;
  new_header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
  subheader->CopyFrom(aaf_peptides_header.peptides_header());
  subheader->set_has_peaks(true);
  source = new_header.add_source();
  source->mutable_header()->CopyFrom(aaf_peptides_header);
  source->set_filename(AbsPath(peptidePbFile));
  HeadedRecordWriter writer(out_peptides, new_header);
  
  // Read peptides protocol buffer file
  vector<const pb::AuxLocation*> locations;
  if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxLocsPbFile)) {
    carp(CARP_FATAL, "Error reading auxlocs file");
  }
  int mass_precision = Params::GetInt("mass-precision");
  int mod_precision  = Params::GetInt("mod-precision");
  
  
  CHECK(reader_->OK());
  CHECK(writer.OK());
  string decoy_peptide_str;
  
  vector<pb::Peptide> pb_peptides; // Used
  set<string> peptide_str_set;
  double last_mass = -1.0;
  pb::Peptide last_pb_peptide;
  const pb::Protein* protein;
  string pepmass_str;
  string pos_str;
  string mod_str;
  int mod_pos_offset;
  
  /* The trick to keep the sets target and decoy peptides disjunt is that:
  One does not need to keep all the unique target peptides in the memory 
  and check every time whether a decoy peptide already exists as a target. 
  It is enought to keep the target in a set (in the memory) peptdes having 
  exactly the same mass. It is because the decoy generation does not chage 
  the mass of the peptide.
  */ 
  if (out_target_decoy_list) {
    *out_target_decoy_list  << "target\t";
    if ( numDecoys > 0)
      *out_target_decoy_list  << "decoy(s)\t";
    *out_target_decoy_list  << "mass\tproteins" << std::endl;
  }
  // Go over the (modified and unmodified) peptides from the protocol buffer and generate decoy peptides 
  bool done = false;
  
  while (!done) {
    while (!done) { // Gather peptides with the same mass if allowDups is false
      done = reader_->Done();
      if (done == true)
        break;
      reader_->Read(&last_pb_peptide);
      if (allowDups) {
        pb_peptides.push_back(last_pb_peptide);
        break;
      }
      if (pb_peptides.empty()) {
        pb_peptides.push_back(last_pb_peptide);
        last_mass = last_pb_peptide.mass();
        continue;
      }
      if (last_mass < last_pb_peptide.mass()) {   // The mass has increased, new set of peptides 
        break;
      }
      pb_peptides.push_back(last_pb_peptide);		
    }  
    
    // Create a set with the unique peptides sequences. The peptides must have the same neutral mass.
    for (vector<pb::Peptide>::iterator pb_pept_itr = pb_peptides.begin(); pb_pept_itr != pb_peptides.end(); ++pb_pept_itr) {
      target_peptide_with_mods = getModifiedPeptideSeq(&(*pb_pept_itr), &vProteinHeaderSequence);
      peptide_str_set.insert(target_peptide_with_mods);
    }
    // For each peptide in the set: 
    // 1. store the target peptide in the index file, 
    // 2. generate decoy peptides and store them in the index file too.
    for (vector<pb::Peptide>::iterator pb_pept_itr=pb_peptides.begin(); pb_pept_itr != pb_peptides.end(); ++pb_pept_itr) {
      // Get the peptide and write it 
      current_pb_peptide_ = (*pb_pept_itr);
      CHECK(writer.Write(&current_pb_peptide_));

      // Get the peptide sequence with modifications
      if (out_target_decoy_list) {
        target_peptide_with_mods = getModifiedPeptideSeq(&(*pb_pept_itr), &vProteinHeaderSequence);
        *out_target_decoy_list  << target_peptide_with_mods;
      }
      
      protein_id = current_pb_peptide_.first_location().protein_id();      
      startLoc = current_pb_peptide_.first_location().pos();
      
      if (numDecoys > 0) {  // Get peptide sequence without mods
        bool first_decoy = true;
        if (out_target_decoy_list) {
          *out_target_decoy_list << '\t';
        }      
        string target_peptide = vProteinHeaderSequence[protein_id]->residues().substr(startLoc, current_pb_peptide_.length());
        string decoy_peptide_str_with_mods;

        //	Generate a decoy peptides:
        protein = vProteinHeaderSequence[protein_id];
        for (int i = 0; i < numDecoys; ++i) {
          
          shuffle = decoy_type == PEPTIDE_SHUFFLE_DECOYS;  

          for (int j = 0; j < generateAttemptsMax; ++j) {
            // Generates a permutation for how generate the decoy peptide from target peptide
            GeneratePeptides::makeDecoyIdx(target_peptide, shuffle, decoy_peptide_idx);
            decoy_peptide_str = target_peptide;
            
            // Create the decoy peptide sequence
            for(int k = 0; k < decoy_peptide_idx.size(); ++k) {
              decoy_peptide_str[decoy_peptide_idx[k]] = target_peptide[k];
            }	
            decoy_peptide_str_with_mods = decoy_peptide_str;
            // Add the modificaitons to the decoy:
            if (current_pb_peptide_.modifications_size() > 0){
              mod_pos_offset = 0;
              vector<double> deltas(decoy_peptide_str.length());
              for (int m = 0; m < current_pb_peptide_.modifications_size(); ++m) {
                mod_code = current_pb_peptide_.modifications(m);
                MassConstants::DecodeMod(mod_code, &mod_index, &delta);
                decoy_index =  decoy_peptide_idx[mod_index];
                deltas[decoy_index] = delta;
              }
              for (int d = 0; d < deltas.size(); ++d){
                if (deltas[d] == 0.0) continue;              
                mod_str = '[' + StringUtils::ToString(deltas[d], mod_precision) + ']';
                decoy_peptide_str_with_mods.insert(d + 1 + mod_pos_offset, mod_str);
                mod_pos_offset += mod_str.length();
                
              }
            }            
            // Check if this modified decoy peptide has not been generated yet.
            if (allowDups) {
              success = true;
              break;
            } else {
              success = peptide_str_set.find(decoy_peptide_str_with_mods) == peptide_str_set.end();
            }
            if (success == true){
              peptide_str_set.insert(decoy_peptide_str_with_mods);
              break;
            }
            shuffle = true; // Failed to generate decoy, so try shuffling in the next attempt.
          }
          if (success == false) {
            carp(CARP_DEBUG, "Failed to generate decoys for sequence %s", target_peptide.c_str());
            ++failedDecoyCnt;
            continue; // it could be a 'break;' too
          }

          // According to the indeces create a decoy protein string,
          decoy_pd_protein = writeDecoyPbProtein(++curProtein, protein, decoy_peptide_str, startLoc, proteinWriter);

          // Create a protocol buffer peptide object for the decoy peptide. Note that the decoy peptide may contain modifications.
          pb::Peptide decoy_current_pb_peptide_ = current_pb_peptide_;
          if (current_pb_peptide_.modifications_size() > 0) {
            decoy_current_pb_peptide_.clear_modifications();
            for (int m = 0; m < current_pb_peptide_.modifications_size(); ++m) {
              mod_code = current_pb_peptide_.modifications(m);

              MassConstants::DecodeMod(mod_code, &mod_index, &unique_delta);
              decoy_index = decoy_peptide_idx[mod_index];
              mod_code = MassConstants::EncodeMod(decoy_index, unique_delta);
              decoy_current_pb_peptide_.add_modifications(mod_code);
            }
          }
          decoy_current_pb_peptide_.set_id(numTargets + decoy_count++);
          decoy_current_pb_peptide_.clear_first_location();
          decoy_current_pb_peptide_.mutable_first_location()->set_protein_id(curProtein);
          decoy_current_pb_peptide_.mutable_first_location()->set_pos((startLoc > 0) ? 1 : 0);
          decoy_current_pb_peptide_.set_decoy_index(i);
          CHECK(writer.Write(&decoy_current_pb_peptide_));
          delete decoy_pd_protein;
          
          //report the decoy peptide if needed.
          if (out_target_decoy_list) {
            if (first_decoy == false)
              *out_target_decoy_list  <<',';
            *out_target_decoy_list << decoy_peptide_str_with_mods.c_str();
            first_decoy = false;
          }
          if (decoy_count % 1000000 == 0) {
            carp(CARP_INFO, "Wrote %u decoy peptides", decoy_count);
          }          
        }
      }
      // Print 1) the peptide neutral mass, 2) protein header of origin and 3) the locations of the target peptides
      if (out_target_decoy_list){
        string pepmass_str = StringUtils::ToString(current_pb_peptide_.mass(), mass_precision);
        *out_target_decoy_list  <<'\t' << pepmass_str;

        pos_str =  StringUtils::ToString(startLoc + 1, 1);
        string proteinNames = vProteinHeaderSequence[protein_id]->name() + '(' + pos_str + ')';
        if (current_pb_peptide_.has_aux_locations_index()) {
          const pb::AuxLocation* aux = locations[current_pb_peptide_.aux_locations_index()];
          for (int i = 0; i < aux->location_size(); ++i) {
            const pb::Location& location = aux->location(i);
            protein = vProteinHeaderSequence[location.protein_id()];
            pos_str =  StringUtils::ToString(location.pos() + 1, 1);
            proteinNames += ',' + protein->name() + '(' + pos_str + ')';
          }
        }
        *out_target_decoy_list  <<'\t' << proteinNames << endl;
      }
      
    }
    pb_peptides.clear();
    peptide_str_set.clear();
    if (!allowDups){
      pb_peptides.push_back(last_pb_peptide);		
      last_mass = last_pb_peptide.mass();
    }	  
  }

  if (out_target_decoy_list){			
    out_target_decoy_list->close();
    delete out_target_decoy_list;
  }  
    if (failedDecoyCnt > 0) {
    carp(CARP_INFO, "Failed to generate decoys for %u low complexity peptides.", failedDecoyCnt);
  }
  carp(CARP_INFO, "Generated %u target peptides.", numTargets);
  carp(CARP_INFO, "Generated %u decoy peptides.", decoy_count);
  carp(CARP_INFO, "Generated %u peptides in total.", numTargets+decoy_count);
    
  // Clean up
/*  for (vector<const pb::Protein*>::iterator i = vProteinHeaderSequence.begin();
       i != vProteinHeaderSequence.end();
       ++i) {
    delete *i;
  }
  */// Recover stderr
  cerr.rdbuf(old);
  // The destructor is explicitly called in order to release the resource to enable deleting. 
  // This was added to resolve the race condition issue which arises on windows.
  aaf_peptide_reader.~HeadedRecordReader(); 
  FileUtils::Remove(modless_peptides);
  FileUtils::Remove(peakless_peptides);
  return 0;

}

string TideIndexApplication::getName() const {
  return "tide-index";
}

string TideIndexApplication::getDescription() const {
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

vector<string> TideIndexApplication::getArgs() const {
  string arr[] = {
    "protein fasta file",
    "index name"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> TideIndexApplication::getOptions() const {
  string arr[] = {
    "allow-dups",
    "clip-nterm-methionine",
    "cterm-peptide-mods-spec",
    "cterm-protein-mods-spec",
    "custom-enzyme",
    "decoy-format",
    "decoy-prefix",
    "digestion",
    "enzyme",
    "isotopic-mass",
    "keep-terminal-aminos",
    "mass-precision",
    "max-length",
    "max-mass",
    "max-mods",
    "min-length",
    "min-mass",
    "min-mods",
    "missed-cleavages",
    "mod-precision",
    "mods-spec",
    "nterm-peptide-mods-spec",
    "nterm-protein-mods-spec",
    "auto-modifications",
    "auto-modifications-spectra",
    "num-decoys-per-target",
    "output-dir",
    "overwrite",
    "parameter-file",
    "peptide-list",
    "seed",
    "sort",
    "temp-dir",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > TideIndexApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("index",
    "A binary index, using the name specified on the command line."));
  outputs.push_back(make_pair("tide-index.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("tide-index.log.txt",
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution."));
  return outputs;
}

bool TideIndexApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T TideIndexApplication::getCommand() const {
  return TIDE_INDEX_COMMAND;
}

FixPt TideIndexApplication::calcPepMassTide(
  GeneratePeptides::PeptideReference* pep,
  MASS_TYPE_T massType,
  string prot
) {
  FixPt mass;
  FixPt aaMass;
  const MassConstants::FixPtTableSet *_tables;

  if (massType == AVERAGE) {
    mass = MassConstants::fixp_avg_h2o;
    _tables = &MassConstants::avg_tables;
  } else if (massType == MONO) {
    mass = MassConstants::fixp_mono_h2o;
    _tables = &MassConstants::mono_tables;
  } else {
    carp(CARP_FATAL, "Invalid mass type");
  }

  for (size_t i = 0; i < pep->length_; ++i) {
    if (i == 0) {
      if(pep->pos_ == 0)  //apply protein terminal mod if this is protein N-terminal
        aaMass = _tables->nprotterm_table[prot.at(0)];
      else //apply peptide N-terminal mod 
        aaMass = _tables->nterm_table[prot.at(pep->pos_)];
    } else if (i == pep->length_ - 1) {
      if((pep->pos_ + pep->length_) == prot.length())  //check if this is protein C-terminal
        aaMass = _tables->cprotterm_table[prot.at(pep->pos_ + i)];
      else
        aaMass = _tables->cterm_table[prot.at(pep->pos_ + i)];
    } else {
      aaMass = _tables->_table[prot.at(pep->pos_ + i)];
    }
    if (aaMass == 0) {
      return 0;
    }
    mass += aaMass;
  }
  return mass;
}

pb::Protein* TideIndexApplication::writePbProtein(
  HeadedRecordWriter& writer,
  int id,
  const string& name,
  const string& residues,
  int targetPos
) {
  pb::Protein* p = new pb::Protein;
  p->Clear();
  p->set_id(id);
  p->set_name(name);
  p->set_residues(residues);
  if (targetPos >= 0) {
    p->set_target_pos(targetPos);
  }
  writer.Write(p);
  return p;
}

/*
 * This is a bit tricky. We are storing decoy peptide sequences as
 * "pseudo-proteins" in the protocol buffer.  To make this work, we
 * have to store some additional bits of information: the identity of
 * the preceding and following amino acids, as well as the identify of
 * the corresponding target sequence.  All of this information gets
 * appended together before getting put into the protocol buffer.
 * Note that the two termini are handled differently: if there is no
 * preceding amino acid, then nothing is prepended; but if there is no
 * succeeding amino acid, then a hyphen is appended.
 */
pb::Protein* TideIndexApplication::writeDecoyPbProtein(
  int id,
  const pb::Protein* protein,
  string decoyPeptideSequence,
  int startLoc,
  HeadedRecordWriter& proteinWriter
) {
  const string proteinSequence = protein->residues();
  const int pepLen = decoyPeptideSequence.length();

  // Add N term to decoySequence, if it exists
  
  if (startLoc > 0) {
    decoyPeptideSequence.insert(0, 1, proteinSequence.at(startLoc - 1));
  }
  
  // Add C term to decoySequence, if it exists, or hyphen otherwise.
  size_t cTermLoc = startLoc + pepLen;
  
  decoyPeptideSequence.push_back((cTermLoc < proteinSequence.length()) ?
    proteinSequence.at(cTermLoc) : '-');
  
  // Append original target sequence
  decoyPeptideSequence.append(proteinSequence.substr(startLoc, pepLen));
  
  return writePbProtein(proteinWriter, id, Params::GetString("decoy-prefix") + protein->name(),
                 decoyPeptideSequence, startLoc);
}

void TideIndexApplication::getPbPeptide(
  int id,
  const TideIndexPeptide& peptide,
  pb::Peptide& outPbPeptide
) {
  outPbPeptide.Clear();
  outPbPeptide.set_id(id);
  outPbPeptide.set_mass(peptide.getMass());
  outPbPeptide.set_length(peptide.getLength());
  outPbPeptide.mutable_first_location()->set_protein_id(peptide.getProteinId());
  outPbPeptide.mutable_first_location()->set_pos(peptide.getProteinPos());
  if (peptide.isDecoy()) {
    outPbPeptide.set_decoy_index(peptide.decoyIdx());
  }
}

void TideIndexApplication::addAuxLoc(
  int proteinId,
  int proteinPos,
  pb::AuxLocation& outAuxLoc
) {
  pb::Location* location = outAuxLoc.add_location();
  location->set_protein_id(proteinId);
  location->set_pos(proteinPos);
}

void TideIndexApplication::processParams() {
  if (Params::GetBool("auto-modifications")) {
    if (!Params::IsDefault("mods-spec")) {
      carp(CARP_FATAL, "Automatic modification inference cannot be used with user specified "
                       "modifications. Please rerun with either auto-modifications set to 'false' "
                       "or with modifications turned off.");
    }
    vector<string> files = StringUtils::Split(Params::GetString("auto-modifications-spectra"), ',');
    for (vector<string>::iterator i = files.begin(); i != files.end(); ) {
      if ((*i = StringUtils::Trim(*i)).empty()) {
        i = files.erase(i);
      } else {
        i++;
      }
    }
    if (files.empty()) {
      carp(CARP_FATAL, "Spectra files must be specified with the 'auto-modifications-spectra' "
                       "parameter when 'auto-modifications' is enabled.");
    }
    vector<ParamMedic::RunAttributeResult> modsResult;
    ParamMedicApplication::processFiles(files, false, true, NULL, &modsResult);
    vector<ParamMedic::Modification> mods = ParamMedic::Modification::GetFromResults(modsResult);
    vector<string> modStrings;
    vector<string> modNStrings;
    vector<string> modCStrings;
    for (vector<ParamMedic::Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
      string location = i->getLocation();
      const double mass = i->getMassDiff();
      const bool variable = i->getVariable();

      vector<string>* modStringVector;
      string modCountStr = variable ? "4" : "";

      if (location == ParamMedic::Modification::LOCATION_NTERM) {
        modStringVector = &modNStrings;
        location = "X";
      } else if (location == ParamMedic::Modification::LOCATION_CTERM) {
        modStringVector = &modCStrings;
        location = "X";
      } else {
        modStringVector = &modStrings;
      }
      modStringVector->push_back(modCountStr + location + (mass >= 0 ? '+' : '-') +
        StringUtils::ToString(mass));
    }
    Params::Set("mods-spec", StringUtils::Join(modStrings, ','));
    Params::Set("nterm-peptide-mods-spec", StringUtils::Join(modNStrings, ','));
    Params::Set("cterm-peptide-mods-spec", StringUtils::Join(modCStrings, ','));
  }

  // Update mods-spec parameter for default cysteine mod
  string default_cysteine = "C+" + StringUtils::ToString(CYSTEINE_DEFAULT);
  string mods_spec = Params::GetString("mods-spec");
  if (mods_spec.find('C') == string::npos) {
    mods_spec = mods_spec.empty() ?
      default_cysteine : default_cysteine + ',' + mods_spec;
    carp(CARP_DETAILED_INFO, "Using default cysteine mod '%s' ('%s')",
         default_cysteine.c_str(), mods_spec.c_str());
  }
  Params::Set("mods-spec", mods_spec);

  // Override enzyme if it is something other than "custom-enzyme"
  // when a custom enzyme is specified
  if (!Params::GetString("custom-enzyme").empty() &&
      Params::GetString("enzyme") != "custom-enzyme") {
    Params::Set("enzyme", "custom-enzyme");
    carp(CARP_WARNING, "'custom-enzyme' was set: setting 'enzyme' to 'custom-enzyme'");
  }
}
// Why is this here? It is not a TideIndexApplication member function. -AKF
string getModifiedPeptideSeq(const pb::Peptide* peptide,
  const ProteinVec* proteins) {
  int mod_index;
  double mod_delta;
  stringstream mod_stream;
  const pb::Location& location = peptide->first_location();
  const pb::Protein* protein = proteins->at(location.protein_id());
  // Get peptide sequence without mods
  string pep_str = protein->residues().substr(location.pos(), peptide->length());

  // Store all mod indices/deltas
  map<int, double> mod_map;
  set<int> mod_indices;
  for (int j = 0; j < peptide->modifications_size(); ++j) {
    //        var_mod_table.DecodeMod(ModCoder::Mod(peptide->modifications(j)), &mod_index, &mod_delta);
    MassConstants::DecodeMod(ModCoder::Mod(peptide->modifications(j)),
      &mod_index, &mod_delta);
    mod_indices.insert(mod_index);
    mod_map[mod_index] = mod_delta;
  }
  int modPrecision = Params::GetInt("mod-precision");
  for (set<int>::const_reverse_iterator j = mod_indices.rbegin();
    j != mod_indices.rend();
    ++j) {
    // Insert the modification string into the peptide sequence
    mod_stream << '[' << StringUtils::ToString(mod_map[*j], modPrecision) << ']';
    pep_str.insert(*j + 1, mod_stream.str());
    mod_stream.str("");
  }
  return pep_str;
}

// Larry's code
TideIndexApplication::TideIndexPeptide* TideIndexApplication::getNextPeptide(ifstream &sortedFile, ProteinVec& vProteinHeaderSequence){
  
  string line;
  vector<std::string> strs;
  if (getline(sortedFile, line)) {

	#ifdef _WIN32
	  // This extra code is added to remove redundant quotes in string.
	  boost::split(strs, line, boost::is_any_of(","));
	  for (int i = 0; i < strs.size(); i++) {
		  strs[i].erase(remove(strs[i].begin(), strs[i].end(), '\"'), strs[i].end());
	  }

	#else
	  boost::split(strs, line, boost::is_any_of(","));
	#endif
    
    FixPt mass = stoul(strs[0]);
	// string peptide sequence is skipped stoi(strs[1]);
    int proteinId = stoi(strs[2]);
    int proteinPos = stoi(strs[3]);
    int length = stoi(strs[4]);
    int decoyIdx = -1;//stoi(strs[6]); // -1 if not a decoy; There are no decoy peptides generated at this point

    const string& proteinSequence = vProteinHeaderSequence[proteinId]->residues();
    
    TideIndexPeptide* pepTarget = new TideIndexPeptide(mass, length, &proteinSequence, proteinId, proteinPos, decoyIdx);
	
    return pepTarget;
  }else{
    return nullptr;
  }
  
}
// Larry's code ends here
/*
* Local Variables:
* mode: c
* c-basic-offset: 2
* End:
*/
