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
#include "app/tide/modifications.h"
#include "app/tide/records_to_vector-inl.h"
#include "ParamMedicApplication.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "residue_stats.pb.h"
#include "crux_version.h"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <regex>
#include <assert.h>
#include <filesystem>

#ifdef _MSC_VER
#include <io.h>
#endif
#define CHECK(x) GOOGLE_CHECK(x)

std::string peptideFile = "pepTarget.txt";
string TideIndexApplication::tide_index_mzTab_filename_ = "tide-index.params.mztab";

extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
                                const string& input_filename,
                                const string& output_filename);
extern unsigned long long AddMods(HeadedRecordReader* reader,
                    string out_file,
                    string tmpDir,                    
                    const pb::Header& header,
                    const vector<const pb::Protein*>& proteins,
                    vector<string>& temp_file_name,
                    unsigned long long memory_limit,
                    VariableModTable* var_mod_table);
DECLARE_int32(max_mods);
DECLARE_int32(min_mods);

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
  bool allowDups = Params::GetBool("allow-dups");
  if (FLAGS_min_mods > FLAGS_max_mods) {
    carp(CARP_FATAL, "The value for 'min-mods' cannot be greater than the value "
                     "for 'max-mods'");
  }
  bool sort_on_disk;
  
  unsigned long long memory_limit = Params::GetInt("memory-limit"); //4; // RAM memory limit in GB to be used in in silico protein cleavage.
  
  memory_limit = memory_limit*1000000000/(sizeof(TideIndexPeptide)); //convert the memory limit to number of peptides.
 
  
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

  ofstream* out_target_decoy_list = NULL;  
  if (Params::GetBool("peptide-list")) {
     out_target_decoy_list = create_stream_in_path(make_file_path(
      "tide-index.peptides.txt").c_str(), NULL, overwrite);
  }
  
/*  TODO: Recover the option to generate decoy protein fasta file.
    ofstream* out_decoy_fasta = GeneratePeptides::canGenerateDecoyProteins() ?
    create_stream_in_path(make_file_path(
      "tide-index.decoy.fasta").c_str(), NULL, overwrite) : NULL;
*/    
  string out_proteins = FileUtils::Join(index, "protix");
  string out_peptides = FileUtils::Join(index, "pepix");
  string out_residue_stats = FileUtils::Join(index, "residue_stat");
  string modless_peptides = out_peptides + ".nomods.tmp";
  string peakless_peptides = out_peptides + ".nopeaks.tmp";
  string pathPeptideFile = FileUtils::Join(index, peptideFile);
  string pathMZTabFile = FileUtils::Join(index, tide_index_mzTab_filename_);


  if (create_output_directory(index.c_str(), overwrite) != 0) {
    carp(CARP_FATAL, "Error creating index directory");
  } else if (FileUtils::Exists(out_proteins) ||
             FileUtils::Exists(out_peptides) ||
             FileUtils::Exists(out_residue_stats)) {
    if (overwrite) {
      carp(CARP_DEBUG, "Removing old index file(s)");
      FileUtils::Remove(out_proteins);
      FileUtils::Remove(out_peptides);
      FileUtils::Remove(out_residue_stats);
      FileUtils::Remove(modless_peptides);
      FileUtils::Remove(peakless_peptides);
      FileUtils::Remove(pathPeptideFile);
      FileUtils::Remove(pathMZTabFile);      
    } else {
      carp(CARP_FATAL, "Index file(s) already exist, use --overwrite T or a "
                       "different index name");
    }
  }
  // Define variables for calculating amino acid frequencies (used in tide-search for exact p-value calculation)
  const unsigned int MaxModifiedAAMassBin = MassConstants::ToFixPt(2000.0);   //2000 is the maximum mass of a modified amino acid
  nvAAMassCounterN_ = new unsigned int[MaxModifiedAAMassBin];   //N-terminal amino acids
  nvAAMassCounterC_ = new unsigned int[MaxModifiedAAMassBin];   //C-terminal amino acids
  nvAAMassCounterI_ = new unsigned int[MaxModifiedAAMassBin];   //inner amino acids in the peptides
  memset(nvAAMassCounterN_, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
  memset(nvAAMassCounterC_, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
  memset(nvAAMassCounterI_, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
  cntTerm_ = 0;
  cntInside_ = 0;
  mod_precision_  = Params::GetInt("mod-precision");

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
    carp(CARP_FATAL, "Cannot generate multiple decoys per target in non-shuffled decoy-format!");
  }
  
  // Set up output paths
  if (!FileUtils::Exists(fasta)) {
    carp(CARP_FATAL, "Fasta file %s does not exist", fasta.c_str());
  }

 // Start tide-index
  carp(CARP_INFO, "Reading %s and computing unmodified target peptides...",
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
  headerSource->set_decoy_prefix(Params::GetString("decoy-prefix"));
  HeadedRecordWriter proteinWriter(out_proteins, proteinPbHeader);


  // Generate peptide sequences via in silico cleavage.     
  // Container for the protein header and protein seuqnces.
  ProteinVec vProteinHeaderSequence;  
  
  string proteinHeader;
  std::string proteinSequence;

  FixPt minMassFixPt = MassConstants::ToFixPt(min_mass);
  FixPt maxMassFixPt = MassConstants::ToFixPt(max_mass);
  ifstream file(fasta.c_str(), ifstream::in);
  boost::iostreams::filtering_istreambuf in;
  if (boost::filesystem::extension(fasta) == ".gz") {
    in.push(boost::iostreams::gzip_decompressor());
  }
  in.push(file);
  istream fastaStream(&in);
  unsigned long long invalidPepCnt = 0;
  unsigned long long failedDecoyCnt = 0;

  unsigned long long targetsGenerated = 0;

  long long curProtein = -1;  
  unsigned int pept_file_idx = 0;
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
        carp(CARP_DEBUG, "Ignoring invalid sequence <%s>", std::string(proteinSequence.data()+i->pos_, i->length_).c_str());  
        ++invalidPepCnt;
        continue;
      } else if (pepMass < minMassFixPt || pepMass > maxMassFixPt) {
        // Skip to next peptide if not in mass range
        continue;
      }
      peptide_list.push_back(TideIndexPeptide(pepMass, i->length_, &(pbProtein->residues()), curProtein, i->pos_, -1));
      
      if (peptide_list.size() >= memory_limit){  //reached the memory limit. dump peptides to disk
        // Peptides are being sorted ...
        sort(peptide_list.begin(), peptide_list.end(), less<TideIndexPeptide>());
        
        // ... and dumped in a binary file.
        string pept_file = pathPeptideFile + to_string(pept_file_idx) + ".txt";

        ++pept_file_idx;
        
        dump_peptides_to_binary_file(&peptide_list, pept_file);
        peptide_list.clear();
        vector<TideIndexPeptide> tmp;
        peptide_list.swap(tmp);
  
      }
      ++targetsGenerated;

    }
    if ((curProtein+1) % 10000 == 0) {
      carp(CARP_INFO, "Processed %ld protein sequences", curProtein+1);
    }
  }
  carp(CARP_INFO, "Cleaved %ld protein sequences in total.", curProtein+1);

  sort_on_disk = true;
  if (pept_file_idx == 0) {  //Peptides fit in memory, no need to use disk, sort them in place
    sort(peptide_list.begin(), peptide_list.end(), less<TideIndexPeptide>());
    sort_on_disk = false;
  } else if (peptide_list.size() > 0){ // Some peptides have been already dump on disk, need to dump the remaining ones in peptide_list.
    sort(peptide_list.begin(), peptide_list.end(), less<TideIndexPeptide>());
    string pept_file = pathPeptideFile + to_string(pept_file_idx) + ".txt";
    ++pept_file_idx;
    dump_peptides_to_binary_file(&peptide_list, pept_file);
    peptide_list.clear();
    vector<TideIndexPeptide> tmp;
    peptide_list.swap(tmp);    
  }
    
  if (targetsGenerated == 0) {
    carp(CARP_FATAL, "No target sequences generated.  Is \'%s\' a FASTA file?",
         fasta.c_str());
  }
  if (invalidPepCnt > 0) {
    carp(CARP_INFO, "Ignoring %lu peptide sequences containing unrecognized characters.", invalidPepCnt);
  }
  carp(CARP_INFO, "Generated %lu targets, including duplicates.", targetsGenerated);

  // Prepare the protocol buffer for the peptides.  
  carp(CARP_INFO, "Writing peptides");

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
  
  if (!settings.has_enzyme() || settings.enzyme().empty()) {
    carp(CARP_FATAL, "Enzyme settings error");
  }

  header_no_mods.set_file_type(pb::Header::PEPTIDES);
  header_no_mods.mutable_peptides_header()->set_has_peaks(false);
  header_no_mods.mutable_peptides_header()->set_decoys(decoy_type);

  pb::Peptide pbPeptide;
  unsigned long long count = 0;
  unsigned long long numTargets = 0;
  unsigned long long numDuplicateTargets = 0;
  unsigned long long peptide_cnt = 0;
  
  if (!sort_on_disk && peptide_list.size() == 0)
    carp(CARP_FATAL, "No peptides were generated.");

  unsigned long long numLines = 0;
  TideIndexPeptide currentPeptide;
  TideIndexPeptide duplicatedPeptide;
  TideIndexPeptide* pept_ptr;
  // Filter peptides and keep the unique target peptides and gather the 
  // location of the peptide in other protein sequences 
  vector<FILE*> sortedFiles;
  if (sort_on_disk) {
    //open each file which contain sorted peptides, read the first peptide from each file and put them in a heap.
    for (int i = 0; i < pept_file_idx; ++i) {
      string pept_file = pathPeptideFile + to_string(i) + ".txt";
      FILE* fp = fopen(pept_file.c_str(), "rb");
      sortedFiles.push_back(fp);
      pept_ptr = readNextPeptide(fp, vProteinHeaderSequence, i);  // get the first peptide  
      if (pept_ptr != nullptr) {
        peptide_list.push_back(*pept_ptr);
        delete pept_ptr;
      }
    }
    // Get the lightest peptide from the heap, and remove it from the heap
    std::make_heap(peptide_list.begin(), peptide_list.end(), greater<TideIndexPeptide>());
    currentPeptide = peptide_list.front();   
    int sourceId = currentPeptide.getSourceId();          
    
    // Read another peptide from the disk in order to replace currentPeptide
    pept_ptr = readNextPeptide(sortedFiles[sourceId], vProteinHeaderSequence, sourceId);  // get a peptide  

    std::pop_heap(peptide_list.begin(), peptide_list.end(), greater<TideIndexPeptide>());
    peptide_list.pop_back();   
   
    if (pept_ptr != nullptr) {
      peptide_list.push_back(*pept_ptr);
      push_heap(peptide_list.begin(), peptide_list.end(), greater<TideIndexPeptide>());
      delete pept_ptr;            
    }
  } else {
    currentPeptide = peptide_list[peptide_cnt++];  // get the first peptide  
  }
  
  if (1 == 1) {  // This is needed because we need to destroy the peptideWriter and pbAuxLoc later. Ugly solution :/
    // Create the auxiliary locations header and writer
    HeadedRecordWriter peptideWriter(peptidePbFile, header_no_mods); // put header in outfile  
    bool finished = false;    

    // Decoy generation stuff
    bool success;
    vector<int> decoy_peptide_idx;    
    vector<pb::Peptide> pb_peptides;  
    set<string> peptide_target_str_set;   
    vector<set<string>> peptide_decoy_str_set(numDecoys);    
    string decoy_peptide_str;
    int generateAttemptsMax = 6;    
    pb::Peptide currentPBPeptide;
    getPbPeptide(count++, currentPeptide, currentPBPeptide);      
    pb::AuxLocation pbAuxLoc;
    if (numDecoys == 0) {
      allowDups = true;
    }
    /* The trick to keep the sets target and decoy peptides disjoint is that:
    One does not need to keep all the unique target peptides in the memory
    and check every time whether a decoy peptide already exists as a target.
    It is enought to keep the target in a set (in the memory) peptdes having
    exactly the same mass. It is because the decoy generation does not chage
    the mass of the peptide.
    */    
    double last_mass = -1.0;
    while (!finished) {
      while (true) {
        
        if (sort_on_disk) {
          if (peptide_list.size() == 0) {
            finished = true;
            break;
          }
          
          duplicatedPeptide = peptide_list.front();   
          std::pop_heap (peptide_list.begin(), peptide_list.end(), greater<TideIndexPeptide>());
          peptide_list.pop_back();   
          int sourceId = duplicatedPeptide.getSourceId();          
          pept_ptr = readNextPeptide(sortedFiles[sourceId], vProteinHeaderSequence, sourceId);  // get a peptide  
          numLines++;
          
          if (pept_ptr != nullptr) {
            peptide_list.push_back(*pept_ptr);
            push_heap(peptide_list.begin(), peptide_list.end(), greater<TideIndexPeptide>());
            delete pept_ptr;            
          }
          if (duplicatedPeptide.getMass() < currentPeptide.getMass()) {  // Check if sorting worked properly.
            carp(CARP_INFO, "peptide mass: %lf, subsequent peptide mass %lf", currentPeptide.getMass(), duplicatedPeptide.getMass());
            carp(CARP_FATAL, "Peptides are not sorted correctly. Sorting seems to be failed. Try again and check the free disk space.");
          }
        } else { // sorting in memory. All peptides in peptide_list (in memory) 
          if (peptide_cnt >= peptide_list.size()) {
            finished = true;          
            break;
          }
          duplicatedPeptide = peptide_list[peptide_cnt++];  // get a peptide  
        }
        
        if( duplicatedPeptide == currentPeptide) {
          numDuplicateTargets++;
          carp(CARP_DEBUG, "Skipping duplicate %s.", currentPeptide.getSequence().c_str());
          pb::Location* location = pbAuxLoc.add_location();
          location->set_protein_id(duplicatedPeptide.getProteinId());
          location->set_pos(duplicatedPeptide.getProteinPos());
        } else {
          break;
        }
      }
      
      if (pbAuxLoc.location_size() > 0) {
        pb::AuxLocation* tempAuxLoc = new pb::AuxLocation(pbAuxLoc);
        currentPBPeptide.set_allocated_aux_loc(tempAuxLoc);
        pbAuxLoc.Clear();
      }       
      // Gather the target peptides of having the same mass.
      pb_peptides.push_back(currentPBPeptide);
      
      if (allowDups == false) {
        string target_peptide = vProteinHeaderSequence[currentPBPeptide.first_location().protein_id()]->residues().substr( currentPBPeptide.first_location().pos(), currentPBPeptide.length());
        peptide_target_str_set.insert(target_peptide);
      }
            
      if (duplicatedPeptide.getMass() > currentPeptide.getMass() || allowDups == true || finished == true) {  // Dump peptides to disk: 1. generate decoy permutation idx, and then write them to disk
      
        for (vector<pb::Peptide>::iterator pb_pept_itr = pb_peptides.begin(); pb_pept_itr != pb_peptides.end(); ++pb_pept_itr) {
          
          string target_peptide = vProteinHeaderSequence[(*pb_pept_itr).first_location().protein_id()]->residues().substr( (*pb_pept_itr).first_location().pos(), (*pb_pept_itr).length());
          
          for (int i = 0; i < numDecoys; ++i) {

            shuffle = decoy_type == PEPTIDE_SHUFFLE_DECOYS;
            for (int j = 0; j < generateAttemptsMax; ++j) {
              // Generates a permutation for how generate the decoy peptide from target peptide
              GeneratePeptides::makeDecoyIdx(target_peptide, shuffle, decoy_peptide_idx);
              decoy_peptide_str = target_peptide;

              // Create the decoy peptide sequence, No modifications yet
              for (int k = 0; k < decoy_peptide_idx.size(); ++k) {
                decoy_peptide_str[decoy_peptide_idx[k]] = target_peptide[k];
              }
              // Check if this modified decoy peptide has not been generated yet.
              if (allowDups) {
                success = true;
                break;
              } else {
                // The decoy peptide string with modications can be found in the set of unique peptides?
                success = peptide_target_str_set.find(decoy_peptide_str) == peptide_target_str_set.end(); // generated decoy not found in target peptides
                if (success == true)
                  success = peptide_decoy_str_set[i].find(decoy_peptide_str) == peptide_decoy_str_set[i].end();  // generated decoy not found in decoy peptides
                if (success == true) {    // add decoy peptides to unique decoy peptide set
                  peptide_decoy_str_set[i].insert(decoy_peptide_str);
                  break;   // break the for generateAttemptsMax loop
                }
              }
              shuffle = true; // Failed to generate decoy, so try shuffling in the next attempt.
            }
            if (success == false) {
              carp(CARP_DEBUG, "Failed to generate decoys for sequence %s", target_peptide.c_str());
              ++failedDecoyCnt;
            } else { // Add the decoy permutation idx to the target peptide
              for (int k = 0; k < decoy_peptide_idx.size(); ++k) {
                (*pb_pept_itr).add_decoy_perm_idx(decoy_peptide_idx[k]);
              }
            }
            (*pb_pept_itr).add_decoy_perm_idx(-1);  // Add -1 as a separator between muptiple decoys per target
          }
          // Write the target peptide to disk
          peptideWriter.Write(&(*pb_pept_itr));

          ++numTargets;
          if (numTargets % 1000000 == 0) {
            carp(CARP_INFO, "Wrote %lu unique target peptides", numTargets);
          }      
        }
        // Clear the sets and lists.
        pb_peptides.clear();
        peptide_target_str_set.clear(); 
        for (int i = 0; i < numDecoys; ++i)
          peptide_decoy_str_set[i].clear();
      }
      
      currentPeptide = duplicatedPeptide;
      getPbPeptide(count++, currentPeptide, currentPBPeptide);            
    }
  }
  carp(CARP_DETAILED_INFO, "%lu peptides in file", numLines);
  
  // Release the memory allocated.
  peptide_list.clear();
  vector<TideIndexPeptide> tmp;
  peptide_list.swap(tmp);

  
  carp(CARP_INFO, "Skipped %lu duplicate targets.",
       numDuplicateTargets);
  
  carp(CARP_INFO, "Generated %lu unique target peptides.", numTargets);

  peptidePbFile = peakless_peptides;

  if (sort_on_disk) {
    //Delete intermediate peptarget files.
    for (int i = 0; i < pept_file_idx; ++i) {
      string pept_file = pathPeptideFile + to_string(i) + ".txt";
      FileUtils::Remove(pept_file);
    }
  }
  vector<string> mod_temp_file_names;
  if (need_mods) {
    carp(CARP_INFO, "Computing modified peptides...");
    HeadedRecordReader reader(modless_peptides, NULL, 1024 << 10); // 1024kb buffer
    numTargets = AddMods(&reader, peakless_peptides, Params::GetString("temp-dir"), header_with_mods, vProteinHeaderSequence, mod_temp_file_names, Params::GetInt("memory-limit"), &var_mod_table);
    carp(CARP_INFO, "Created %lu modified and unmodified target peptides.", numTargets);
  } 
  // If no modified peptides are created, then mod_temp_file_names is empty and read the peptides from peptidePbFile

  if (numDecoys > 0) {
      carp(CARP_INFO, "Generating %d decoy(s) per target peptide", numDecoys);
  } else {
      carp(CARP_INFO, "No decoy peptides will be generated");
  }
  unsigned long long decoy_count = 0;
  
  if (numDecoys == 0 && out_target_decoy_list == NULL && need_mods == false) {
    if (rename(peptidePbFile.c_str(), out_peptides.c_str()) != 0)
      carp(CARP_FATAL, "Error creating index files");
    else 
      carp(CARP_INFO, "Pepix file created successfully");
    
  } else {
    
    bool success;
    vector<int> decoy_peptide_idx;
    int startLoc;
    int protein_id;
    int mod_code;
    int decoy_index;
    int mod_index;
    int unique_delta;
    double delta;
    double mass;

    string target_peptide_with_mods;
    string decoy_peptide_with_mods;
    int prot_id, pos, len;

    pb::Header new_header;
    new_header.set_file_type(pb::Header::PEPTIDES);
    pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
    subheader->CopyFrom(header_with_mods.peptides_header());
    subheader->set_has_peaks(true);
    source = new_header.add_source();
    source->mutable_header()->CopyFrom(header_with_mods);
    HeadedRecordWriter writer(out_peptides, new_header);

    // Read peptides protocol buffer file
    int mass_precision = Params::GetInt("mass-precision");
    int mod_precision = Params::GetInt("mod-precision");

    pb::Peptide current_pb_peptide_;
    string decoy_peptide_str;
    pb::Peptide temp_pb_peptide;
    const pb::Protein* protein;
    string pepmass_str;
    string pos_str;
    string mod_str;
    int mod_pos_offset;

    if (out_target_decoy_list) {
      *out_target_decoy_list << "target\t";
      if (numDecoys > 0)
        *out_target_decoy_list << "decoy(s)\t";
      *out_target_decoy_list << "mass\tproteins" << std::endl;
    }
    // Go over the (modified and unmodified) peptides from the protocol buffer and generate decoy peptides 
    bool done = false;
    
    // The duplicated target peptides have already been filtered out, no need to check it again iff decoys are not gerenated.
    if (numDecoys == 0) {
      allowDups = true;
    }
    
    // Prepare a queue (pool) to merge the modified peptide files
    vector<pb::Peptide> pb_peptide_pool;
    vector<RecordReader*> readers;
    int source_id = 0;
    pb::Header aaf_peptides_header;
    HeadedRecordReader aaf_peptide_reader(peptidePbFile, &aaf_peptides_header);
    if (aaf_peptides_header.file_type() != pb::Header::PEPTIDES ||
        !aaf_peptides_header.has_peptides_header()) {
      carp(CARP_FATAL, "Error reading index (%s)", peptidePbFile.c_str());
      }

    RecordReader* reader_;
    reader_ = aaf_peptide_reader.Reader();
    
    if (!mod_temp_file_names.empty()) {
      for (vector<string>::iterator i = mod_temp_file_names.begin(); i != mod_temp_file_names.end(); ++i) {
        RecordReader* reader = new RecordReader(*i, 1024 << 10);
        CHECK(reader->OK());
        readers.push_back(reader);
        if (!reader->Done()) {
           reader->Read(&current_pb_peptide_);
           CHECK(reader->OK());         
           current_pb_peptide_.set_decoy_index(source_id);   //use this field to temporarily indicate the origin file of a peptide. 
           pb_peptide_pool.push_back(current_pb_peptide_);         
         }
        source_id++;
        carp(CARP_DEBUG, "temp modification file %s", (*i).c_str());
      }
    } else {
      readers.push_back(reader_);
      CHECK(reader_->OK());
      if (!reader_->Done()) {
        reader_->Read(&current_pb_peptide_);
        CHECK(reader_->OK());         
        current_pb_peptide_.set_decoy_index(0);   //use this field to temporarily indicate the origin file of a peptide. 
        pb_peptide_pool.push_back(current_pb_peptide_);         
      }
    }
    // Get the lightest peptide from the heap to the front
    std::make_heap(pb_peptide_pool.begin(), pb_peptide_pool.end(), PbPeptideSortGreater());    
    
    CHECK(writer.OK());
    
    peptide_cnt = 0;
    while (!done) {
      
      // Check if there is still a peptide in the pool.
      if (pb_peptide_pool.size() == 0) {
        break;
      }
      // Here we do the modified peptide merge.
      // Get the peptide from the pool with the smallest mass
      current_pb_peptide_ = pb_peptide_pool.front();  
      std::pop_heap(pb_peptide_pool.begin(), pb_peptide_pool.end(), PbPeptideSortGreater());
      pb_peptide_pool.pop_back();         
      
      // Load another peptide into the pool from the porotocol files.
      source_id = current_pb_peptide_.decoy_index();
      if ( !readers[source_id]->Done() ) {
        readers[source_id]->Read(&temp_pb_peptide);
        CHECK(readers[source_id]->OK());
        temp_pb_peptide.set_decoy_index(source_id);  // We use the decoy index in order to keep track the source file ID of the peptide
        pb_peptide_pool.push_back(temp_pb_peptide);
        std::push_heap(pb_peptide_pool.begin(), pb_peptide_pool.end(), PbPeptideSortGreater());  // maintain heap
      }
      current_pb_peptide_.set_decoy_index(-1);  //restore the source id and use the decoy index as planned

      // Get the amino acid frequencies from the peptides
      getAAFrequencies(current_pb_peptide_, vProteinHeaderSequence);

      // Get the peptide sequence with modifications
      if (out_target_decoy_list) {
        target_peptide_with_mods = getModifiedPeptideSeq(&current_pb_peptide_, &vProteinHeaderSequence);
        *out_target_decoy_list << target_peptide_with_mods;
      }

      protein_id = current_pb_peptide_.first_location().protein_id();
      startLoc = current_pb_peptide_.first_location().pos();

      if (numDecoys > 0) {  // Get peptide sequence without mods
        if (out_target_decoy_list) {
          *out_target_decoy_list << '\t';
        }
        string target_peptide = vProteinHeaderSequence[protein_id]->residues().substr(startLoc, current_pb_peptide_.length());
        string decoy_peptide_str_with_mods;

        //  Generate a decoy peptide:
        protein = vProteinHeaderSequence[protein_id];
        int decoy_permutation_idx = 0;
        int perm_idx = 0;
        for (int i = 0; i < numDecoys; ++i) {
          if (i > 0 && out_target_decoy_list) {
            *out_target_decoy_list << ',';           
          }
          decoy_peptide_idx.clear();
          while (true) {
            perm_idx = current_pb_peptide_.decoy_perm_idx(decoy_permutation_idx++);
            if (perm_idx == -1){
              break;
            }
            decoy_peptide_idx.push_back(perm_idx);
          }
          if (decoy_peptide_idx.empty() == true)
            continue;
          decoy_peptide_str = target_peptide;

          // Create the decoy peptide sequence without modifications
          for (int k = 0; k < decoy_peptide_idx.size(); ++k) {
            decoy_peptide_str[decoy_peptide_idx[k]] = target_peptide[k];
          }

          // Create a protocol buffer peptide object for the decoy peptide. Note that the decoy peptide may contain modifications.
          pb::Peptide decoy_current_pb_peptide_ = current_pb_peptide_;
          if (current_pb_peptide_.modifications_size() > 0) {
            decoy_current_pb_peptide_.clear_modifications();
            for (int m = 0; m < current_pb_peptide_.modifications_size(); ++m) {
              mod_code = current_pb_peptide_.modifications(m);

              var_mod_table.DecodeMod(mod_code, &mod_index, &unique_delta);
              decoy_index = decoy_peptide_idx[mod_index];
              mod_code = var_mod_table.EncodeMod(decoy_index, unique_delta);
              decoy_current_pb_peptide_.add_modifications(mod_code);
            }
          }
          decoy_current_pb_peptide_.set_id(numTargets + decoy_count++);
          decoy_current_pb_peptide_.clear_decoy_sequence();
          decoy_current_pb_peptide_.set_decoy_sequence(decoy_peptide_str);
          decoy_current_pb_peptide_.set_decoy_index(i);
          decoy_current_pb_peptide_.clear_decoy_perm_idx();
          CHECK(writer.Write(&decoy_current_pb_peptide_));

          //report the decoy peptide if needed.
          if (out_target_decoy_list) {
            decoy_peptide_str_with_mods = getModifiedPeptideSeq(&decoy_current_pb_peptide_,  &vProteinHeaderSequence);
            *out_target_decoy_list << decoy_peptide_str_with_mods.c_str();
          }
        }
      }
      // Print 1) the peptide neutral mass, 2) protein header of origin and 3) the locations of the target peptides
      if (out_target_decoy_list) {
        string pepmass_str = StringUtils::ToString(current_pb_peptide_.mass(), mass_precision);
        *out_target_decoy_list << '\t' << pepmass_str;

        pos_str = StringUtils::ToString(startLoc + 1, 1);
        string proteinNames = vProteinHeaderSequence[protein_id]->name() + '(' + pos_str + ')';
        if (current_pb_peptide_.has_aux_loc() == true) {
          const pb::AuxLocation& aux_loc = current_pb_peptide_.aux_loc();
          for (int i = 0; i < aux_loc.location_size(); ++i) {
            const pb::Location& location = aux_loc.location(i);
            protein = vProteinHeaderSequence[location.protein_id()];
            pos_str = StringUtils::ToString(location.pos() + 1, 1);
            proteinNames += ',' + protein->name() + '(' + pos_str + ')';
          }
        }
        *out_target_decoy_list << '\t' << proteinNames << endl;
      }
      current_pb_peptide_.clear_decoy_perm_idx();
      CHECK(writer.Write(&current_pb_peptide_));      
      ++peptide_cnt;
      if (peptide_cnt % 10000000 == 0) {
        carp(CARP_INFO, "Wrote %lu target and their corresponding decoy peptides", peptide_cnt);
      }
    }
    for (vector<string>::iterator i = mod_temp_file_names.begin(); i != mod_temp_file_names.end(); ++i) {
      unlink((*i).c_str());
    }

    if (out_target_decoy_list) {
      out_target_decoy_list->close();
      delete out_target_decoy_list;
    }
    if (failedDecoyCnt > 0) {
      carp(CARP_INFO, "Failed to generate decoys for %lu low complexity peptides.", failedDecoyCnt);
    }
  }
  // Write the amino acid frequencies
  vector<double> dAAFreqN;
  vector<double> dAAFreqI;
  vector<double> dAAFreqC;
  vector<double> dAAMass;

  unsigned int uiUniqueMasses = 0;
  for (int i = 0; i < MaxModifiedAAMassBin; ++i) {
    if (nvAAMassCounterN_[i] || nvAAMassCounterI_[i] || nvAAMassCounterC_[i]) {
      ++uiUniqueMasses;
      dAAMass.push_back(MassConstants::ToDouble(i));
      dAAFreqN.push_back((double)nvAAMassCounterN_[i] / cntTerm_);
      dAAFreqI.push_back((double)nvAAMassCounterI_[i] / cntInside_);
      dAAFreqC.push_back((double)nvAAMassCounterC_[i] / cntTerm_);
    }
  }
  RecordWriter residue_stat_wirter = RecordWriter(out_residue_stats);
  CHECK(residue_stat_wirter.OK());  
  for (int i = 0; i < dAAMass.size(); ++i){
    pb::ResidueStats last_residue_stat;
    last_residue_stat.set_aamass(dAAMass[i]);
    last_residue_stat.set_aafreqn(dAAFreqN[i]);
    last_residue_stat.set_aafreqi(dAAFreqI[i]);
    last_residue_stat.set_aafreqc(dAAFreqC[i]);
    string aa_str = mMass2AA_[dAAMass[i]];
    last_residue_stat.set_aa_str(aa_str);
    CHECK(residue_stat_wirter.Write(&last_residue_stat));
    // printf("%lf, %lf, %lf, %lf, %s\n", dAAMass[i], dAAFreqN[i], dAAFreqI[i], dAAFreqC[i], aa_str.c_str());

  }

  carp(CARP_INFO, "Generated %lu target peptides.", peptide_cnt);
  carp(CARP_INFO, "Generated %lu decoy peptides.", decoy_count);
  carp(CARP_INFO, "Generated %lu peptides in total.", peptide_cnt + decoy_count);
  
  // Recover stderr
  cerr.rdbuf(old);
 
  FileUtils::Remove(modless_peptides);
  FileUtils::Remove(peakless_peptides);
  
  delete nvAAMassCounterN_;   //N-terminal amino acids
  delete nvAAMassCounterC_;   //C-terminal amino acids
  delete nvAAMassCounterI_;   //inner amino acids in the peptides

  // Dump the parameters in mzTAB format
  try {
    int cnt = 1;
    ofstream mzTabStream(pathMZTabFile);
   
    mzTabStream << "MTD\tsoftware[1]\t[MS, MS:1002575, tide-index, " << CRUX_VERSION << "]\n";    
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tauto-modifications-spectra = " << Params::GetString("auto-modifications-spectra") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tclip-nterm-methionine = " << Params::GetString("clip-nterm-methionine") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tisotopic-mass = " << Params::GetString("isotopic-mass") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmax-length = " << Params::GetInt("max-length") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmax-mass = " << Params::GetDouble("max-mass") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmin-length = " << Params::GetInt("min-length") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmin-mass = " << Params::GetDouble("min-mass") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tcterm-peptide-mods-spec = " << Params::GetString("cterm-peptide-mods-spec") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tcterm-protein-mods-spec = " << Params::GetString("cterm-protein-mods-spec") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmax-mods = " << Params::GetInt("max-mods") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmin-mods = " << Params::GetInt("min-mods") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmod-precision = " << Params::GetInt("mod-precision") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmods-spec = " << Params::GetString("mods-spec") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tnterm-peptide-mods-spec = " << Params::GetString("nterm-peptide-mods-spec") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tnterm-protein-mods-spec = " << Params::GetString("nterm-protein-mods-spec") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tauto-modifications = " << Params::GetString("auto-modifications") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tallow-dups = " << Params::GetString("allow-dups") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tdecoy-format = " << Params::GetString("decoy-format") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tkeep-terminal-aminos = " << Params::GetString("keep-terminal-aminos") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tnum-decoys-per-target = " << numDecoys <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tseed = " << Params::GetString("seed") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tcustom-enzyme = " << Params::GetString("custom-enzyme") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tdigestion = " << Params::GetString("digestion") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tenzyme = " << Params::GetString("enzyme") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmissed-cleavages = " << Params::GetInt("missed-cleavages") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tdecoy-prefix = " << Params::GetString("decoy-prefix") <<"\n";
    mzTabStream << "MTD\tsoftware[1]-setting[" << cnt++ << "]\tmass-precision = " << Params::GetInt("mass-precision") <<"\n";

    mzTabStream.close();

  } catch (...){
    carp(CARP_INFO, "mzTab file was not created");
  }
  
  // Recover stderr
  cerr.rdbuf(old);


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
    "keep-terminal-aminos",  //TODO: remove this option. handled in GeneratePeptides.Cpp
    "mass-precision",
    "max-length",
    "max-mass",
    "max-mods",
    "memory-limit",
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
      carp(CARP_FATAL, "Spectrum files must be specified with the 'auto-modifications-spectra' "
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
  string default_cysteine = "C[Unimod:4]"; //+ StringUtils::ToString(CYSTEINE_DEFAULT);
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
  // stringstream mod_stream;
  int mod_pos_offset = 0;
  int index;
  double delta;
  int modPrecision = Params::GetInt("mod-precision");

  const pb::Location& location = peptide->first_location();
  const pb::Protein* protein = proteins->at(location.protein_id());
  string mod_str;
  string seq_with_mods ;
  // Get peptide sequence without mods, 
  if (peptide->has_decoy_sequence()){  // decoy or target
    seq_with_mods = peptide->decoy_sequence();
  } else {
    seq_with_mods = protein->residues().substr(location.pos(), peptide->length());
  }

  if (peptide->has_nterm_mod()){ // Handle N-terminal modifications
    MassConstants::DecodeMod(ModCoder::Mod(peptide->nterm_mod()), &index, &delta);
    mod_str = "[" + StringUtils::ToString(delta, modPrecision) + "]-";
    seq_with_mods.insert(0, mod_str);
    mod_pos_offset += mod_str.length();
  }

  int num_mods = peptide->modifications_size();
  if (num_mods > 0) {
    vector<int> mod;
    
    for (int i = 0; i < num_mods; ++i) {
      mod.push_back(peptide->modifications(i));
    }
    
    sort(mod.begin(), mod.end());
    
    for (int i = 0; i < num_mods; ++i) {
      int index;
      double delta;
      MassConstants::DecodeMod(mod[i], &index, &delta);
      mod_str = "[" + StringUtils::ToString(delta, modPrecision) + "]";
      seq_with_mods.insert(index + 1 + mod_pos_offset, mod_str);
      mod_pos_offset += mod_str.length();
    }
  }
  if (peptide->has_cterm_mod()){  // Handle C-terminal modifications
    MassConstants::DecodeMod(ModCoder::Mod(peptide->cterm_mod()), &index, &delta);
    mod_str = "-[" + StringUtils::ToString(delta, modPrecision) + "]";
    seq_with_mods.insert(index + 1 + mod_pos_offset, mod_str);
    mod_pos_offset += mod_str.length();
  }

  return seq_with_mods;  
}

TideIndexApplication::TideIndexPeptide* TideIndexApplication::readNextPeptide(FILE* fp, ProteinVec& vProteinHeaderSequence, int sourceId) {
  
  FixPt pepMass;
  int prot_id;
  int pos;
  int len;
  int ret;
  ret = fread(&pepMass, sizeof(FixPt), 1, fp);  
  if (ret == 0)
    return nullptr; 
  ret = fread(&prot_id, sizeof(int), 1, fp);  
  if (ret == 0)
    return nullptr; 
  ret = fread(&pos, sizeof(int), 1, fp);  
  if (ret == 0)
    return nullptr; 
  ret = fread(&len, sizeof(int), 1, fp);  
  if (ret == 0)
    return nullptr; 
  
  int decoyIdx = -1; // -1 if not a decoy; There are no decoy peptides generated at this point

  const string& proteinSequence = vProteinHeaderSequence[prot_id]->residues();
    
  TideIndexPeptide* pepTarget = new TideIndexPeptide(pepMass, len, &proteinSequence, prot_id, pos, decoyIdx, sourceId);
  
  return pepTarget;
}

void TideIndexApplication::dump_peptides_to_binary_file(vector<TideIndexPeptide> *peptide_list, string pept_file) {
        
  FILE* fp = fopen(pept_file.c_str(), "wb");  // Peptides stored in this file to be sorted on disk.
  FixPt pepMass;
  int prot_id;
  int len;
  int pos;
  int ret;
  for (vector<TideIndexPeptide>::iterator pept_itr = peptide_list->begin(); pept_itr != peptide_list->end(); ++pept_itr) {
    pepMass = (*pept_itr).getFixPtMass();
    prot_id = (*pept_itr).getProteinId();
    len = (*pept_itr).getLength();
    pos = (*pept_itr).getProteinPos();
    
    ret = fwrite(&pepMass, sizeof(FixPt), 1, fp);
    if (ret == 0) {
      carp(CARP_FATAL, "Error while writting to disk. ");
    }
    ret = fwrite(&prot_id, sizeof(int), 1, fp);
    if (ret == 0) {
      carp(CARP_FATAL, "Error while writting to disk. ");
    }
    ret = fwrite(&pos, sizeof(int), 1, fp);
    if (ret == 0) {
      carp(CARP_FATAL, "Error while writting to disk. ");
    }
    ret = fwrite(&len, sizeof(int), 1, fp);
    if (ret == 0) {
      carp(CARP_FATAL, "Error while writting to disk. ");
    }
  }   
  fclose(fp);    

}

void TideIndexApplication::getAAFrequencies(pb::Peptide& current_pb_peptide, ProteinVec& vProteinHeaderSequence){
  unsigned int len;
  unsigned int i;
  unsigned int residue_bin;  
  string tempAA;

  Peptide peptide(current_pb_peptide, vProteinHeaderSequence);
  vector<double> residue_masses = peptide.getAAMasses(); //retrieves the amino acid masses, modifications included
  string peptide_seq = peptide.Seq();
  len = current_pb_peptide.length();
  vector<double> residue_mods(len, 0);  // Initialize a vecotr of peptide length  with zeros.   

  // Handle variable modifications
  if (current_pb_peptide.has_nterm_mod()){ // Handle N-terminal modifications
    int index;
    double delta;
    MassConstants::DecodeMod(ModCoder::Mod(current_pb_peptide.nterm_mod()), &index, &delta);
    residue_mods[index] = delta;
  }

  for (i = 0; i < current_pb_peptide.modifications_size(); ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(current_pb_peptide.modifications(i), &index, &delta);
    residue_mods[index] = delta;
  }
  
  if (current_pb_peptide.has_cterm_mod()){  // Handle C-terminal modifications
    int index;
    double delta;
    MassConstants::DecodeMod(ModCoder::Mod(current_pb_peptide.cterm_mod()), &index, &delta);
    residue_mods[index] = delta;
  }

  // count AA masses
  residue_bin = MassConstants::ToFixPt(residue_masses[0]);  
  ++nvAAMassCounterN_[residue_bin];  // N-temrianl
  if (nvAAMassCounterN_[residue_bin] == 1){
    tempAA = peptide_seq[0];
    if (residue_mods[0] != 0) {
      tempAA += "[" + StringUtils::ToString(residue_mods[0], mod_precision_) + ']';
    }
    mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
  }
  for (i = 1; i < len-1; ++i) {
    residue_bin = MassConstants::ToFixPt(residue_masses[i]);
    ++nvAAMassCounterI_[residue_bin];  // non-terminal
    if (nvAAMassCounterI_[residue_bin] == 1){
      tempAA = peptide_seq[i];
      if (residue_mods[i] != 0) {
        tempAA += "[" + StringUtils::ToString(residue_mods[i], mod_precision_) + ']';
      }
      mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
    }		
    ++cntInside_;
  }
  residue_bin = MassConstants::ToFixPt(residue_masses[len - 1]);
  ++nvAAMassCounterC_[residue_bin];  // C-temrinal
  if (nvAAMassCounterC_[residue_bin] == 1){
    tempAA = peptide_seq[len - 1];
    if (residue_mods[len - 1] != 0) {
      tempAA += "[" + StringUtils::ToString(residue_mods[len - 1], mod_precision_) + ']';
    }
    mMass2AA_[MassConstants::ToDouble(residue_bin)] = tempAA;
  }
  ++cntTerm_;
}

