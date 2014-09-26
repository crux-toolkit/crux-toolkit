#include <cstdio>
#include <fstream>
#include "carp.h"
#include "CarpStreamBuf.h"
#include "GenerateDecoys.h"
#include "TideIndexApplication.h"
#include "TideMatchSet.h"
#include "tide/modifications.h"
#include "tide/records_to_vector-inl.h"

#ifdef _MSC_VER
#include <io.h>
#endif

extern HASH_T* parameters;
extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
                                const string& input_filename,
                                const string& output_filename);
//extern void AddMods(HeadedRecordReader* reader,
//                    string out_file,
//                    const pb::Header& header,
//                    const vector<const pb::Protein*>& proteins);
extern void AddMods(HeadedRecordReader* reader,
                    string out_file,
                    const pb::Header& header,
                    const vector<const pb::Protein*>& proteins, VariableModTable& var_mod_table);
DECLARE_int32(max_mods);
DECLARE_int32(min_mods);
DECLARE_string(tmpfile_prefix);

TideIndexApplication::TideIndexApplication() {
}

TideIndexApplication::~TideIndexApplication() {
}

int TideIndexApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "decoy-format",
    "keep-terminal-aminos",
    "decoy-prefix",
    "enzyme",
    "custom-enzyme",
    "digestion",
    "missed-cleavages",
    "max-length",
    "max-mass",
    "min-length",
    "min-mass",
    "monoisotopic-precursor",
    "mods-spec",
    "cterm-peptide-mods-spec",
    "nterm-peptide-mods-spec",
    "cterm-protein-mods-spec",
    "nterm-protein-mods-spec",
    "max-mods",
    "min-mods",
    "output-dir",
    "overwrite",
    "peptide-list",
    "parameter-file",
    "seed",
    "clip-nterm-methionine",
    "verbosity"
  };

  // Crux command line parsing
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "protein fasta file",
    "index name"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-index...");

  // Reroute stderr
  CarpStreamBuf buffer;
  streambuf* old = cerr.rdbuf();
  cerr.rdbuf(&buffer);

  // Build command line string
  string cmd_line = "crux tide-index";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  FLAGS_tmpfile_prefix = make_file_path("modified_peptides_partial_");

  // Get options
  double min_mass = get_double_parameter("min-mass");
  double max_mass = get_double_parameter("max-mass");
  int min_length = get_int_parameter("min-length");
  int max_length = get_int_parameter("max-length");
  bool monoisotopic_precursor = get_boolean_parameter("monoisotopic-precursor");
  FLAGS_max_mods = get_int_parameter("max-mods");
  FLAGS_min_mods = get_int_parameter("min-mods");
  if (FLAGS_min_mods > FLAGS_max_mods) {
		carp(CARP_FATAL, "The value for 'min-mods' cannot be greater than the value "
                     "for 'max-mods'");
	}
  MASS_TYPE_T mass_type = (monoisotopic_precursor) ? MONO : AVERAGE;
  int missed_cleavages = get_int_parameter("missed-cleavages");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  ENZYME_T enzyme_t = get_enzyme_type_parameter("enzyme");
  char* enzymePtr = enzyme_type_to_string(enzyme_t);
  string enzyme(enzymePtr);
  free(enzymePtr);
  if (enzyme == "no-enzyme") {
    enzyme = "none";
  } else if (digestion != FULL_DIGEST && digestion != PARTIAL_DIGEST) {
    carp(CARP_FATAL, "'digestion' must be 'full-digest' or 'partial-digest'");
  }

  VariableModTable var_mod_table;
  var_mod_table.ClearTables();
  //parse regular amino acid modifications
  string mods_spec = get_string_parameter_pointer("mods-spec");
  carp(CARP_DEBUG, "mods_spec='%s'", mods_spec.c_str());
  if (!var_mod_table.Parse(mods_spec.c_str())) {
    carp(CARP_FATAL, "Error parsing mods");
  }
  //parse terminal modifications
  mods_spec = get_string_parameter_pointer("cterm-peptide-mods-spec");
  if (!mods_spec.empty())
   if (!var_mod_table.Parse(mods_spec.c_str(), CTPEP)) {
    carp(CARP_FATAL, "Error parsing c-terminal peptide mods");
  }
  mods_spec = get_string_parameter_pointer("nterm-peptide-mods-spec");
  if (!mods_spec.empty())
   if (!var_mod_table.Parse(mods_spec.c_str(), NTPEP)) {
    carp(CARP_FATAL, "Error parsing n-terminal peptide mods");
  }
  mods_spec = get_string_parameter_pointer("cterm-protein-mods-spec");
  if (!mods_spec.empty())
   if (!var_mod_table.Parse(mods_spec.c_str(), CTPRO)) {
    carp(CARP_FATAL, "Error parsing c-terminal protein mods");
  }
  mods_spec = get_string_parameter_pointer("nterm-protein-mods-spec");
  if (!mods_spec.empty())
   if (!var_mod_table.Parse(mods_spec.c_str(), NTPRO)) {
    carp(CARP_FATAL, "Error parsing n-terminal protein mods");
  }

  var_mod_table.SerializeUniqueDeltas();

  if (!MassConstants::Init(var_mod_table.ParsedModTable(), 0, 0)) {
    carp(CARP_FATAL, "Error in MassConstants::Init");
  }

  DECOY_TYPE_T decoy_type = get_tide_decoy_type_parameter("decoy-format");
  string decoyPrefix = get_string_parameter_pointer("decoy-prefix");

  // Set up output paths
  string fasta = get_string_parameter_pointer("protein fasta file");
  string index = get_string_parameter_pointer("index name");
  bool overwrite = get_boolean_parameter("overwrite");

  if (!file_exists(fasta)) {
    carp(CARP_FATAL, "Fasta file %s does not exist", fasta.c_str());
  }

  string out_proteins = index + "/" + "protix";
  string out_peptides = index + "/" + "pepix";
  string out_aux = index + "/" + "auxlocs";
  string modless_peptides = out_peptides + ".nomods.tmp";
  string peakless_peptides = out_peptides + ".nopeaks.tmp";
  ofstream* out_target_list = NULL;
  ofstream* out_decoy_list = NULL;
  if (get_boolean_parameter("peptide-list")) {
    out_target_list = create_stream_in_path(make_file_path(
      "tide-index.peptides.target.txt").c_str(), NULL, overwrite);
    if (decoy_type != NO_DECOYS) {
      out_decoy_list = create_stream_in_path(make_file_path(
        "tide-index.peptides.decoy.txt").c_str(), NULL, overwrite);
    }
  }
  ofstream* out_decoy_fasta = GenerateDecoys::canGenerateDecoyProteins() ?
    create_stream_in_path(make_file_path(
      "tide-index.decoy.fasta").c_str(), NULL, overwrite) : NULL;

  if (create_output_directory(index.c_str(), overwrite) != 0) {
    carp(CARP_FATAL, "Error creating index directory");
  } else if (file_exists(out_proteins) ||
             file_exists(out_peptides) ||
             file_exists(out_aux)) {
    if (overwrite) {
      carp(CARP_DEBUG, "Cleaning old index file(s)");
      remove(out_proteins.c_str());
      remove(out_peptides.c_str());
      remove(out_aux.c_str());
      remove(modless_peptides.c_str());
      remove(peakless_peptides.c_str());
    } else {
      carp(CARP_FATAL, "Index file(s) already exist, use --overwrite T or a "
                       "different index name");
    }
  }

  // Start tide-index
  carp(CARP_INFO, "Reading %s and computing unmodified peptides...",
       fasta.c_str());
  pb::Header proteinPbHeader;
  vector<TideIndexPeptide> peptideHeap;
  vector<string*> proteinSequences;
  fastaToPb(cmd_line, enzyme_t, digestion, missed_cleavages, min_mass, max_mass,
            min_length, max_length, mass_type, decoy_type, fasta, out_proteins,
            proteinPbHeader, peptideHeap, proteinSequences, out_decoy_fasta);

  pb::Header header_with_mods;

  // Set up peptides header
  pb::Header_PeptidesHeader& pep_header = *(header_with_mods.mutable_peptides_header());
  pep_header.Clear();
  pep_header.set_min_mass(min_mass);
  pep_header.set_max_mass(max_mass);
  pep_header.set_min_length(min_length);
  pep_header.set_max_length(max_length);
  pep_header.set_monoisotopic_precursor(monoisotopic_precursor);
  pep_header.set_enzyme(enzyme);
  if (enzyme != "none") {
    pep_header.set_full_digestion(digestion == FULL_DIGEST);
    pep_header.set_max_missed_cleavages(missed_cleavages);
  }
  pep_header.mutable_mods()->CopyFrom(*(var_mod_table.ParsedModTable()));

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

  bool need_mods = var_mod_table.Unique_delta_size()>0;

  string basic_peptides = need_mods ? modless_peptides : peakless_peptides;
  carp(CARP_DETAILED_DEBUG, "basic_peptides=%s", basic_peptides.c_str());

  writePeptidesAndAuxLocs(peptideHeap, basic_peptides, out_aux, header_no_mods);
  // Do some clean up
  for (vector<string*>::iterator i = proteinSequences.begin();
       i != proteinSequences.end();
       ++i) {
    delete *i;
  }
  vector<TideIndexPeptide>().swap(peptideHeap);

  vector<const pb::Protein*> proteins;
  if (!ReadRecordsToVector<pb::Protein>(&proteins, out_proteins)) {
    carp(CARP_FATAL, "Error reading proteins file");
  }

  if (need_mods) {
    carp(CARP_INFO, "Computing modified peptides...");
    HeadedRecordReader reader(modless_peptides, NULL, 1024 << 10); // 1024kb buffer
    AddMods(&reader, peakless_peptides, header_with_mods, proteins, var_mod_table);
  }

  if (out_target_list) {
    // Write peptide lists
    carp(CARP_INFO, "Writing peptide lists...");
    // Initialize modification variables
    stringstream mod_stream;
    mod_stream.precision(get_int_parameter("mod-precision"));
    int mod_index;
    double mod_delta;
    // This set holds target peptide strings
    set<string> targetPepStrs;
    // This vector holds decoy peptide strings and masses
    vector< pair<string, double> > decoyPepStrs;
    // Read peptides pb file
    vector<const pb::Peptide*> peptides;
//    if (!ReadRecordsToVector<pb::Peptide>(&peptides, peakless_peptides)) {
//      carp(CARP_FATAL, "Error reading peptides file");
//    }
    vector<const pb::AuxLocation*> locations;
    if (!ReadRecordsToVector<pb::AuxLocation>(&locations, out_aux)) {
      carp(CARP_FATAL, "Error reading auxlocs file");
    }
    cout << "read all peptides with location" << endl;
    // Iterate over all pb peptides
    unsigned int writeCountTargets = 0, writeCountDecoys = 0;

    HeadedRecordReader reader(peakless_peptides, NULL);
    while (!reader.Done()) {
      pb::Peptide* protobuf = new pb::Peptide;
      reader.Read(protobuf);
//      vec->push_back(protobuf);
//    }
//    for (vector<const pb::Peptide*>::const_iterator i = peptides.begin();
//         i != peptides.end();
//         ++i) {
//      const pb::Peptide* peptide = *i;
      pb::Peptide* peptide = protobuf;
      const pb::Location& location = peptide->first_location();
      const pb::Protein* protein = proteins[location.protein_id()];
      bool writeTarget = true;
      bool writeDecoy = false;

      if (out_decoy_list) {
        writeDecoy = peptide->is_decoy();
        writeTarget = !writeDecoy;
      }
      // Get peptide sequence without mods
      string pep_str = protein->residues().substr(
        location.pos(), peptide->length());
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
      // Iterate over mod indices in reverse order
      for (set<int>::const_reverse_iterator j = mod_indices.rbegin();
           j != mod_indices.rend();
           ++j) {
        // Insert the modification string into the peptide sequence
        mod_stream << '[' << mod_map[*j] << ']';
        pep_str.insert(*j + 1, mod_stream.str());
        mod_stream.str("");
      }
      if (writeTarget) {
        // This is a target, output it
        targetPepStrs.insert(pep_str);
        *out_target_list << pep_str << '\t' << peptide->mass() << endl;
        ++writeCountTargets;
      }
      if (writeDecoy) {
        // This is a decoy, save it to output later
        decoyPepStrs.push_back(make_pair(pep_str, peptide->mass()));
        ++writeCountDecoys;
      }
      delete peptide;
    }

    // Iterate over saved decoys and output them
    for (vector< pair<string, double> >::iterator i = decoyPepStrs.begin();
         i != decoyPepStrs.end();
         ++i) {
      *out_decoy_list << i->first << '\t' << i->second;
      if (targetPepStrs.find(i->first) != targetPepStrs.end()) {
        *out_decoy_list << "\t*";
      }
      *out_decoy_list << endl;
    }

    // Clean up peptides
    for (vector<const pb::Peptide*>::iterator i = peptides.begin();
         i != peptides.end();
         ++i) {
      delete *i;
    }

    // Close and clean up streams
    if (out_decoy_list) {
      out_decoy_list->close();
      delete out_decoy_list;
    }
    out_target_list->close();
    delete out_target_list;

    carp(CARP_DEBUG, "Wrote %d targets and %d decoys to peptide list",
         writeCountTargets, writeCountDecoys);
  }

  carp(CARP_INFO, "Precomputing theoretical spectra...");
  AddTheoreticalPeaks(proteins, peakless_peptides, out_peptides);

  // Clean up
  for (vector<const pb::Protein*>::iterator i = proteins.begin();
       i != proteins.end();
       ++i) {
    delete *i;
  }

  // Recover stderr
  cerr.rdbuf(old);

  return 0;
}

string TideIndexApplication::getName() {
  return "tide-index";
}

string TideIndexApplication::getDescription() {
  return "Create an index for all peptides in a fasta file, for use in "
         "subsequent calls to tide-search.";
}

bool TideIndexApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideIndexApplication::getCommand() {
  return TIDE_INDEX_COMMAND;
}

void TideIndexApplication::fastaToPb(
  const string& commandLine,
  const ENZYME_T enzyme,
  const DIGEST_T digestion,
  int missedCleavages,
  FLOAT_T minMass,
  FLOAT_T maxMass,
  int minLength,
  int maxLength,
  MASS_TYPE_T massType,
  DECOY_TYPE_T decoyType,
  const string& fasta,
  const string& proteinPbFile,
  pb::Header& outProteinPbHeader,
  vector<TideIndexPeptide>& outPeptideHeap,
  vector<string*>& outProteinSequences,
  ofstream* decoyFasta
) {
  string decoyPrefix = get_string_parameter_pointer("decoy-prefix");
  outProteinPbHeader.Clear();
  outProteinPbHeader.set_file_type(pb::Header::RAW_PROTEINS);
  outProteinPbHeader.set_command_line(commandLine);
  pb::Header_Source* headerSource = outProteinPbHeader.add_source();
  headerSource->set_filename(AbsPath(fasta));
  headerSource->set_filetype("fasta");
  unsigned int invalidPepCnt = 0;
  unsigned int failedDecoyCnt = 0;

  outPeptideHeap.clear();
  outProteinSequences.clear();

  HeadedRecordWriter proteinWriter(proteinPbFile, outProteinPbHeader);
  ifstream fastaStream(fasta.c_str(), ifstream::in);
  pb::Protein pbProtein;
  string proteinName;
  string* proteinSequence = new string;
  int curProtein = -1;
  vector< pair< ProteinInfo, vector<PeptideInfo> > > cleavedPeptideInfo;
  set<string> setTargets, setDecoys;
  map<const string*, TargetInfo> targetInfo;

  // Iterate over all proteins in FASTA file
  unsigned int targetsGenerated = 0, decoysGenerated = 0;
  while (GenerateDecoys::getNextProtein(fastaStream, proteinName, *proteinSequence)) {
    outProteinSequences.push_back(proteinSequence);
    // Trim protein name, only take first word
    size_t endFirstWord = proteinName.find_first_of(" \t\v\r\n");
    if (endFirstWord != string::npos) {
      proteinName.erase(endFirstWord);
    }
    cleavedPeptideInfo.push_back(make_pair(
      ProteinInfo(proteinName, proteinSequence), vector<PeptideInfo>()));
    const ProteinInfo& proteinInfo = cleavedPeptideInfo.back().first;
    vector<PeptideInfo>& cleavedPeptides = cleavedPeptideInfo.back().second;
    // Write pb::Protein
    getPbProtein(++curProtein, proteinName, *proteinSequence, pbProtein);
    proteinWriter.Write(&pbProtein);
    GenerateDecoys::cleaveProtein(*proteinSequence, enzyme, digestion,
      missedCleavages, minLength, maxLength, cleavedPeptides);
    // Iterate over all generated peptides for this protein
    for (vector<PeptideInfo>::iterator i = cleavedPeptides.begin();
         i != cleavedPeptides.end(); ) {
      const string& cleavedSequence = i->first;
      FLOAT_T pepMass = calcPepMassTide(cleavedSequence, massType);
      if (pepMass < 0.0) {
        // Sequence contained some invalid character
        carp(CARP_DEBUG, "Ignoring invalid sequence <%s>",
             cleavedSequence.c_str());
        ++invalidPepCnt;
        i = cleavedPeptides.erase(i);
        continue;
      } else if (pepMass < minMass || pepMass > maxMass) {
        // Skip to next peptide if not in mass range
        ++i;
        continue;
      }
      const int startLoc = i->second;
      const int pepLen = cleavedSequence.length();
      // Add target to heap
      TideIndexPeptide pepTarget(
        pepMass, pepLen, proteinSequence, curProtein, startLoc, false);
      outPeptideHeap.push_back(pepTarget);
      push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
        greater<TideIndexPeptide>());
      if (decoyType != NO_DECOYS) {
        const string* setTarget = &*(setTargets.insert(cleavedSequence).first);
        targetInfo.insert(make_pair(
          setTarget, TargetInfo(proteinInfo, startLoc, pepMass))).first->second;
      }
      ++targetsGenerated;
      ++i;
    }
    proteinSequence = new string;
  }
  delete proteinSequence;

  // Generate decoys
  map<const string*, const string*> targetToDecoy;
  if (decoyType == PROTEIN_REVERSE_DECOYS) {
    if (decoyFasta) {
      carp(CARP_INFO, "Writing reverse-protein fasta and decoys...");
    }
    vector<PeptideInfo> cleavedReverse;
    for (vector< pair< ProteinInfo, vector<PeptideInfo> > >::const_iterator i =
         cleavedPeptideInfo.begin(); i != cleavedPeptideInfo.end(); ++i) {
      string decoyProtein = *(i->first.sequence);
      reverse(decoyProtein.begin(), decoyProtein.end());
      if (decoyFasta) {
        (*decoyFasta) << ">"<< decoyPrefix << i->first.name << endl
                      << decoyProtein << endl;
      }
      GenerateDecoys::cleaveProtein(decoyProtein, enzyme, digestion,
        missedCleavages, minLength, maxLength, cleavedReverse);
      // Iterate over all generated peptides for this protein
      for (vector<PeptideInfo>::iterator j = cleavedReverse.begin();
           j != cleavedReverse.end();
           ++j) {
        const string& cleavedSequence = j->first;
        FLOAT_T pepMass = calcPepMassTide(cleavedSequence, massType);
        if (pepMass < 0.0) {
          // Sequence contained some invalid character
          carp(CARP_DEBUG, "Ignoring invalid sequence in decoy fasta <%s>",
               cleavedSequence.c_str());
          ++invalidPepCnt;
          continue;
        } else if (pepMass < minMass || pepMass > maxMass) {
          // Skip to next peptide if not in mass range
          continue;
        } else if (setTargets.find(cleavedSequence) != setTargets.end()) {
          // Sequence already exists as a target
          continue;
        }
        const int startLoc = j->second;
        const int pepLen = cleavedSequence.length();
        string* decoySequence = new string(j->first);
        outProteinSequences.push_back(decoySequence);

        // Write pb::Protein
        getDecoyPbProtein(++curProtein, ProteinInfo(i->first.name, &decoyProtein),
                          *decoySequence, startLoc, pbProtein);
        proteinWriter.Write(&pbProtein);
        // Add decoy to heap
        TideIndexPeptide pepDecoy(
          pepMass, cleavedSequence.length(), decoySequence,
          curProtein, (startLoc > 0) ? 1 : 0, true);
        outPeptideHeap.push_back(pepDecoy);
        push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
          greater<TideIndexPeptide>());
        ++decoysGenerated;
      }
    }
  } else {
    for (set<string>::const_iterator i = setTargets.begin();
         i != setTargets.end();
         ++i) {
      const string* setTarget = &*i;
      const map<const string*, TargetInfo>::iterator targetLookup =
        targetInfo.find(setTarget);
      const ProteinInfo& proteinInfo = targetLookup->second.proteinInfo;
      const int startLoc = targetLookup->second.start;
      const map<const string*, const string*>::const_iterator decoyCheck =
        targetToDecoy.find(setTarget);
      string* decoySequence = new string;
      if (decoyCheck != targetToDecoy.end()) {
        // Decoy already generated for this sequence
        *decoySequence = *(decoyCheck->second);
      } else {
        // Try to generate decoy
        if (!GenerateDecoys::makeDecoy(*i, setTargets, setDecoys,
                                       decoyType == PEPTIDE_SHUFFLE_DECOYS,
                                       *decoySequence)) {
        carp(CARP_DEBUG, "Failed to generate decoy for sequence %s",
             i->c_str());
        ++failedDecoyCnt;
        delete decoySequence;
        continue;
        } else {
         targetToDecoy[setTarget] = &*(setDecoys.insert(*decoySequence).first);
        }
      }
      outProteinSequences.push_back(decoySequence);

      // Write pb::Protein
      getDecoyPbProtein(++curProtein, proteinInfo, *decoySequence,
                        startLoc, pbProtein);
      proteinWriter.Write(&pbProtein);
      // Add decoy to heap
      FLOAT_T pepMass = targetLookup->second.mass;
      TideIndexPeptide pepDecoy(
        pepMass, i->length(), decoySequence, curProtein, (startLoc > 0) ? 1 : 0, true);
      outPeptideHeap.push_back(pepDecoy);
      push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
        greater<TideIndexPeptide>());
      ++decoysGenerated;
    }
  }
  if (invalidPepCnt > 0) {
    carp(CARP_INFO, "Ignoring %d peptide sequences containing unrecognized characters", invalidPepCnt);
  }
  if (failedDecoyCnt > 0) {
    carp(CARP_INFO, "Failed to generate decoys for %d low complexity peptides", failedDecoyCnt);
  }
  carp(CARP_DEBUG, "FASTA produced %d targets and %d decoys",
       targetsGenerated, decoysGenerated);

  // Write to decoy fasta if necessary (if protein-reverse, we already wrote it)
  if (decoyFasta && decoyType != PROTEIN_REVERSE_DECOYS) {
    carp(CARP_INFO, "Writing decoy fasta...");
    // Iterate over all (protein, peptides from that protein)
    for (vector< pair< ProteinInfo, vector<PeptideInfo> > >::const_iterator i =
         cleavedPeptideInfo.begin(); i != cleavedPeptideInfo.end(); ++i) {
      string decoyProtein = *(i->first.sequence);
      // Iterate over all peptides from the protein
      for (vector<PeptideInfo>::const_iterator j = i->second.begin();
           j != i->second.end();
           ++j) {
        // In the protein sequence, replace the target peptide with its decoy
        const string& targetSequence = j->first;
        const string* setTarget = &*(setTargets.find(targetSequence));
        const map<const string*, const string*>::const_iterator decoyCheck =
          targetToDecoy.find(setTarget);
        if (decoyCheck != targetToDecoy.end()) {
          decoyProtein.replace(j->second, targetSequence.length(),
                               *(decoyCheck->second));
        }
      }
      // Write out the final protein
      (*decoyFasta) << ">" << decoyPrefix << i->first.name << endl
                    << decoyProtein << endl;
    }
  }
}

void TideIndexApplication::writePeptidesAndAuxLocs(
  vector<TideIndexPeptide>& peptideHeap,
  const string& peptidePbFile,
  const string& auxLocsPbFile,
  pb::Header& pbHeader
) {
  // Check header
  if (pbHeader.source_size() != 1) {
    carp(CARP_FATAL, "pbHeader had a number of sources other than 1");
  }
  pb::Header_Source& headerSource = *(pbHeader.mutable_source(0));
  if (!headerSource.has_filename() || headerSource.has_filetype()) {
    carp(CARP_FATAL, "pbHeader source invalid");
  }

  string proteinsFile = headerSource.filename();
  vector<const pb::Protein*> proteins;
  pb::Header proteinsHeader;
  carp(CARP_INFO, "Reading proteins");
  if (!ReadRecordsToVector<pb::Protein>(&proteins, proteinsFile,
                                        &proteinsHeader)) {
    carp(CARP_FATAL, "Error reading proteins from %s", proteinsFile.c_str());
  } else if (!proteinsHeader.file_type() == pb::Header::RAW_PROTEINS) {
    carp(CARP_FATAL, "Proteins file %s had invalid type", proteinsFile.c_str());
  }
  // Clean up
  for (vector<const pb::Protein*>::iterator i = proteins.begin();
       i != proteins.end();
       ++i) {
    delete *i;
  }
  // The raw proteins file is read in. It's a valid source file;
  // remember it as such:
  headerSource.mutable_header()->CopyFrom(proteinsHeader);

  // Now check other desired settings
  if (!pbHeader.has_peptides_header()) {
    carp(CARP_FATAL, "!pbHeader->has_peptideHeapheader()");
  }
  const pb::Header_PeptidesHeader& settings = pbHeader.peptides_header();
  //if (!Peptide::SetMinMaxMassAndLength(settings)) {
  //  carp(CARP_FATAL, "Error setting min/max mass/length");
  if (!settings.has_enzyme() || settings.enzyme().empty()) {
    carp(CARP_FATAL, "Enzyme settings error");
  }

  pbHeader.set_file_type(pb::Header::PEPTIDES);
  pbHeader.mutable_peptides_header()->set_has_peaks(false);
  pbHeader.mutable_peptides_header()->set_decoys(
    get_tide_decoy_type_parameter("decoy-format"));
  HeadedRecordWriter peptideWriter(peptidePbFile, pbHeader); // put header in outfile

  // Create the auxiliary locations header and writer
  pb::Header auxLocsHeader;
  auxLocsHeader.set_file_type(pb::Header::AUX_LOCATIONS);
  pb::Header_Source* auxLocsSource = auxLocsHeader.add_source();
  auxLocsSource->set_filename(peptidePbFile);
  auxLocsSource->mutable_header()->CopyFrom(pbHeader);
  HeadedRecordWriter auxLocWriter(auxLocsPbFile, auxLocsHeader);

  pb::Peptide pbPeptide;
  pb::AuxLocation pbAuxLoc;
  int auxLocIdx = -1;
  carp(CARP_DEBUG, "%d peptides in heap", peptideHeap.size());
  int count = 0;
  sort_heap(peptideHeap.begin(), peptideHeap.end(),
            greater<TideIndexPeptide>());
  while (!peptideHeap.empty()) {
    TideIndexPeptide curPeptide(peptideHeap.back());
    peptideHeap.pop_back();
    // For duplicate peptides we only record the location
    while (!peptideHeap.empty() && peptideHeap.back() == curPeptide) {
      pb::Location* location = pbAuxLoc.add_location();
      location->set_protein_id(peptideHeap.back().getProteinId());
      location->set_pos(peptideHeap.back().getProteinPos());
      peptideHeap.pop_back();
    }
    getPbPeptide(count, curPeptide, pbPeptide);

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

    if (++count % 100000 == 0) {
      carp(CARP_INFO, "Wrote %d peptides", count);
    }
  }
}

FLOAT_T TideIndexApplication::calcPepMassTide(
  const string& sequence,
  MASS_TYPE_T massType
) {
  FixPt mass;
  FixPt aaMass;
  if (massType == AVERAGE) {
    mass = MassConstants::fixp_avg_h2o;
    for (size_t i = 0; i < sequence.length(); ++i) {
      aaMass = MassConstants::fixp_avg_table[sequence[i]];
      if (aaMass == 0) {
        return -1;
      }
      mass += aaMass;
    }
  } else if (massType == MONO) {
    mass = MassConstants::fixp_mono_h2o;
    for (size_t i = 0; i < sequence.length(); ++i) {
      aaMass = MassConstants::fixp_mono_table[sequence[i]];
      if (aaMass == 0) {
        return -1;
      }
      mass += aaMass;
    }
  } else {
    carp(CARP_FATAL, "Invalid mass type");
  }
  return MassConstants::ToDouble(mass);
}

void TideIndexApplication::getPbProtein(
  int id,
  const string& name,
  const string& residues,
  pb::Protein& outPbProtein
) {
  outPbProtein.Clear();
  outPbProtein.set_id(id);
  outPbProtein.set_name(name);
  outPbProtein.set_residues(residues);
}

void TideIndexApplication::getDecoyPbProtein(
  int id,
  const ProteinInfo& targetProteinInfo,
  string decoyPeptideSequence,
  int startLoc,
  pb::Protein& outPbProtein
) {
  const string* proteinSequence = targetProteinInfo.sequence;
  const int pepLen = decoyPeptideSequence.length();

  // Add N term to decoySequence, if it exists
  if (startLoc > 0) {
    decoyPeptideSequence.insert(0, 1, proteinSequence->at(startLoc - 1));
  }
  // Add C term to decoySequence, if it exists
  size_t cTermLoc = startLoc + pepLen;
  decoyPeptideSequence.push_back((cTermLoc < proteinSequence->length()) ?
    proteinSequence->at(cTermLoc) : '-');
  // Append original target sequence, unless using protein level decoys
  if (get_tide_decoy_type_parameter("decoy-format") != PROTEIN_REVERSE_DECOYS) {
    decoyPeptideSequence.append(targetProteinInfo.sequence->substr(startLoc, pepLen));
  }

  getPbProtein(id, get_string_parameter_pointer("decoy-prefix") + targetProteinInfo.name,
               decoyPeptideSequence, outPbProtein);
  outPbProtein.set_target_pos(startLoc);
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
  outPbPeptide.mutable_first_location()->set_pos(peptide.getProteinPos());;
  outPbPeptide.set_is_decoy(peptide.isDecoy());
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

void TideIndexApplication::writeParamFile() {
  char* param_file_name = cat_string(getFileStem().c_str(), ".params.txt");

  // Update mods-spec parameter for default cysteine mod
  const string default_cysteine = "C+57.02146";
  string mods_spec = get_string_parameter_pointer("mods-spec");
  if (mods_spec.find('C') == string::npos) {
    mods_spec = mods_spec.empty() ?
      default_cysteine : default_cysteine + ',' + mods_spec;
    carp(CARP_DEBUG, "Using default cysteine mod '%s' ('%s')",
         default_cysteine.c_str(), mods_spec.c_str());
  }
  add_or_update_hash(parameters, "mods-spec", mods_spec.c_str());

  print_parameter_file(&param_file_name);
  free(param_file_name);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
