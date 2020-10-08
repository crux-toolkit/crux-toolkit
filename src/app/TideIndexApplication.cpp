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

#ifdef _MSC_VER
#include <io.h>
#endif

extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
                                const string& input_filename,
                                const string& output_filename);
extern void AddMods(HeadedRecordReader* reader,
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
    var_mod_table.ParsedCtproModTable(), 0, 0)) {
    carp(CARP_FATAL, "Error in MassConstants::Init");
  }

  DECOY_TYPE_T decoy_type = get_tide_decoy_type_parameter("decoy-format");
  string decoyPrefix = Params::GetString("decoy-prefix");

  // Set up output paths
  bool overwrite = Params::GetBool("overwrite");

  if (!FileUtils::Exists(fasta)) {
    carp(CARP_FATAL, "Fasta file %s does not exist", fasta.c_str());
  }
  else{
    carp(CARP_INFO, "Fasta file %s is found.", fasta.c_str());
  }

  string out_proteins = FileUtils::Join(index, "protix");
  string out_peptides = FileUtils::Join(index, "pepix");
  string out_aux = FileUtils::Join(index, "auxlocs");
  string modless_peptides = out_peptides + ".nomods.tmp";
  string peakless_peptides = out_peptides + ".nopeaks.tmp";
  ofstream* out_target_list = NULL;
  ofstream* out_decoy_list = NULL;
  if (Params::GetBool("peptide-list")) {
    out_target_list = create_stream_in_path(make_file_path(
      "tide-index.peptides.target.txt").c_str(), NULL, overwrite);
    if (decoy_type != NO_DECOYS) {
      out_decoy_list = create_stream_in_path(make_file_path(
        "tide-index.peptides.decoy.txt").c_str(), NULL, overwrite);
    }
  }
  ofstream* out_decoy_fasta = GeneratePeptides::canGenerateDecoyProteins() ?
    create_stream_in_path(make_file_path(
      "tide-index.decoy.fasta").c_str(), NULL, overwrite) : NULL;

  if (create_output_directory(index.c_str(), overwrite) != 0) {
    carp(CARP_FATAL, "Error creating index directory");
  } else if (FileUtils::Exists(out_proteins) ||
             FileUtils::Exists(out_peptides) ||
             FileUtils::Exists(out_aux)) {
    if (overwrite) {
      carp(CARP_DEBUG, "Cleaning old index file(s)");
      FileUtils::Remove(out_proteins);
      FileUtils::Remove(out_peptides);
      FileUtils::Remove(out_aux);
      FileUtils::Remove(modless_peptides);
      FileUtils::Remove(peakless_peptides);
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
            min_length, max_length, allowDups, mass_type, decoy_type, fasta, out_proteins,
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
  if (enzyme != "no-enzyme") {
    pep_header.set_full_digestion(digestion == FULL_DIGEST);
    pep_header.set_max_missed_cleavages(missed_cleavages);
  }
  pep_header.mutable_mods()->CopyFrom(*(var_mod_table.ParsedModTable()));
  pep_header.mutable_nterm_mods()->CopyFrom(*(var_mod_table.ParsedNtpepModTable()));
  pep_header.mutable_cterm_mods()->CopyFrom(*(var_mod_table.ParsedCtpepModTable()));
  pep_header.mutable_nprotterm_mods()->CopyFrom(*(var_mod_table.ParsedNtproModTable()));
  pep_header.mutable_cprotterm_mods()->CopyFrom(*(var_mod_table.ParsedCtproModTable()));
  int numDecoys;
  switch (decoy_type) {
    case NO_DECOYS:
      numDecoys = 0;
      break;
    default:
      numDecoys = 1;
      break;
    case PEPTIDE_SHUFFLE_DECOYS:
      numDecoys = Params::GetInt("num-decoys-per-target");
      break;
  }
  switch (numDecoys) {
    case 0:
      carp(CARP_INFO, "No decoys will be generated");
      break;
    case 1:
      carp(CARP_INFO, "Generating 1 decoy per target");
      break;
    default:
      carp(CARP_INFO, "Generating %d decoys per target", numDecoys);
  }
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

  string basic_peptides = need_mods ? modless_peptides : peakless_peptides;

  writePeptidesAndAuxLocs(peptideHeap, basic_peptides, out_aux, header_no_mods);
  // Do some clean up
  for (vector<string*>::iterator i = proteinSequences.begin();
       i != proteinSequences.end();
       ++i) {
    delete *i;
  }
  vector<TideIndexPeptide>().swap(peptideHeap);
  ProteinVec proteins;
  if (!ReadRecordsToVector<pb::Protein>(&proteins, out_proteins)) {
    carp(CARP_FATAL, "Error reading proteins file");
  }

  if (need_mods) {
    carp(CARP_INFO, "Computing modified peptides...");
    HeadedRecordReader reader(modless_peptides, NULL, 1024 << 10); // 1024kb buffer
    AddMods(&reader, peakless_peptides, Params::GetString("temp-dir"), header_with_mods, proteins, &var_mod_table);
  }

  if (out_target_list) {
    // Write peptide lists
    carp(CARP_INFO, "Writing peptide lists...");

    // This set holds target peptide strings
    set<string> targetPepStrs;
    // This vector holds decoy peptide strings and masses
    vector< pair<string, double> > decoyPepStrs;
    // Read peptides protocol buffer file
    vector<const pb::AuxLocation*> locations;
    if (!ReadRecordsToVector<pb::AuxLocation>(&locations, out_aux)) {
      carp(CARP_FATAL, "Error reading auxlocs file");
    }
    int mass_precision = Params::GetInt("mass-precision");
    // Iterate over all protocol buffer peptides
    unsigned int writeCountTargets = 0, writeCountDecoys = 0;
    HeadedRecordReader reader(peakless_peptides, NULL);
    while (!reader.Done()) {
      pb::Peptide* protobuf = new pb::Peptide;
      reader.Read(protobuf);
      pb::Peptide* peptide = protobuf;
      bool writeTarget = true;
      bool writeDecoy = false;

      if (out_decoy_list) {
        writeDecoy = peptide->has_decoy_index();
        writeTarget = !writeDecoy;
      }
      string pep_str = getModifiedPeptideSeq(peptide, &proteins);

      if (writeTarget) {
        // This is a target, output it
        targetPepStrs.insert(pep_str);
        *out_target_list << pep_str << '\t'
                         << StringUtils::ToString(peptide->mass(),
                                                  mass_precision) << endl;
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
    for (vector< pair<string, double> >::const_iterator i = decoyPepStrs.begin();
         i != decoyPepStrs.end();
         ++i) {
      *out_decoy_list << i->first << '\t'
                      << StringUtils::ToString(i->second, mass_precision);
      if (targetPepStrs.find(i->first) != targetPepStrs.end()) {
        *out_decoy_list << "\t*";
      }
      *out_decoy_list << endl;
    }

    // Close and clean up streams
    if (out_decoy_list) {
      out_decoy_list->close();
      delete out_decoy_list;
    }
    out_target_list->close();
    delete out_target_list;

    carp(CARP_DETAILED_INFO, "Wrote %d targets and %d decoys to peptide list",
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

  // clean up out_decoy_fasta
  if (out_decoy_fasta) {
    delete out_decoy_fasta;
  }

  // Recover stderr
  cerr.rdbuf(old);
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
    "temp-dir",
    "verbosity",
    "aws-profile"
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

void TideIndexApplication::fastaToPb(
  const string& commandLine,
  const ENZYME_T enzyme,
  const DIGEST_T digestion,
  int missedCleavages,
  FLOAT_T minMass,
  FLOAT_T maxMass,
  int minLength,
  int maxLength,
  bool allowDups,
  MASS_TYPE_T massType,
  DECOY_TYPE_T decoyType,
  const string& fasta,
  const string& proteinPbFile,
  pb::Header& outProteinPbHeader,
  vector<TideIndexPeptide>& outPeptideHeap,
  vector<string*>& outProteinSequences,
  ofstream* decoyFasta
) {
  typedef GeneratePeptides::CleavedPeptide PeptideInfo;

  string decoyPrefix = Params::GetString("decoy-prefix");
  outProteinPbHeader.Clear();
  outProteinPbHeader.set_file_type(pb::Header::RAW_PROTEINS);
  outProteinPbHeader.set_command_line(commandLine);
  pb::Header_Source* headerSource = outProteinPbHeader.add_source();
  headerSource->set_filename(FileUtils::AbsPath(fasta));
  headerSource->set_filetype("fasta");
  unsigned int invalidPepCnt = 0;
  unsigned int failedDecoyCnt = 0;

  outPeptideHeap.clear();
  outProteinSequences.clear();

  HeadedRecordWriter proteinWriter(proteinPbFile, outProteinPbHeader);

  istream& fastaStream = FileUtils::GetReadStream(fasta); 
  string proteinName;
  string* proteinSequence = new string;
  int curProtein = -1;
  vector< pair< ProteinInfo, vector<PeptideInfo> > > cleavedPeptideInfo;
  set<string> setTargets, setDecoys;
  map<const string*, TargetInfo> targetInfo;

  // Iterate over all proteins in FASTA file
  unsigned int targetsGenerated = 0, decoysGenerated = 0;
  while (GeneratePeptides::getNextProtein(fastaStream, &proteinName, proteinSequence)) {
    outProteinSequences.push_back(proteinSequence);
    cleavedPeptideInfo.push_back(make_pair(
      ProteinInfo(proteinName, proteinSequence), vector<PeptideInfo>()));
    const ProteinInfo& proteinInfo = cleavedPeptideInfo.back().first;
    vector<PeptideInfo>& cleavedPeptides = cleavedPeptideInfo.back().second;
    // Write pb::Protein
    writePbProtein(proteinWriter, ++curProtein, proteinName, *proteinSequence);
    cleavedPeptides = GeneratePeptides::cleaveProtein(
      *proteinSequence, enzyme, digestion, missedCleavages, minLength, maxLength);
    // Iterate over all generated peptides for this protein
    for (vector<PeptideInfo>::iterator i = cleavedPeptides.begin();
         i != cleavedPeptides.end(); ) {
      FLOAT_T pepMass = calcPepMassTide(&(*i), massType, &proteinInfo);
      if (pepMass < 0.0) {
        // Sequence contained some invalid character
        carp(CARP_DEBUG, "Ignoring invalid sequence <%s>", i->Sequence().c_str());
        ++invalidPepCnt;
        i = cleavedPeptides.erase(i);
        continue;
      } else if (pepMass < minMass || pepMass > maxMass) {
        // Skip to next peptide if not in mass range
        ++i;
        continue;
      }
      // Add target to heap
      TideIndexPeptide pepTarget(pepMass, i->Length(), proteinSequence, curProtein, i->Position());
      outPeptideHeap.push_back(pepTarget);
      push_heap(outPeptideHeap.begin(), outPeptideHeap.end(), greater<TideIndexPeptide>());
      if (!allowDups && decoyType != NO_DECOYS) {
        const string* setTarget = &*(setTargets.insert(i->Sequence()).first);
        targetInfo.insert(make_pair(setTarget, TargetInfo(proteinInfo, i->Position(), pepMass)));
      }
      ++targetsGenerated;
      ++i;
    }
    proteinSequence = new string;
  }
  delete proteinSequence;
  FileUtils::CloseStream(fastaStream);
  if (targetsGenerated == 0) {
    carp(CARP_FATAL, "No target sequences generated.  Is \'%s\' a FASTA file?",
         fasta.c_str());
  }
  if (invalidPepCnt > 0) {
    carp(CARP_INFO, "Ignoring %d peptide sequences containing unrecognized characters.", invalidPepCnt);
  }
  carp(CARP_INFO, "Generated %d targets, including duplicates.", targetsGenerated);

  // Generate decoys
  map< const string, vector<const string*> > targetToDecoy;
  int numDecoys = Params::GetInt("num-decoys-per-target");
  if (decoyType == PROTEIN_REVERSE_DECOYS) {
    if (decoyFasta) {
      carp(CARP_INFO, "Writing reverse-protein fasta and decoys...");
    }
    for (vector< pair< ProteinInfo, vector<PeptideInfo> > >::const_iterator i =
         cleavedPeptideInfo.begin(); i != cleavedPeptideInfo.end(); ++i) {
      string decoyProtein = *(i->first.sequence);
      reverse(decoyProtein.begin(), decoyProtein.end());
      if (decoyFasta) {
        (*decoyFasta) << ">"<< decoyPrefix << i->first.name << endl
                      << decoyProtein << endl;
      }
      vector<PeptideInfo> cleavedReverse = GeneratePeptides::cleaveProtein(
        decoyProtein, enzyme, digestion, missedCleavages, minLength, maxLength);
      // Iterate over all generated peptides for this protein
      for (vector<PeptideInfo>::iterator j = cleavedReverse.begin();
           j != cleavedReverse.end();
           ++j) {
        FLOAT_T pepMass = calcPepMassTide(&(*j), massType, &(i->first));
        if (pepMass < 0.0) {
          // Sequence contained some invalid character
          carp(CARP_DEBUG, "Ignoring invalid sequence in decoy fasta <%s>",
               j->Sequence().c_str());
          ++invalidPepCnt;
          continue;
        } else if (pepMass < minMass || pepMass > maxMass) {
          // Skip to next peptide if not in mass range
          continue;
        } else if (!allowDups && setTargets.find(j->Sequence()) != setTargets.end()) {
          // Sequence already exists as a target
          continue;
        }
        string* decoySequence = new string(j->Sequence());
        outProteinSequences.push_back(decoySequence);

        // Write pb::Protein
        writeDecoyPbProtein(++curProtein, ProteinInfo(i->first.name, &decoyProtein),
                            *decoySequence, j->Position(), proteinWriter);
        // Add decoy to heap
        TideIndexPeptide pepDecoy(pepMass, j->Length(), decoySequence,
          curProtein, (j->Position() > 0) ? 1 : 0, 0);
        outPeptideHeap.push_back(pepDecoy);
        push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
          greater<TideIndexPeptide>());
        ++decoysGenerated;
      }
    }
  } else if (!allowDups) {
    for (set<string>::const_iterator i = setTargets.begin();
         i != setTargets.end();
         ++i) {
      const string* setTarget = &*i;
      const map<const string*, TargetInfo>::iterator targetLookup =
        targetInfo.find(setTarget);
      const ProteinInfo& proteinInfo = (targetLookup->second.proteinInfo);
      const int startLoc = targetLookup->second.start;
      FLOAT_T pepMass = targetLookup->second.mass;
      generateDecoys(numDecoys, *setTarget, targetToDecoy, &setTargets, &setDecoys, decoyType, allowDups,
                     failedDecoyCnt, decoysGenerated, curProtein, proteinInfo, startLoc, proteinWriter,
                     pepMass, outPeptideHeap, outProteinSequences);
    }
  } else { // allow dups
    for (vector<pair<ProteinInfo, vector<PeptideInfo> > >::const_iterator i = cleavedPeptideInfo.begin();
         i != cleavedPeptideInfo.end();
         ++i) {
      const ProteinInfo& proteinInfo = i->first;
      for (vector<PeptideInfo>::const_iterator j = i->second.begin();
           j != i->second.end();
           ++j) {
        const string setTarget = j->Sequence();
        const int startLoc = j->Position();
        FLOAT_T pepMass = calcPepMassTide(&(*j), massType, &proteinInfo);
        generateDecoys(numDecoys, setTarget, targetToDecoy, NULL, NULL, decoyType, allowDups, failedDecoyCnt,
                       decoysGenerated, curProtein, proteinInfo, startLoc, proteinWriter,
                       pepMass, outPeptideHeap, outProteinSequences);
      }
    }
  }
  if (failedDecoyCnt > 0) {
    carp(CARP_INFO, "Failed to generate decoys for %d low complexity peptides.", failedDecoyCnt);
  }
  carp(CARP_INFO, "Generated %d decoys.", decoysGenerated);

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
        const string setTarget = j->Sequence();
        const map< const string, vector<const string*> >::const_iterator decoyCheck = targetToDecoy.find(setTarget);
        if (decoyCheck != targetToDecoy.end() && !decoyCheck->second.empty()) {
          decoyProtein.replace(j->Position(), j->Length(), *(decoyCheck->second.front()));
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
  } else if (proteinsHeader.file_type() != pb::Header::RAW_PROTEINS) {
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
  carp(CARP_DETAILED_INFO, "%d peptides in heap", peptideHeap.size());
  int count = 0;
  int numTargets = 0;
  int numDecoys = 0;
  int numDuplicateTargets = 0;
  int numDuplicateDecoys = 0;
  sort_heap(peptideHeap.begin(), peptideHeap.end(),
            greater<TideIndexPeptide>());
  while (!peptideHeap.empty()) {
    TideIndexPeptide curPeptide(peptideHeap.back());
    peptideHeap.pop_back();
    // For duplicate peptides we only record the location
    while (!peptideHeap.empty() && peptideHeap.back() == curPeptide) {
      if (peptideHeap.back().isDecoy()) {
        numDuplicateDecoys++;
      } else {
        numDuplicateTargets++;
      }        
      carp(CARP_DEBUG, "Skipping duplicate %s.", curPeptide.getSequence().c_str());
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

    if (curPeptide.isDecoy()) {
      numDecoys++;
    } else {
      numTargets++;
    }
    if (++count % 100000 == 0) {
      carp(CARP_INFO, "Wrote %d peptides", count);
    }
  }
  carp(CARP_INFO, "Skipped %d duplicate targets and %d duplicate decoys.",
       numDuplicateTargets, numDuplicateDecoys);
  carp(CARP_INFO, "Wrote %d targets and %d decoys.", numTargets, numDecoys);
}

FLOAT_T TideIndexApplication::calcPepMassTide(
  const GeneratePeptides::CleavedPeptide* pep,
  MASS_TYPE_T massType,
  const TideIndexApplication::ProteinInfo* prot
) {
  FixPt mass;
  FixPt aaMass;
  const string sequence = pep->Sequence();
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

  for (size_t i = 0; i < sequence.length(); ++i) {
    if (i == 0) {
      if(pep->Position() == 0)  //apply protein terminal mod if this is protein N-terminal
        aaMass = _tables->nprotterm_table[sequence[i]];
      else
        aaMass = _tables->nterm_table[sequence[i]];
    } else if (i == sequence.length() - 1) {
      if((pep->Position() + pep->Length()) == prot->sequence->length())  //check if this is protein C-terminal
        aaMass = _tables->cprotterm_table[sequence[i]];
      else
        aaMass = _tables->cterm_table[sequence[i]];
    } else {
      aaMass = _tables->_table[sequence[i]];
    }
    if (aaMass == 0) {
      return -1;
    }
    mass += aaMass;
  }
  return MassConstants::ToDouble(mass);
}

void TideIndexApplication::writePbProtein(
  HeadedRecordWriter& writer,
  int id,
  const string& name,
  const string& residues,
  int targetPos
) {
  static pb::Protein p;
  p.Clear();
  p.set_id(id);
  p.set_name(name);
  p.set_residues(residues);
  if (targetPos >= 0) {
    p.set_target_pos(targetPos);
  }
  writer.Write(&p);
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
void TideIndexApplication::writeDecoyPbProtein(
  int id,
  const ProteinInfo& targetProteinInfo,
  string decoyPeptideSequence,
  int startLoc,
  HeadedRecordWriter& proteinWriter
) {
  const string* proteinSequence = targetProteinInfo.sequence;
  const int pepLen = decoyPeptideSequence.length();

  // Add N term to decoySequence, if it exists
  if (startLoc > 0) {
    decoyPeptideSequence.insert(0, 1, proteinSequence->at(startLoc - 1));
  }
  // Add C term to decoySequence, if it exists, or hyphen otherwise.
  size_t cTermLoc = startLoc + pepLen;
  decoyPeptideSequence.push_back((cTermLoc < proteinSequence->length()) ?
    proteinSequence->at(cTermLoc) : '-');

  // Append original target sequence, unless using protein level decoys
  if (get_tide_decoy_type_parameter("decoy-format") != PROTEIN_REVERSE_DECOYS) {
    decoyPeptideSequence.append(targetProteinInfo.sequence->substr(startLoc, pepLen));
  }
  writePbProtein(proteinWriter, id, Params::GetString("decoy-prefix") + targetProteinInfo.name,
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

void TideIndexApplication::generateDecoys(
  int numDecoys,
  const string& setTarget,
  std::map< const string, vector<const string*> >& targetToDecoy,
  set<string>* setTargets,
  set<string>* setDecoys,
  DECOY_TYPE_T decoyType,
  bool allowDups,
  unsigned int& failedDecoyCnt,
  unsigned int& decoysGenerated,
  int& curProtein,
  const ProteinInfo& proteinInfo,
  const int startLoc,
  HeadedRecordWriter& proteinWriter,
  FLOAT_T pepMass,
  vector<TideIndexPeptide>& outPeptideHeap,
  vector<string*>& outProteinSequences
) {
  vector<string*> decoySequences;
  int generateAttemptsMax = 6;
  const map< const string, vector<const string*> >::const_iterator decoyCheck = targetToDecoy.find(setTarget);
  if (decoyCheck != targetToDecoy.end()) {
    // Decoys already generated for this sequence
    decoySequences = vector<string*>(decoyCheck->second.size(), NULL);
    for (size_t i = 0; i < decoyCheck->second.size(); i++) {
      decoySequences[i] = new string(*(decoyCheck->second[i]));
    }
  } else {
    // Try to generate decoys
    bool shuffle = decoyType == PEPTIDE_SHUFFLE_DECOYS;
    if (!shuffle) {
      numDecoys = 1;
      generateAttemptsMax = 1;
    }
    set<string> targets, decoys;
    if (allowDups) {
      setTargets = &targets;
      setDecoys = &decoys;
    }
    set<string> generatedDecoys;
    set<string> dummy;
    for (int i = 0; i < numDecoys; i++) {
      string* outSeq = new string;
      bool success = false;
      for (int j = 0; j < generateAttemptsMax; j++) {
        success = GeneratePeptides::makeDecoy(setTarget, *setTargets, dummy, shuffle, *outSeq);
        if (success) {
          break;
        }
      }
      if (!success) {
        carp(CARP_DEBUG, "Failed to generate decoys for sequence %s", setTarget.c_str());
        delete outSeq;
        ++failedDecoyCnt;
        return;
      }
      decoySequences.push_back(outSeq);
    }
    map< const string, vector<const string*> >::iterator i =
      targetToDecoy.insert(make_pair(setTarget, vector<const string*>())).first;
    if (!allowDups) {
      for (vector<string*>::const_iterator j = decoySequences.begin(); j != decoySequences.end(); j++) {
        set<string>::iterator k = setDecoys->insert(**j).first;
        i->second.push_back(&*k);
      }
    } else {
      for (vector<string*>::const_iterator j = decoySequences.begin(); j != decoySequences.end(); j++) {
        i->second.push_back(*j);
      }
    }
  }

  for (int i = 0; i < numDecoys; i++) {
    string* seq = decoySequences[i];
    carp(CARP_DETAILED_DEBUG, "Got decoy sequence %d: %s.", i, seq->c_str());
    outProteinSequences.push_back(seq);
    // Write pb::Protein In this subroutine, the startLoc is used to
    // construct a longer sequence containing N- and C-term residues,
    // plus the target.
    writeDecoyPbProtein(++curProtein, proteinInfo, *seq, startLoc, proteinWriter);
    // Add decoy to heap
    TideIndexPeptide pepDecoy(pepMass, setTarget.length(), seq, curProtein, (startLoc > 0) ? 1 : 0, i);
    outPeptideHeap.push_back(pepDecoy);
    push_heap(outPeptideHeap.begin(), outPeptideHeap.end(), greater<TideIndexPeptide>());
  }
  decoysGenerated += decoySequences.size();
 }


/*
* Local Variables:
* mode: c
* c-basic-offset: 2
* End:
*/
