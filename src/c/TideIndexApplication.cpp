#include <cstdio>
#include <fstream>
#include "carp.h"
#include "CarpStreamBuf.h"
#include "GenerateDecoys.h"
#include "TideIndexApplication.h"

#include "tide/modifications.h"
#include "tide/records_to_vector-inl.h"

extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
                        				const string& input_filename,
                        				const string& output_filename);
extern void AddMods(HeadedRecordReader* reader,
                    string out_file,
            		    const pb::Header& header,
            		    const vector<const pb::Protein*>& proteins);

TideIndexApplication::TideIndexApplication() {
}

TideIndexApplication::~TideIndexApplication() {
}

int TideIndexApplication::main(int argc, char** argv) {

  const char* option_list[] = {
    "decoy-format",
    "enzyme",
    "digestion",
    "missed-cleavages",
    "max-length",
    "max-mass",
    "min-length",
    "min-mass",
    "monoisotopic-precursor",
    "mods-spec",
    "output-dir",
    "overwrite",
    "peptide-list",
    "parameter-file",
    "verbosity"
  };

  const string default_cysteine = "C+57.0214637206";

  // Crux command line parsing
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "protein fasta file",
    "index name"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running tide-index...");

  // Build command line string
  string cmd_line = "crux tide-index";
  for (int i = 1; i < argc; ++i) {
    cmd_line += " ";
    cmd_line += argv[i];
  }

  // Get options
  double min_mass = get_double_parameter("min-mass");
  double max_mass = get_double_parameter("max-mass");
  int min_length = get_int_parameter("min-length");
  int max_length = get_int_parameter("max-length");
  bool monoisotopic_precursor = get_boolean_parameter("monoisotopic-precursor");
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
  string mods_spec = get_string_parameter_pointer("mods-spec");
  if (mods_spec.find('C') == string::npos) {
    mods_spec = (mods_spec.empty()) ?
      default_cysteine : default_cysteine + ',' + mods_spec;
    carp(CARP_DEBUG, "Using default cysteine mod '%s' ('%s')",
         default_cysteine.c_str(), mods_spec.c_str());
  }
  VariableModTable var_mod_table;
  if (!var_mod_table.Parse(mods_spec.c_str())) {
    carp(CARP_FATAL, "Error parsing mods");
  }
  if (!MassConstants::Init(var_mod_table.ParsedModTable())) {
    carp(CARP_FATAL, "Error in MassConstants::Init");
  }

  string decoy_format = get_string_parameter_pointer("decoy-format");
  DECOY_TYPE decoy_type;
  if (decoy_format == "shuffle") {
    decoy_type = SHUFFLE;
    carp(CARP_DEBUG, "Using shuffled decoys");
  } else if (decoy_format == "reverse") {
    decoy_type = REVERSE;
    carp(CARP_DEBUG, "Using reversed decoys");
  } else if (decoy_format == "none") {
    decoy_type = NONE;
    carp(CARP_DEBUG, "Not using decoys");
  } else {
    carp(CARP_FATAL, "Invalid decoy type %s", decoy_format.c_str());
  }

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
    if (decoy_type != NONE) {
      out_decoy_list = create_stream_in_path(make_file_path(
        "tide-index.peptides.decoy.txt").c_str(), NULL, overwrite);
    }
  }

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

  // Reroute stderr
  CarpStreamBuf buffer;
  streambuf* old = cerr.rdbuf();
  cerr.rdbuf(&buffer);

  // Start tide-index
  carp(CARP_INFO, "Reading %s and computing unmodified peptides...",
       fasta.c_str());
  pb::Header proteinPbHeader;
  vector<TideIndexPeptide> peptideHeap;
  vector<string*> proteinSequences;
  fastaToPb(cmd_line, enzyme_t, digestion, missed_cleavages, min_mass, max_mass,
            min_length, max_length, mass_type, decoy_type, fasta, out_proteins,
            proteinPbHeader, peptideHeap, proteinSequences);

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

  bool need_mods = header_with_mods.peptides_header().mods().variable_mod_size() > 0;
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
    AddMods(&reader, peakless_peptides, header_with_mods, proteins);
  }

  if (out_target_list) {
    // Write peptide lists
    carp(CARP_INFO, "Writing peptide lists...");
    // Initialize modification variables
    stringstream mod_stream;
    mod_stream.precision(get_int_parameter("mod-precision"));
    int mod_index;
    double mod_delta;
    // Initialize decoy variables
    string decoyPrefix = (decoy_type == SHUFFLE) ? "rand_" : "rev_";
    size_t decoyPrefixLen = decoyPrefix.size();
    // This set holds target peptide strings
    set<string> targetPepStrs;
    // This vector holds decoy peptide strings and masses
    vector< pair<string, double> > decoyPepStrs;
    // Read peptides pb file
    vector<const pb::Peptide*> peptides;
    if (!ReadRecordsToVector<pb::Peptide>(&peptides, peakless_peptides)) {
      carp(CARP_FATAL, "Error reading peptides file");
    }
    // Iterate over all pb peptides
    while (!peptides.empty()) {
      const pb::Peptide* peptide = peptides.back();
      const pb::Location& location = peptide->first_location();
      const pb::Protein* protein = proteins[location.protein_id()];
      // Get peptide sequence without mods
      string pep_str = protein->residues().substr(
        location.pos(), peptide->length());
      // Store all mod indices/deltas
      map<int, double> mod_map;
      set<int> mod_indices;
      for (int j = 0; j < peptide->modifications_size(); ++j) {
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
      // Write the modified sequence to the appropriate file stream
      if (!out_decoy_list ||
          (!protein->name().empty() && protein->name()[0] != DecoyMagicByte)) {
        // This is a target, output it
        targetPepStrs.insert(pep_str);
        *out_target_list << pep_str << '\t' << peptide->mass() << endl;
      } else {
        // This is a decoy, save it to output later
        decoyPepStrs.push_back(make_pair(pep_str, peptide->mass()));
      }
      peptides.pop_back();
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
  return "Create an index for all peptides in a fasta file.";
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
  DECOY_TYPE decoyType,
  const string& fasta,
  const string& proteinPbFile,
  pb::Header& outProteinPbHeader,
  vector<TideIndexPeptide>& outPeptideHeap,
  vector<string*>& outProteinSequences
) {
  outProteinPbHeader.Clear();
  outProteinPbHeader.set_file_type(pb::Header::RAW_PROTEINS);
  outProteinPbHeader.set_command_line(commandLine);
  pb::Header_Source* headerSource = outProteinPbHeader.add_source();
  headerSource->set_filename(AbsPath(fasta));
  headerSource->set_filetype("fasta");

  outPeptideHeap.clear();
  outProteinSequences.clear();

  HeadedRecordWriter proteinWriter(proteinPbFile, outProteinPbHeader);
  ifstream fastaStream(fasta.c_str(), ifstream::in);
  pb::Protein pbProtein;
  string proteinName;
  string* proteinSequence = new string;
  int curProtein = -1;
  int curTargetProtein;
  vector< pair<string, int> > cleavedPeptides;
  set<string> dummySet;
  string decoyPrefix;
  if (decoyType == SHUFFLE) {
    decoyPrefix = "rand_";
  } else if (decoyType == REVERSE) {
    decoyPrefix = "rev_";
  }
  // Iterate over all proteins in FASTA file
  while (GenerateDecoys::getNextProtein(fastaStream, proteinName, *proteinSequence)) {
    curTargetProtein = ++curProtein;
    outProteinSequences.push_back(proteinSequence);
    // Trim protein name, only take first word
    size_t endFirstWord = proteinName.find_first_of(" \t\v\r\n");
    if (endFirstWord != string::npos) {
      proteinName.erase(endFirstWord);
    }
    // Write pb::Protein
    getPbProtein(curTargetProtein, proteinName, *proteinSequence, pbProtein);
    proteinWriter.Write(&pbProtein);
    GenerateDecoys::cleaveProtein(*proteinSequence, enzyme, digestion,
      missedCleavages, minLength, maxLength, cleavedPeptides);
    // Iterate over all generated peptides for this protein
    for (vector< pair<string, int> >::iterator i = cleavedPeptides.begin();
         i != cleavedPeptides.end();
         ++i) {
      const string& cleavedSequence = i->first;
      FLOAT_T pepMass = calcPepMassTide(cleavedSequence, massType);
      if (pepMass < 0.0) {
        // Sequence contained some invalid character
        carp(CARP_WARNING, "Ignoring invalid sequence %s",
             cleavedSequence.c_str());
        continue;
      }
      else if (pepMass < minMass || pepMass > maxMass) {
        // Skip to next peptide if not in mass range
        continue;
      }
      const int startLoc = i->second;
      int pepLen = cleavedSequence.length();
      // Add target to heap
      TideIndexPeptide pepTarget(
        pepMass, pepLen, proteinSequence, curTargetProtein, startLoc);
      outPeptideHeap.push_back(pepTarget);
      push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
        greater<TideIndexPeptide>());
      // Skip to next peptide if not generating decoys
      if (decoyType == NONE) {
        continue;
      }
      // Try to generate decoy
      string* decoySequence = new string;
      if (GenerateDecoys::makeDecoy(cleavedSequence, dummySet, dummySet,
                                    decoyType == SHUFFLE, *decoySequence)) {
        // Successfully generated decoy
        outProteinSequences.push_back(decoySequence);
        // Add N term to decoySequence, if it exists
        char nTerm = (startLoc > 0) ? proteinSequence->at(startLoc - 1) : '\0';
        if (nTerm != '\0') {
          decoySequence->insert(0, 1, nTerm);
        }
        // Add C term to decoySequence, if it exists
        size_t cTermLoc = startLoc + pepLen;
        if (cTermLoc < proteinSequence->length()) {
          decoySequence->push_back(proteinSequence->at(cTermLoc));
        }
        // Append unshuffled sequence
        decoySequence->append(cleavedSequence);
        // Write pb::Protein
        stringstream decoyStream;
        decoyStream << DecoyMagicByte << (startLoc + 1) << '.'
                    << decoyPrefix << proteinName;
        getPbProtein(++curProtein, decoyStream.str(), *decoySequence, pbProtein);
        proteinWriter.Write(&pbProtein);
        // Add decoy to heap
        TideIndexPeptide pepDecoy(
          pepMass, pepLen, decoySequence, curProtein, (nTerm != '\0') ? 1 : 0);
        outPeptideHeap.push_back(pepDecoy);
        push_heap(outPeptideHeap.begin(), outPeptideHeap.end(),
          greater<TideIndexPeptide>());
      } else {
        carp(CARP_WARNING, "Failed to generate decoy for sequence %s",
             cleavedSequence.c_str());
        delete decoySequence;
      }
    }
    proteinSequence = new string;
  }
  delete proteinSequence;
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
  HeadedRecordWriter peptideWriter(peptidePbFile, pbHeader); // put header in outfile

  // Create the auxiliary locations header and writer
  pb::Header auxLocsHeader;
  auxLocsHeader.set_file_type(pb::Header::AUX_LOCATIONS);
  pb::Header_Source* auxLocsSource = auxLocsHeader.add_source();
  auxLocsSource->set_filename(peptidePbFile);
  auxLocsSource->mutable_header()->CopyFrom(pbHeader);
  HeadedRecordWriter auxLocWriter(auxLocsPbFile, auxLocsHeader);

  pb::Peptide pbPeptide;
  pb::AuxLocation* pbAuxLoc = new pb::AuxLocation();
  int auxLocIdx = 0;
  carp(CARP_DEBUG, "%d peptides in heap", peptideHeap.size());
  int count = 0;
  sort_heap(peptideHeap.begin(), peptideHeap.end(),
            greater<TideIndexPeptide>());
  while (!peptideHeap.empty()) {
    TideIndexPeptide curPeptide(peptideHeap.back());
    peptideHeap.pop_back();
    // For duplicate peptides we only record the location
    while (!peptideHeap.empty() && peptideHeap.back() == curPeptide) {
      pb::Location* location = pbAuxLoc->add_location();
      location->set_protein_id(peptideHeap.back().getProteinId());
      location->set_pos(peptideHeap.back().getProteinPos());
      peptideHeap.pop_back();
    }
    getPbPeptide(count, curPeptide, pbPeptide);

    // Not all peptides have aux locations associated with them. Check to see
    // if GetGroup added any locations to aux_location. If yes, only then
    // assign the corresponding array index to the peptide and write it out.
    if (pbAuxLoc->location_size() > 0) {
      pbPeptide.set_aux_locations_index(auxLocIdx++);
      auxLocWriter.Write(pbAuxLoc);
      delete pbAuxLoc;
      pbAuxLoc = new pb::AuxLocation();
    }

    // Write the peptide AFTER the aux_locations check, in case we added an
    // aux_locations_index to the peptide.
    peptideWriter.Write(&pbPeptide);

    if (++count % 100000 == 0) {
      carp(CARP_INFO, "Wrote %d peptides", count);
    }
  }
  delete pbAuxLoc;
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
      aaMass = MassConstants::fixp_avg_table[sequence[i]];
      if (aaMass == 0) {
        return numeric_limits<double>::signaling_NaN();
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

void TideIndexApplication::getPbPeptide(
  int id,
  const TideIndexPeptide& peptide,
  pb::Peptide& outPbPeptide
) {
  outPbPeptide.set_id(id);
  outPbPeptide.set_mass(peptide.getMass());
  outPbPeptide.set_length(peptide.getLength());
  outPbPeptide.mutable_first_location()->set_protein_id(peptide.getProteinId());
  outPbPeptide.mutable_first_location()->set_pos(peptide.getProteinPos());;
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

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
