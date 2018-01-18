#include "GeneratePeptides.h"
#include "parameter.h"
#include "model/ProteinPeptideIterator.h"
#include "util/mass.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/StringUtils.h"

using namespace std;

MASS_TYPE_T GeneratePeptides::massType_ = AVERAGE;

GeneratePeptides::GeneratePeptides() {
}

GeneratePeptides::~GeneratePeptides() {
}

int GeneratePeptides::main(int argc, char** argv) {
  string fastaPath = Params::GetString("protein fasta file");
  if (!FileUtils::Exists(fastaPath)) {
    carp(CARP_FATAL, "File does not exist: '%s'", fastaPath.c_str());
  }

  DECOY_TYPE_T decoyType = get_tide_decoy_type_parameter("decoy-format");
  massType_ = get_mass_type_parameter("isotopic-mass");
  bool overwrite = Params::GetBool("overwrite");

  string targetsFile = make_file_path("generate-peptides.target.txt");
  string decoysFile = make_file_path("generate-peptides.decoy.txt");

  ofstream* targetPeptides = create_stream_in_path(targetsFile.c_str(), NULL, overwrite);
  ofstream* decoyPeptides = decoyType != NO_DECOYS
    ? create_stream_in_path(decoysFile.c_str(), NULL, overwrite)
    : NULL;

  string decoyFastaPath = canGenerateDecoyProteins()
    ? make_file_path("generate-peptides.proteins.decoy.txt")
    : "";

  processFasta(fastaPath, targetPeptides, decoyFastaPath, decoyPeptides, decoyType);
  return 0;
}

void GeneratePeptides::processFasta(
  const string& fastaPath,
  ofstream* targetList,
  const string& decoyFastaPath,
  ofstream* decoyList,
  DECOY_TYPE_T decoyType
) {
  double minMass = Params::GetDouble("min-mass");
  double maxMass = Params::GetDouble("max-mass");

  carp(CARP_INFO, "Reading %s", fastaPath.c_str());
  ifstream* fasta = new ifstream(fastaPath.c_str(), ifstream::in);

  ofstream* decoyFasta = decoyFastaPath.empty() ? NULL
    : create_stream_in_path(decoyFastaPath.c_str(), NULL, Params::GetBool("overwrite"));

  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  if (!Params::GetString("custom-enzyme").empty()) {
    enzyme = CUSTOM_ENZYME;
  }
  DIGEST_T digest = get_digest_type_parameter("digestion");
  int missed = Params::GetInt("missed-cleavages");
  int minLen = Params::GetInt("min-length");
  int maxLen = Params::GetInt("max-length");
  string decoyPrefix = Params::GetString("decoy-prefix");
  bool proteinReverse = decoyType == PROTEIN_REVERSE_DECOYS;
  bool peptideShuffle = decoyType == PEPTIDE_SHUFFLE_DECOYS;
  bool peptideReverse = decoyType == PEPTIDE_REVERSE_DECOYS;

  set<string> targets, decoys;
  map<const string*, const string*> targetToDecoy;
  map< OrderedPeptide, vector<string> > peptideToProtein;
  int proteinTotal = 0, peptideTotal = 0;

  // Iterate over all proteins from this FASTA
  while (true) {
    string id, proteinSequence;
    if (!getNextProtein(*fasta, &id, &proteinSequence)) {
      break;
    }
    ++proteinTotal;
    carp(CARP_DEBUG, "Read %s", id.c_str());

    // Get peptides
    vector<CleavedPeptide> peptides =
      cleaveProtein(proteinSequence, enzyme, digest, missed, minLen, maxLen);
    peptideTotal += peptides.size();

    // Write reversed protein to decoy FASTA, if protein reverse
    if (decoyFasta && proteinReverse) {
      reverse(proteinSequence.begin(), proteinSequence.end());
      *decoyFasta << '>' << decoyPrefix << id << endl
                  << proteinSequence << endl;
    }

    // Iterate over all peptides from this protein
    for (vector<CleavedPeptide>::const_iterator i = peptides.begin();
         i != peptides.end();
         i++) {
      FLOAT_T mass = i->Mass();
      if (mass < minMass || mass > maxMass) {
        carp(CARP_DETAILED_DEBUG, "Skipping peptide with mass %f", mass);
        continue;
      }
      const string& sequence = i->Sequence();
      pair<set<string>::iterator, bool> insert = targets.insert(sequence);
      if (insert.second) {
        peptideToProtein[*i] = vector<string>(1, id);
      } else {
        peptideToProtein[*i].push_back(id);
      }
    }
  }
  delete fasta;
  carp(CARP_DEBUG, "Read %d proteins and %d peptides", proteinTotal, peptideTotal);

  // Once we have all targets, try to generate a decoy for each one
  if (peptideShuffle || peptideReverse) {
    int decoyFailures = 0;
    for (set<string>::const_iterator i = targets.begin(); i != targets.end(); i++) {
      string decoy;
      if (makeDecoy(*i, targets, decoys, peptideShuffle, decoy)) {
        targetToDecoy[&*i] = &*(decoys.insert(decoy).first);
      } else {
        ++decoyFailures;
      }
    }
    if (decoyFailures > 0) {
      carp(CARP_WARNING, "Failed to generate decoys for %d targets", decoyFailures);
    }
  }

  // Re-read FASTA and generate decoy FASTA
  if (decoyFasta && (peptideShuffle || peptideReverse)) {
    ifstream* fasta = new ifstream(fastaPath.c_str(), ifstream::in);
    while (true) {
      string id, proteinSequence;
      if (!getNextProtein(*fasta, &id, &proteinSequence)) {
        break;
      }

      *decoyFasta << '>' << decoyPrefix << id << endl;

      // Get peptides
      vector<CleavedPeptide> peptides =
        cleaveProtein(proteinSequence, enzyme, digest, missed, 0, 1000000);

      // Iterate over all peptides from this protein
      for (vector<CleavedPeptide>::const_iterator i = peptides.begin();
           i != peptides.end();
           i++) {
        const string& sequence = i->Sequence();
        set<string>::const_iterator j = targets.find(sequence);
        if (j != targets.end()) {
          map<const string*, const string*>::const_iterator k = targetToDecoy.find(&*j);
          if (k != targetToDecoy.end()) {
            *decoyFasta << *k->second;
            continue;
          }
        }
        *decoyFasta << sequence;
      }
      *decoyFasta << endl;
    }
    delete fasta;
  }

  int precision = Params::GetInt("mass-precision");
  for (map< OrderedPeptide, vector<string> >::const_iterator i = peptideToProtein.begin();
       i != peptideToProtein.end();
       i++) {
    const string& sequence = i->first.Sequence();
    FLOAT_T mass = i->first.Mass();
    vector<string>::const_iterator j = i->second.begin();
    *targetList << sequence << '\t'
                << StringUtils::ToString(mass + MASS_PROTON, precision) << '\t'
                << *j;
    for (j = j + 1; j != i->second.end(); j++) {
      string proteinId = *j;
      *targetList << ',' << proteinId;
    }
    *targetList << endl;

    if ((peptideShuffle || peptideReverse) && decoyList) {
      set<string>::const_iterator target = targets.find(sequence);
      map<const string*, const string*>::const_iterator k = targetToDecoy.find(&*target);
      if (k != targetToDecoy.end()) {
        j = i->second.begin();
        *decoyList << *k->second << '\t'
                   << StringUtils::ToString(mass + MASS_PROTON, precision) << '\t'
                   << decoyPrefix << *j;
        for (j = j + 1; j != i->second.end(); j++) {
          *decoyList << ',' << decoyPrefix << *j;
        }
        *decoyList << endl;
      }
    }
  }

  if (decoyFasta) {
    delete decoyFasta;
  }

  delete targetList;

  if (proteinReverse && !decoyFastaPath.empty()) {
    processFasta(decoyFastaPath, decoyList, "", NULL, NO_DECOYS);
  } else if (decoyList) {
    delete decoyList;
  }
}

bool GeneratePeptides::canGenerateDecoyProteins() {
  const string decoyFormat = Params::GetString("decoy-format");

  // Can never write decoy proteins if not making decoys
  if (decoyFormat == "none") {
    return false;
  }

  // Can always write decoy proteins if making protein-level decoys
  if (decoyFormat == "protein-reverse") {
    return true;
  }

  // If making peptide-level decoys, we can only write decoy proteins if:
  // Using an enzyme, full digestion, no missed cleavages
  bool customEnzyme = !Params::GetString("custom-enzyme").empty();
  bool useEnzyme = get_enzyme_type_parameter("enzyme") != NO_ENZYME;
  bool fullDigest = get_digest_type_parameter("digestion") == FULL_DIGEST;
  bool noMissedCleavages = Params::GetInt("missed-cleavages") == 0;

  return (customEnzyme || useEnzyme) && fullDigest && noMissedCleavages;
}

/**
 * Reads the next protein ID and corresponding sequence from the FASTA stream
 * Returns false if no more proteins in stream
 */
bool GeneratePeptides::getNextProtein(
  ifstream& fasta,  ///< FASTA stream
  string* outId,  ///< string to store protein ID
  string* outSequence ///< string to store sequence
) {
  outId->clear();
  outSequence->clear();
  if (!fasta.good()) {
    return false;
  }

  while (fasta.good()) {
    string line;
    getline(fasta, line);
    line = StringUtils::Trim(line);
    if (outId->empty()) {
      // Reading id
      if (StringUtils::StartsWith(line, ">")) {
        bool idStart = false;
        string::const_iterator begin = line.begin() + 1;
        string::const_iterator end = line.end();
        for (string::const_iterator i = begin; i != line.end(); i++) {
          bool space = isspace(*i);
          if (!idStart && !space) {
            begin = i;
            idStart = true;
          } else if (idStart && space) {
            end = i;
            break;
          }
        }
        *outId = string(begin, end);
      }
    } else {
      // Reading sequence
      *outSequence += line;
      if (fasta.eof() || fasta.peek() == '>') {
        break;
      }
    }
  }

  if (StringUtils::EndsWith(*outSequence, "*")) {
    // Remove the last character of the sequence if it is an asterisk
    outSequence->erase(outSequence->length() - 1);
  }

  if (outSequence->empty()) {
    carp(CARP_WARNING, "Found protein ID without sequence: %s", outId->c_str());
    outId->clear();
    return false;
  }

  if (fasta.fail()) {
    carp(CARP_FATAL, "Error reading FASTA file");
  }

  return true;
}

/**
 * Cleave protein sequence using specified enzyme and store results in vector
 * Vector also contains start location of each peptide within the protein
 */
vector<GeneratePeptides::CleavedPeptide> GeneratePeptides::cleaveProtein(
  const string& sequence, ///< Protein sequence to cleave
  ENZYME_T enzyme,  ///< Enzyme to use for cleavage
  DIGEST_T digest,  ///< Digestion to use for cleavage
  int missedCleavages,  ///< Maximum allowed missed cleavages
  int minLength,  //< Min length of peptides to return
  int maxLength  //< Max length of peptides to return
) {
  vector<CleavedPeptide> peptides;
  // No enzyme
  // Get all substrings min <= length <= max
  if (enzyme == NO_ENZYME) {
    for (int i = 0; i < sequence.length(); ++i) {
      for (int j = minLength; i + j <= sequence.length() && j <= maxLength; ++j) {
        peptides.push_back(CleavedPeptide(sequence.substr(i, j), i));
      }
    }
    return peptides;
  }

  size_t pepStart = 0, nextPepStart = 0;
  int cleaveSites = 0;
  for (int i = 0; i < sequence.length(); ++i) {
    // Determine if this is a valid cleavage position
    bool cleavePos = i != sequence.length() - 1 &&
      ProteinPeptideIterator::validCleavagePosition(sequence.c_str() + i, enzyme);
    if (digest == PARTIAL_DIGEST && i != sequence.length() - 1 && !cleavePos) {
      // Partial digestion (not last AA or cleavage position), add this peptide
      peptides.push_back(CleavedPeptide(sequence.substr(pepStart, i + 1 - pepStart), pepStart));
    } else if (cleavePos) {
      // Cleavage position, add this peptide
      peptides.push_back(CleavedPeptide(sequence.substr(pepStart, i + 1 - pepStart), pepStart));
      if (Params::GetBool("clip-nterm-methionine") && sequence[0] == 'M' &&
          pepStart == 0 && digest != PARTIAL_DIGEST) {
        peptides.push_back(CleavedPeptide(sequence.substr(1, i), 1));
      }
      if (++cleaveSites == 1) {
        // This is the first cleavage position, remember it
        nextPepStart = i + 1;
      }
      if (digest == PARTIAL_DIGEST) {
        // For partial digest, add peptides ending at this cleavage position
        for (int j = pepStart + 1; j < nextPepStart; ++j) {
          peptides.push_back(CleavedPeptide(sequence.substr(j, i - j + 1), j));
        }
      }
      if (cleaveSites > missedCleavages) {
        // We have missed the allowed amount of cleavages
        // Move iterator+pepStart to the first cleavage position
        pepStart = nextPepStart;
        i = pepStart - 1;
        cleaveSites = 0;
      }
    } else if (i == sequence.length() - 1 &&
               cleaveSites > 0 && cleaveSites <= missedCleavages) {
      // Last AA in sequence and we haven't missed the allowed amount yet
      // Add this peptide and move iterator+pepStart to first cleavage position
      peptides.push_back(CleavedPeptide(sequence.substr(pepStart), pepStart));
      if (digest == PARTIAL_DIGEST) {
        // For partial digest, add peptides ending at last AA
        for (int j = pepStart + 1; j < nextPepStart; ++j) {
          peptides.push_back(CleavedPeptide(sequence.substr(j, i - j + 1), j));
        }
      }
      pepStart = nextPepStart;
      i = pepStart - 1;
      cleaveSites = 0;
    }
  }
  //For peptides that do not have an internal cleavage 
  if (Params::GetBool("clip-nterm-methionine") && sequence[0] == 'M' &&
      pepStart == 0 && digest != PARTIAL_DIGEST) {
    peptides.push_back(CleavedPeptide(sequence.substr(1), 1));
  }
  // Add the last peptide
  peptides.push_back(CleavedPeptide(sequence.substr(nextPepStart), nextPepStart));
  if (digest == PARTIAL_DIGEST) {
    // For partial digest, add peptides ending at last AA
    for (int j = pepStart + 1; j < sequence.length(); ++j) {
      peptides.push_back(CleavedPeptide(sequence.substr(j), j));
    }
  }
  // Erase peptides that don't meet length requirement
  for (vector<CleavedPeptide>::iterator i = peptides.begin(); i != peptides.end(); ) {
    if (i->Length() < minLength || i->Length() > maxLength) {
      i = peptides.erase(i);
    } else {
      ++i;
    }
  }
  return peptides;
}

/**
 * Makes a decoy from the sequence.
 * Returns false on failure, and decoyOut will be the same as seq.
 */
bool GeneratePeptides::makeDecoy(
  const string& seq,  ///< sequence to make decoy from
  const set<string>& targetSeqs,  ///< targets to check against
  const set<string>& decoySeqs,  ///< decoys to check against
  bool shuffle, ///< shuffle (if false, reverse)
  string& decoyOut  ///< string to store decoy
) {
  const string keepTerminal = Params::GetString("keep-terminal-aminos");
  string decoyPre, decoyPost;
  if (keepTerminal == "N") {
    if (seq.length() <= 2) {
      decoyOut = seq;
      return false;
    }
    decoyPre = seq[0];
    decoyOut = seq.substr(1);
  } else if (keepTerminal == "C") {
    if (seq.length() <= 2) {
      decoyOut = seq;
      return false;
    }
    decoyPost = seq[seq.length() - 1];
    decoyOut = seq.substr(0, seq.length() - 1);
  } else if (keepTerminal == "NC") {
    if (seq.length() <= 3) {
      decoyOut = seq;
      return false;
    }
    decoyPre = seq[0];
    decoyPost = seq[seq.length() - 1];
    decoyOut = seq.substr(1, seq.length() - 2);
  } else {
    decoyOut = seq;
    if (seq.length() <= 1) {
      return false;
    }
  }

  if (!shuffle) {
    // Reverse
    if (reversePeptide(decoyOut)) {
      // Re-add n/c
      string decoyCheck = decoyPre + decoyOut + decoyPost;
      // Check in sets
      if (targetSeqs.find(decoyCheck) == targetSeqs.end() &&
          decoySeqs.find(decoyCheck) == decoySeqs.end()) {
        decoyOut = decoyCheck;
        return true;
      }
    }
    carp(CARP_DEBUG, "Failed reversing %s, shuffling", seq.c_str());
  }

  // Shuffle
  if (shufflePeptide(decoyOut)) {
    // Re-add n/c
    string decoyCheck = decoyPre + decoyOut + decoyPost;
    // Check in sets
    if (targetSeqs.find(decoyCheck) == targetSeqs.end() &&
        decoySeqs.find(decoyCheck) == decoySeqs.end()) {
      decoyOut = decoyCheck;
      return true;
    }
  }

  decoyOut = seq;
  return false;
}

/**
 * Shuffles the peptide randomly.
 * Returns false if no different sequence was generated
 */
bool GeneratePeptides::shufflePeptide(
  string& seq, ///< Peptide sequence to shuffle
  unsigned int maxShuffleAttempts ///< Maximum number of shuffle attempts
) {
  switch (seq.length()) {
  case 0:
  case 1:
    return false;
  case 2: {
    char second = seq[1];
    seq[1] = seq[0];
    seq[0] = second;
    return true;
  }
  default:
    string originalSeq(seq);
    for (int i = 0; i < maxShuffleAttempts; i++) {
      random_shuffle(seq.begin(), seq.end(), myrandom_limit);
      if (seq != originalSeq) {
        return true;
      }
    }
    return false;
  }
}

/**
 * Reverses the peptide sequence.
 * Returns false if no different sequence was generated
 */
bool GeneratePeptides::reversePeptide(
  string& seq ///< Peptide sequence to reverse
) {
  string originalSeq(seq);
  reverse(seq.begin(), seq.end());
  return seq != originalSeq;
}

string GeneratePeptides::getName() const {
  return "generate-peptides";
}

string GeneratePeptides::getDescription() const {
  return
    "[[nohtml:Extract from a given set of protein sequences a list of target "
    "and decoy peptides fitting the specified criteria.]]"
    "[[html:<p>This command takes as input a protein FASTA file and outputs "
    "the corresponding list of peptides, as well as a matched list of decoy "
    "peptides and decoy proteins. Decoys are generated either by reversing or "
    "shuffling the non-terminal amino acids of each peptide. The program will "
    "shuffle each peptide multiple times to attempt to ensure that there is no "
    "overlap between the target and decoy peptides. For homopolymers, this is "
    "not possible. In this case, the occurrence of these target/decoy overlaps "
    "is recorded in the log file.</p><p>The program considers only the "
    "standard set of 20 amino acids. Peptides containing non-amino acid "
    "alphanumeric characters (BXZ) are skipped. Non-alphanumeric characters "
    "are ignored completely.</p>]]";
}

vector<string> GeneratePeptides::getArgs() const {
  string arr[] = {
    "protein fasta file"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> GeneratePeptides::getOptions() const {
  string arr[] = {
    "min-mass",
    "max-mass",
    "min-length",
    "max-length",
    "enzyme",
    "custom-enzyme",
    "digestion",
    "missed-cleavages",
    "isotopic-mass",
    "seed",
    "clip-nterm-methionine",
    "decoy-format",
    "decoy-prefix",
    "keep-terminal-aminos",
    "mod-precision",
    "overwrite",
    "fileroot",
    "output-dir",
    "parameter-file",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > GeneratePeptides::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("generate-peptides.target.txt",
    "A text file containing the target peptides, one per line. Each line has "
    "three tab-delimited columns, containing the peptide sequence, the m+h "
    "mass of the unmodified peptide, and a comma-delimited list of protein IDs "
    "in which the peptide occurs."));
  outputs.push_back(make_pair("generate-peptides.decoy.txt",
    "A text file containing the decoy peptides, one per line. Each line has "
    "three tab-delimited columns, containing the peptide sequence, the m+h "
    "mass of the unmodified peptide, and a comma-delimited list of protein IDs "
    "in which the peptide occurs. "
    "There is a one-to-one correspondence between targets and decoys."));
  outputs.push_back(make_pair("generate-peptides.proteins.decoy.txt",
    "a FASTA format file containing decoy proteins, in which all of the "
    "peptides have been replaced with their shuffled or reversed counterparts. "
    "Note that this file will only be created if the enzyme specificity is "
    "\"full-digest\" and no missed cleavages are allowed."));
  outputs.push_back(make_pair("generate-peptides.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("generate-peptides.log.txt",
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution."));
  return outputs;
}

bool GeneratePeptides::needsOutputDirectory() const {
  return true;
}

COMMAND_T GeneratePeptides::getCommand() const {
  return GENERATE_PEPTIDES_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
