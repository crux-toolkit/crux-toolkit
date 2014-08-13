#include "GenerateDecoys.h"
#include "parameter.h"
#include "ProteinPeptideIterator.h"

using namespace std;

MASS_TYPE_T GenerateDecoys::massType_ = AVERAGE;

GenerateDecoys::GenerateDecoys() {
}

GenerateDecoys::~GenerateDecoys() {
}

int GenerateDecoys::main(int argc, char** argv) {

  const char* option_list[] = {
    "min-mass",
    "max-mass",
    "min-length",
    "max-length",
    "enzyme",
    "custom-enzyme",
    "digestion",
    "missed-cleavages",
    "isotopic-mass",
    "clip-nterm-methionine",
    "decoy-format",
    "keep-terminal-aminos",
    "overwrite",
    "fileroot",
    "output-dir",
    "parameter-file",
    "verbosity"
  };

  int num_options = sizeof(option_list) / sizeof(char*);

  const char* argument_list[] = { "protein fasta file" };
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize(argument_list, num_arguments, option_list, num_options, argc, argv);

  const int MAX_SHUFFLE_ATTEMPTS = 10;

  // Get decoy type
  const DECOY_TYPE_T decoyType = get_tide_decoy_type_parameter("decoy-format");
  bool proteinReverse = decoyType == PROTEIN_REVERSE_DECOYS;

  // Get options
  double minMass = get_double_parameter("min-mass");
  double maxMass = get_double_parameter("max-mass");
  massType_ = get_mass_type_parameter("isotopic-mass");

  bool overwrite = get_boolean_parameter("overwrite");

  string targetsFile = make_file_path("peptides.target.txt");
  string decoysFile = make_file_path("peptides.decoy.txt");
  string proteinDecoysFile = make_file_path("proteins.decoy.txt");
  ofstream* targetsStream = create_stream_in_path(targetsFile.c_str(), NULL, overwrite);
  ofstream* decoysStream = decoyType != NO_DECOYS ?
    create_stream_in_path(decoysFile.c_str(), NULL, overwrite) : NULL;
  ofstream* proteinDecoysStream = canGenerateDecoyProteins() ?
    create_stream_in_path(proteinDecoysFile.c_str(), NULL, overwrite) : NULL;

  // Read fasta
  string fastaFile = get_string_parameter_pointer("protein fasta file");
  carp(CARP_INFO, "Reading %s", fastaFile.c_str());
  map< string, vector<string> > proteins;
  set<string> targetSeqs;
  set<string> decoySeqs;
  readFasta(fastaFile, proteins, targetSeqs,
            proteinReverse ? proteinDecoysStream : NULL,
            proteinReverse ? &decoySeqs : NULL);

  // Make decoys from targets and write to peptides files
  if (decoysStream) {
    carp(CARP_INFO, "Making decoys and writing peptides files");
  } else {
    carp(CARP_INFO, "Writing peptides file");
  }
  map<string, const string*> targetToDecoy;
  for (set<string>::const_iterator i = targetSeqs.begin();
       i != targetSeqs.end();
       ++i) {
    const string& targetSeq = *i;
    // Don't need to check peptide length, it is done in cleaveProtein
    // Check peptide mass
    FLOAT_T pepMass = Crux::Peptide::calcSequenceMass(i->c_str(), massType_);
    if (pepMass < minMass || pepMass > maxMass) {
      carp(CARP_DETAILED_DEBUG, "Skipping peptide with mass %f", pepMass);
      continue;
    }
    (*targetsStream) << *i << endl;
    if (decoysStream && !proteinReverse) {
      // Try to make decoy
      string decoySeq;
      if (makeDecoy(targetSeq, targetSeqs, decoySeqs,
                    decoyType == PEPTIDE_SHUFFLE_DECOYS, decoySeq)) {
        // Success
        pair<set<string>::iterator, bool> decoyInsert = decoySeqs.insert(decoySeq);
        targetToDecoy[targetSeq] = &(*(decoyInsert.first));
      } else {
        carp(CARP_WARNING, "Could not make decoy from %s", targetSeq.c_str());
      }
      (*decoysStream) << decoySeq << endl;
    }
  }

  // Write decoy peptides for protein shuffle
  if (proteinReverse) {
    for (set<string>::const_iterator i = decoySeqs.begin();
         i != decoySeqs.end();
         ++i) {
      FLOAT_T pepMass = Crux::Peptide::calcSequenceMass(i->c_str(), massType_);
      if (pepMass < minMass || pepMass > maxMass) {
        carp(CARP_DETAILED_DEBUG, "Skipping peptide with mass %f", pepMass);
        continue;
      }
      (*decoysStream) << *i << endl;
    }
  }

  targetsStream->close();
  delete targetsStream;
  if (decoysStream) {
    decoysStream->close();
    delete decoysStream;
    // Write decoy proteins (unless protein-shuffle, because we already did)
    if (proteinDecoysStream && !proteinReverse) {
      carp(CARP_INFO, "Writing decoy proteins");
      for (map< string, vector<string> >::iterator proteinIter = proteins.begin();
           proteinIter != proteins.end();
           ++proteinIter) {
        (*proteinDecoysStream) << ">decoy_" << proteinIter->first << endl;
        for (vector<string>::iterator pepIter = proteinIter->second.begin();
             pepIter != proteinIter->second.end();
             ++pepIter) {
          map<string, const string*>::const_iterator lookup = targetToDecoy.find(*pepIter);
          const string* toOutput = (lookup != targetToDecoy.end()) ?
            lookup->second : &(*pepIter);
          (*proteinDecoysStream) << *toOutput;
        }
        (*proteinDecoysStream) << endl;
      }
      carp(CARP_DEBUG, "Printed %d decoy proteins", proteins.size());
      proteinDecoysStream->close();
      delete proteinDecoysStream;
    }
  }

  return 0;
}

bool GenerateDecoys::canGenerateDecoyProteins() {
  const string decoyFormat = get_string_parameter_pointer("decoy-format");

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
  bool customEnzyme =
    string(get_string_parameter_pointer("custom-enzyme")) == "__NULL_STR";
  bool useEnzyme = get_enzyme_type_parameter("enzyme") != NO_ENZYME;
  bool fullDigest = get_digest_type_parameter("digestion") == FULL_DIGEST;
  bool noMissedCleavages = get_int_parameter("missed-cleavages") == 0;

  return (customEnzyme || useEnzyme) && fullDigest && noMissedCleavages;
}

/**
 * Given a FASTA file, read in all protein IDs/sequences and cleave them.
 * Return a map of protein IDs to digested peptides from that protein
 */
void GenerateDecoys::readFasta(
  const string& fastaName,  ///< FASTA file name
  map< string, vector<string> >& outProteins, ///< map to store proteins
  set<string>& outPeptides,  ///< set of unique peptides
  ofstream* reversedFasta, ///< optional stream to write reversed proteins
  set<string>* outReversedPeptides  ///< optional set of peptides from rev fasta
) {
  // Open FASTA
  ifstream fasta(fastaName.c_str(), ifstream::in);

  // Read all proteins
  outProteins.clear();
  outPeptides.clear();
  string id, sequence;
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digest = get_digest_type_parameter("digestion");
  int missedCleavages = get_int_parameter("missed-cleavages");
  int minLength = get_int_parameter("min-length");
  int maxLength = get_int_parameter("max-length");
  if (string(get_string_parameter_pointer("custom-enzyme")) != "__NULL_STR") {
    enzyme = CUSTOM_ENZYME;
  }
  vector<string> trypticPeptides;
  vector<string> reversedFastaPeptides;
  int proteinTotal = 0;
  int peptideTotal = 0;
  while (getNextProtein(fasta, id, sequence)) {
    ++proteinTotal;
    carp(CARP_DEBUG, "Read %s", id.c_str());
    cleaveProtein(sequence, enzyme, digest, missedCleavages,
                  minLength, maxLength, trypticPeptides);
    peptideTotal += trypticPeptides.size();
    copy(trypticPeptides.begin(), trypticPeptides.end(),
         inserter(outPeptides, outPeptides.end()));
    outProteins[id] = trypticPeptides;
    if (reversedFasta) {
      reverse(sequence.begin(), sequence.end());
      (*reversedFasta) << ">decoy_" << id << '\n' << sequence << endl;
      if (outReversedPeptides) {
        cleaveProtein(sequence, enzyme, digest, missedCleavages,
                      minLength, maxLength, reversedFastaPeptides);
        copy(reversedFastaPeptides.begin(), reversedFastaPeptides.end(),
             inserter(*outReversedPeptides, outReversedPeptides->end()));
      }
    }
  }
  fasta.close();

  carp(CARP_DEBUG, "Read %d proteins and %d peptides", proteinTotal, peptideTotal);
}

/**
 * Reads the next protein ID and corresponding sequence from the FASTA stream
 * Returns false if no more proteins in stream
 */
bool GenerateDecoys::getNextProtein(
  ifstream& fasta,  ///< FASTA stream
  string& outId,  ///< string to store protein ID
  string& outSequence ///< string to store sequence
) {
  const string whitespace = " \t\n\v\f\r";
  outId.clear();
  outSequence.clear();
  if (!fasta.good()) {
    return false;
  }

  while (fasta.good()) {
    string line;
    getline(fasta, line);
    // Trim whitespace
    size_t first = line.find_first_not_of(whitespace);
    if (first != string::npos) {
      line.erase(0, first);
      size_t last = line.find_last_not_of(whitespace);
      line.erase(last + 1);
    } else {
      line.clear();
    }
    if (outId.empty()) {
      // Reading id
      if (line.length() > 0 && line[0] == '>') {
        outId = line.substr(1);
      }
    } else {
      // Reading sequence
      outSequence += line;
      if (fasta.eof() || fasta.peek() == '>') {
        break;
      }
    }
  }

  if (!outSequence.empty() && outSequence[outSequence.length() - 1] == '*') {
    // Remove the last character of the sequence if it is an asterisk
    outSequence.erase(outSequence.length() - 1);
  }

  if (outSequence.empty()) {
    carp(CARP_WARNING, "Found protein ID without sequence: %s", outId.c_str());
    outId.clear();
    return false;
  }

  if (fasta.fail()) {
    carp(CARP_FATAL, "Error reading FASTA file");
  }

  return true;
}

/**
 * Cleave protein sequence using specified enzyme and store results in vector
 */
  void GenerateDecoys::cleaveProtein(
  const string& sequence, ///< Protein sequence to cleave
  ENZYME_T enzyme,  ///< Enzyme to use for cleavage
  DIGEST_T digest,  ///< Digestion to use for cleavage
  int missedCleavages,  ///< Maximum allowed missed cleavages
  int minLength,  //< Min length of peptides to return
  int maxLength,  //< Max length of peptides to return
  vector<string>& outPeptides ///< vector to store peptides
) {
  outPeptides.clear();
  bool clip_nterm_methionine = get_boolean_parameter("clip-nterm-methionine");
  if (enzyme != NO_ENZYME) {
    // Enzyme
    size_t pepStart = 0, nextPepStart = 0;
    int cleaveSites = 0;
    for (int i = 0; i < sequence.length(); ++i) {
      // Determine if this is a valid cleavage position
      bool cleavePos =
        ProteinPeptideIterator::validCleavagePosition(sequence.c_str() + i, enzyme);
      if (i != sequence.length() - 1 && !cleavePos && digest == PARTIAL_DIGEST) {
        // Partial digestion (not last AA or cleavage position), add this peptide
        outPeptides.push_back(sequence.substr(pepStart, i + 1 - pepStart));
      } else if (cleavePos) {
        // Cleavage position, add this peptide
        outPeptides.push_back(sequence.substr(pepStart, i + 1 - pepStart));
        if (clip_nterm_methionine && sequence[0] == 'M' && pepStart == 0 && digest != PARTIAL_DIGEST) {
          outPeptides.push_back(sequence.substr(pepStart + 1, i + 1 - pepStart - 1));
        }
        if (++cleaveSites == 1) {
          // This is the first cleavage position, remember it
          nextPepStart = i + 1;
        }
        if (digest == PARTIAL_DIGEST) {
          // For partial digest, add peptides ending at this cleavage position
          for (int j = pepStart + 1; j < nextPepStart; ++j) {
            outPeptides.push_back(sequence.substr(j, i - j + 1));
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
        outPeptides.push_back(sequence.substr(pepStart));
        if (digest == PARTIAL_DIGEST) {
          // For partial digest, add peptides ending at last AA
          for (int j = pepStart + 1; j < nextPepStart; ++j) {
            outPeptides.push_back(sequence.substr(j, i - j + 1));
          }
        }
        pepStart = nextPepStart;
        i = pepStart - 1;
        cleaveSites = 0;
      }
    }
    // Add the last peptide
    outPeptides.push_back(sequence.substr(nextPepStart));
    if (digest == PARTIAL_DIGEST) {
      // For partial digest, add peptides ending at last AA
      for (int j = pepStart + 1; j < sequence.length(); ++j) {
        outPeptides.push_back(sequence.substr(j));
      }
    }
    // Erase peptides that don't meet length requirement
    for (vector<string>::reverse_iterator i = outPeptides.rbegin();
         i != outPeptides.rend();
         ++i) {
      if (i->length() < minLength || i->length() > maxLength) {
        outPeptides.erase((i + 1).base());
      }
    }
  } else {
    // No enzyme
    // Get all substrings min <= length <= max
    for (int i = 0; i < sequence.length(); ++i) {
      for (int j = minLength; i + j <= sequence.length() && j <= maxLength; ++j) {
        outPeptides.push_back(sequence.substr(i, j));
      }
    }
  }
}

/**
 * Cleave protein sequence using specified enzyme and store results in vector
 * Vector also contains start location of each peptide within the protein
 */
void GenerateDecoys::cleaveProtein(
  const string& sequence, ///< Protein sequence to cleave
  ENZYME_T enzyme,  ///< Enzyme to use for cleavage
  DIGEST_T digest,  ///< Digestion to use for cleavage
  int missedCleavages,  ///< Maximum allowed missed cleavages
  int minLength,  //< Min length of peptides to return
  int maxLength,  //< Max length of peptides to return
  vector< pair<string, int> >& outPeptides ///< vector to store peptides
) {
  outPeptides.clear();
  bool clip_nterm_methionine = get_boolean_parameter("clip-nterm-methionine");
  if (enzyme != NO_ENZYME) {
    // Enzyme
    size_t pepStart = 0, nextPepStart = 0;
    int cleaveSites = 0;
    for (int i = 0; i < sequence.length(); ++i) {
      // Determine if this is a valid cleavage position
      bool cleavePos =
        ProteinPeptideIterator::validCleavagePosition(sequence.c_str() + i, enzyme);
      if (i != sequence.length() - 1 && !cleavePos && digest == PARTIAL_DIGEST) {
        // Partial digestion (not last AA or cleavage position), add this peptide
        outPeptides.push_back(
          make_pair(sequence.substr(pepStart, i + 1 - pepStart), pepStart));
      } else if (cleavePos) {
        // Cleavage position, add this peptide
        outPeptides.push_back(
          make_pair(sequence.substr(pepStart, i + 1 - pepStart), pepStart));
        if (clip_nterm_methionine && sequence[0] == 'M' && pepStart == 0 && digest != PARTIAL_DIGEST) {
          outPeptides.push_back(
            make_pair(sequence.substr(pepStart + 1, i + 1 - pepStart - 1), pepStart + 1));
        }
        if (++cleaveSites == 1) {
          // This is the first cleavage position, remember it
          nextPepStart = i + 1;
        }
        if (digest == PARTIAL_DIGEST) {
          // For partial digest, add peptides ending at this cleavage position
          for (int j = pepStart + 1; j < nextPepStart; ++j) {
            outPeptides.push_back(make_pair(sequence.substr(j, i - j + 1), j));
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
        outPeptides.push_back(make_pair(sequence.substr(pepStart), pepStart));
        if (digest == PARTIAL_DIGEST) {
          // For partial digest, add peptides ending at last AA
          for (int j = pepStart + 1; j < nextPepStart; ++j) {
            outPeptides.push_back(make_pair(sequence.substr(j, i - j + 1), j));
          }
        }
        pepStart = nextPepStart;
        i = pepStart - 1;
        cleaveSites = 0;
      }
    }
    // Add the last peptide
    outPeptides.push_back(make_pair(sequence.substr(nextPepStart), nextPepStart));
    if (digest == PARTIAL_DIGEST) {
      // For partial digest, add peptides ending at last AA
      for (int j = pepStart + 1; j < sequence.length(); ++j) {
        outPeptides.push_back(make_pair(sequence.substr(j), j));
      }
    }
    // Erase peptides that don't meet length requirement
    for (vector< pair<string, int> >::iterator i = outPeptides.begin();
         i != outPeptides.end();
    ) {
      if (i->first.length() < minLength || i->first.length() > maxLength) {
        i = outPeptides.erase(i);
      } else {
        ++i;
      }
    }
  } else {
    // No enzyme
    // Get all substrings min <= length <= max
    for (int i = 0; i < sequence.length(); ++i) {
      for (int j = minLength; i + j <= sequence.length() && j <= maxLength; ++j) {
        outPeptides.push_back(make_pair(sequence.substr(i, j), i));
      }
    }
  }
}

/**
 * Makes a decoy from the sequence.
 * Returns false on failure, and decoyOut will be the same as seq.
 */
bool GenerateDecoys::makeDecoy(
  const string& seq,  ///< sequence to make decoy from
  const set<string>& targetSeqs,  ///< targets to check against
  const set<string>& decoySeqs,  ///< decoys to check against
  bool shuffle, ///< shuffle (if false, reverse)
  string& decoyOut  ///< string to store decoy
) {
  const int MAX_SHUFFLE_ATTEMPTS = 6;

  const string keepTerminal = get_string_parameter_pointer("keep-terminal-aminos");
  string decoyPre;
  string decoyPost;
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
  for (int i = 0; i < MAX_SHUFFLE_ATTEMPTS; ++i) {
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
  }

  decoyOut = seq;
  return false;
}

/**
 * Shuffles the peptide randomly.
 * Returns false if no different sequence was generated
 */
bool GenerateDecoys::shufflePeptide(
  string& seq ///< Peptide sequence to shuffle
) {
  // Special case, length 2
  if (seq.length() == 2) {
    char second = seq[1];
    seq[1] = seq[0];
    seq[0] = second;
    return true;
  }

  string originalSeq(seq);
  random_shuffle(seq.begin(), seq.end(), myrandom_limit);
  return seq != originalSeq;
}

/**
 * Reverses the peptide sequence.
 * Returns false if no different sequence was generated
 */
bool GenerateDecoys::reversePeptide(
  string& seq ///< Peptide sequence to reverse
) {
  string originalSeq(seq);
  reverse(seq.begin(), seq.end());
  return seq != originalSeq;
}

string GenerateDecoys::getName() {
  return "generate-decoys";
}

string GenerateDecoys::getDescription() {
  return "Generates a corresponding list of peptides, as well as a matched "
         "list of decoy peptides and decoy proteins from a FASTA file";
}

bool GenerateDecoys::needsOutputDirectory() {
  return true;
}

COMMAND_T GenerateDecoys::getCommand() {
  return GENERATE_DECOYS_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
