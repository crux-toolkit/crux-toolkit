/*************************************************************************//**
 * \file Peptide.cpp
 * \brief Object for representing a single peptide.
 ****************************************************************************/
#include "Peptide.h"
#include "PeptideSrc.h"
#include "PostProcessProtein.h"
#include <string.h>

#include <numeric>
#include <set>
#include <vector>
#include "io/MatchFileReader.h"
#include "util/AminoAcidUtil.h"
#include "util/GlobalParams.h"
#include "util/StringUtils.h"

using namespace std;
using namespace Crux;

/**
 * \struct residue_iterator
 * \brief Object to iterate over the residues in a peptide, starting at the
 * first residue of the peptide, and proceeding in order.
 */
struct residue_iterator {
  Peptide*  peptide; ///< The peptide whose residues to iterate over.
  char*   sequence;    ///< The peptide sequence
  int     residue_idx; ///< The index of the current peak
};

/* Private functions */

/* Public functions--Allocators/Deallocators */

/**
 * \returns An (empty) peptide object.
 */
Peptide::Peptide() : length_(0), decoy_modified_seq_(NULL) {
}

Peptide::Peptide(string sequence)
  : sequence_(sequence), length_(sequence.length()), decoy_modified_seq_(NULL) {
}

Peptide::Peptide(string sequence, vector<Modification> mods)
  : sequence_(sequence), varMods_(mods), length_(sequence.length()),
    decoy_modified_seq_(NULL) {
}

// FIXME association part might be need to change
/**
 * \returns A new peptide object, populated with the user specified parameters.
 */
Peptide::Peptide(
  unsigned char length,     ///< The length of the peptide -in
  Protein* parent_protein, ///< the parent_protein of this peptide -in
  int start_idx ///< the start index of this peptide in the protein sequence -in
  )
  : length_(length), decoy_modified_seq_(NULL) {
  // FIXME: find the level of digest for this specific protein
  peptide_srcs_.push_back(new PeptideSrc(NON_SPECIFIC_DIGEST, parent_protein, start_idx));
}

/**
 * \brief Allocates a new peptide giving it the values of the source
 * peptide.
 * \returns A newly allocated peptide identical to the source.
 */
Peptide::Peptide(
  Peptide* src ///< source peptide -in
) {
  if (!src) {
    carp(CARP_ERROR, "Cannot copy null peptide!");
  } else {
    length_ = src->length_;
    PeptideSrc::copy(src->peptide_srcs_,peptide_srcs_);
    varMods_ = src->varMods_;
    if( src->decoy_modified_seq_ == NULL ){
      decoy_modified_seq_ = NULL;
    } else {
      decoy_modified_seq_ = copy_mod_aa_seq(src->decoy_modified_seq_, src->length_);
    }
  }
}

/**
 * Merges two identical peptides by adding the peptide_src of the
 * second to the first.  The second peptide remains unchanged.
 * Does not comfirm identity of peptides.
 * \returns true if merge is successfull.
 */
bool Peptide::mergePeptidesCopySrc(
  Peptide* peptide_dest,
  Peptide* peptide_giver
  ){

  vector<PeptideSrc*>& dest_srcs = peptide_dest->peptide_srcs_;
  vector<PeptideSrc*>& giver_srcs = peptide_giver->peptide_srcs_;

  // do both peptides have at least one peptide_src?
  if(dest_srcs.empty() || giver_srcs.empty()){
    carp(CARP_ERROR, "failed to merge two peptides");
    return false;
  }
 
  for(vector<PeptideSrc*>::iterator iter=giver_srcs.begin();
    iter != giver_srcs.end();
    iter++
  ){
    PeptideSrc* new_src = new PeptideSrc(*(*iter));
    dest_srcs.push_back(new_src);
  }

  return true;
}

/**
 * Frees an allocated peptide object.
 * Depending on peptide_src implementation determines how to free srcs
 */
Peptide::~Peptide() {
  for (vector<PeptideSrc*>::iterator i = peptide_srcs_.begin();
       i != peptide_srcs_.end();
       i++) {
    delete *i;
  }
  if(decoy_modified_seq_){
    freeModSeq(decoy_modified_seq_);
  }
}

// Public functions--Getters and Setters 

/* source-related getters and setters */

/**
 * sets the peptide_src field in the peptide
 * this method should be ONLY used when the peptide has no existing list of peptide_src
 * must pass on a heap allocated peptide_srcs_ object
 * does not copy in the object, just the pointer to the object.
 */
void Peptide::setPeptideSrc(
  PeptideSrc*  new_association ///< new peptide_src -in
  ) {
  assert(peptide_srcs_.empty());
  peptide_srcs_.push_back(new_association);
}

/**
 * this method adds the new_association to the end of the existing peptide's 
 * if no prior existing list, adds it at the front
 * must pass on a heap allocated peptide_srcs_ object
 * does not copy in the object, just the pointer to the object.
 */
void Peptide::addPeptideSrc(
  PeptideSrc* new_association ///< new peptide_src -in
  ) {
  peptide_srcs_.push_back(new_association);
}

// TODO: why do we need both of these?
/**
 * returns a pointer to the peptide_protein_association field of the peptide
 */
PeptideSrc* Peptide::getPeptideSrc() const {
  return !peptide_srcs_.empty() ? peptide_srcs_[0] : NULL;
}

/**
 *returns the pepide_srcs_
 */
vector<PeptideSrc*>& Peptide::getPeptideSrcVector() {
  return peptide_srcs_;
}

/**
 *Return the begining of the peptide_srcs_
 */
PeptideSrcIterator Peptide::getPeptideSrcBegin() {
  return peptide_srcs_.begin();
}

/**
 *Return the end of the peptide_srcs_
 */
PeptideSrcIterator Peptide::getPeptideSrcEnd() {
  return peptide_srcs_.end();
}

/**
 * \returns The number of peptide sources (i.e. proteins) the peptide has.
 */
int Peptide::getNumPeptideSrc(){
  return peptide_srcs_.size();
}

/**
 * get the peptide->first peptide_src->parent protein->database
 */
Database* Peptide::getFirstSrcDatabase() {
  return peptide_srcs_.at(0)->getParentProtein()->getDatabase();
}

// set by peptide_src?
/**
 * returns a pointer to the peptide's first parent protein field of the peptide
 */
Protein* Peptide::getParentProtein() {
  return peptide_srcs_.at(0)->getParentProtein();
}

/**
 * sets the sequence length of the peptide
 * length maximum of 255
 */
void Peptide::setLength(
  unsigned char length  ///< the length of sequence -in
  ) {
  length_ = length;
}

/* sequence-related getters and setters */
/**
 *\returns the sequence length of the peptide
 */
unsigned char Peptide::getLength() const {
  return length_;
}

/**
 * \brief Check whether a given sequence is equal to a given peptide.
 * \returns A Boolean indicating equality or not.
 */
static bool equal_peptides(
 char* peptide_sequence, ///< peptide sequence -in
 Peptide* peptide_object ///< peptide object -in
 )
{
  char* parent_sequence = 
    peptide_object->getPeptideSrc()->getParentProtein()->getSequencePointer();
  int start_idx = peptide_object->getPeptideSrc()->getStartIdx();

  int result = strncmp(peptide_sequence, 
                       &(parent_sequence[start_idx-1]), 
                       peptide_object->getLength());

  // Return true if strncmp returns 0.
  return((bool)(!result));
}

/**
 * \brief Get a string representation of the peptide sequence with no
 * added modification symbols.
 * Returns decoy sequence for decoy peptides.
 * \returns The newly-allocated sequence of peptide
 */
char* Peptide::getSequence() const {
  //Check that the peptide has parent(s)
  if (!sequence_.empty()) {
    return copy_string_part(sequence_.c_str(), sequence_.length());
  }

  if (decoy_modified_seq_ != NULL) {
    return modified_aa_to_unmodified_string(decoy_modified_seq_, length_);
  }
  string unshuffled = getUnshuffledSequence();
  if (unshuffled.empty()) {
    return NULL;
  }
  return copy_string_part(unshuffled.c_str(), unshuffled.length());
}

/**
 * \brief Get a string representation of the target (unshuffled)
 * peptide sequence with no added modification symbols.
 * For target peptides, returns the same as get_peptide_sequence.
 * \returns The newly-allocated sequence of peptide
 */
string Peptide::getUnshuffledSequence() const {
  if (!sequence_.empty()) {
    return sequence_;
  }
  //check parent(s) protein of the peptide not be empty 
  PeptideSrc* src = getPeptideSrc();
  if (src == NULL) {
    carp(CARP_ERROR, "Cannot get sequence from peptide with no peptide src.");
    return "";
  }
  char* parent_sequence = src->getParentProtein()->getSequencePointer(src->getStartIdx()-1);
  return string(parent_sequence, length_);
}

/**
 * \brief Get a pointer to the peptide sequence that is NOT null
 * terminated.
 * USE WITH CAUTION.  Pointer is to the parent protein sequence and
 * thus is not null-terminated until the end of the protein.  Parent
 * protein is taken from the first protein source.
 * 
 * \returns A pointer to an existing peptide sequence.
 */
char* Peptide::getSequencePointer() {
  if (peptide_srcs_.empty()) {
    carp(CARP_FATAL, "ERROR: no peptide_src to retrieve peptide sequence pointer\n");
  }
  PeptideSrc* src = getPeptideSrc();
  return src->getParentProtein()->getSequencePointer(src->getStartIdx()-1);
}

/**
 * \returns The sequence of peptide as used in sqt files, namely with
 * each flanking AA and any modifications 
 * 
 * Format is [AA|-].[peptide_sequence].[AA|-] where AA is a flanking
 * amino acid and - indicates this is the end of the protein sequence
 * Gets flanking AAs from the first peptide_src, thus must have at
 * least one peptide src 
 *
 * \returns A newly allocated char* with formated peptide sequence
 */
string Peptide::getSequenceSqt() {
  return string(1, getNTermFlankingAA()) + "." +
         getModifiedSequenceWithSymbols() +
         "." + string(1, getCTermFlankingAA());
}

/**
 * \brief Return a char for the amino acid n-terminal to the peptide
 * in the peptide src at the given index.
 *
 * \returns A char (A-Z) or - if peptide is the first in the protein.
 */
char Peptide::getNTermFlankingAA() {
  PeptideSrc* src = getPeptideSrc();
  if (src == NULL) {
    return '-';
  }

  // get protein seq
  Protein* protein = src->getParentProtein();

  // get peptide start idx, protein index starts at 1
  int start_index = src->getStartIdx();

  if (protein->isPostProcess()) {
    return ((PostProcessProtein*)protein)->getNTermFlankingAA(start_index - 1);
  }
  char* protein_seq = protein->getSequencePointer();

  char aa = '-';
  // if not at beginning, return char
  if( start_index > 1 ){
    aa = protein_seq[start_index - 2]; // -1 for 1-based shift
                                       // -1 for aa before start
  }
  return aa;
}

/**
 * \brief Return a char for the amino acid c-terminal to the peptide
 * in the peptide src at the given index.
 *
 * \returns A char (A-Z) or - if peptide is the last in the protein.
 */
char Peptide::getCTermFlankingAA() {
  PeptideSrc* src = getPeptideSrc();
  if (src == NULL) {
    return '-';
  }

  // get protein seq and length
  Protein* protein = src->getParentProtein();

  // get peptide start idx, protein index starts at 1
  int start_index = src->getStartIdx();

  if (protein->isPostProcess()) {
    return ((PostProcessProtein*)protein)->getCTermFlankingAA(start_index - 1);
  }
  char* protein_seq = protein->getSequencePointer();
  int protein_length = protein->getLength();

  // get peptide end idx
  int end_index = start_index + length_ - 1;

  char aa = '-';
  // if not at end, return char
  if( end_index < protein_length ){
    aa = protein_seq[end_index]; // -1 for 1-based shift, +1 for aa after end
  } 
  return aa;
}

void Peptide::addMod(const ModificationDefinition* mod, unsigned char index) {
  varMods_.push_back(Modification(mod, index));
}

void Peptide::setMods(const vector<Modification>& mods) {
  varMods_.clear();
  for (vector<Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
    if (!i->Static()) {
      varMods_.push_back(*i);
    }
  }
}

vector<Modification> Peptide::getMods() const {
  vector<Modification> mods = varMods_;
  vector<Modification> staticMods = getStaticMods();
  mods.insert(mods.end(), staticMods.begin(), staticMods.end());
  return mods;
}

vector<Modification> Peptide::getVarMods() const {
  return varMods_;
}

vector<Modification> Peptide::getStaticMods() const {
  char* seq = getSequence();
  vector<Modification> mods;
  for (unsigned char i = 0; seq[i] != '\0'; i++) {
    const vector<const ModificationDefinition*>& staticMods =
      ModificationDefinition::StaticMods(seq[i]);
    for (vector<const ModificationDefinition*>::const_iterator j = staticMods.begin();
         j != staticMods.end();
         j++) {
      if ((*j)->Position() == ANY ||
          ((*j)->Position() == PEPTIDE_N && i == 0) ||
          ((*j)->Position() == PEPTIDE_C && seq[i + 1] == '\0')) {
        mods.push_back(Modification(*j, (unsigned char)i));
      }
    }
  }
  free(seq);
  return mods;
}

bool Peptide::hasMonoLink() const {
  for (vector<Modification>::const_iterator iter = varMods_.begin(); iter != varMods_.end(); iter++) {
    if (iter->MonoLink()) {
      return true;
    }
  }
  return false;
}


/**
 * \brief Add a modification to a peptide.
 *
 * Adds the modified sequence to the peptide and changes the peptide
 * mass based on the mass change in the peptide_mod.
 * \returns void
 */
void Peptide::setMod(
  MODIFIED_AA_T* mod_seq, ///< modified seq to add
  PEPTIDE_MOD_T* pep_mod  ///< mod that made the seq
) {
  if (!mod_seq || !pep_mod) {
    carp(CARP_ERROR, "Cannot modify peptide. mod, or seq is NULL.");
    return;
  }
  vector<Modification> newMods;
  Modification::FromSeq(mod_seq, length_, &sequence_, &newMods);
  for (vector<Modification>::const_iterator i = newMods.begin(); i != newMods.end(); i++) {
    varMods_.push_back(*i);
  }
}

string Peptide::getModsString() const {
  vector<Modification> allMods = getStaticMods();
  for (vector<Modification>::const_iterator i = varMods_.begin(); i != varMods_.end(); i++) {
    allMods.push_back(*i);
  }
  std::sort(allMods.begin(), allMods.end(), Modification::SortFunction);
  vector<string> modStrings;
  for (vector<Modification>::const_iterator i = allMods.begin(); i != allMods.end(); i++) {
    modStrings.push_back(i->String());
  }
  return StringUtils::Join(modStrings, ',');
}

bool Peptide::isModified() {
  return !varMods_.empty();
}

bool Peptide::isDecoy() {
  return decoy_modified_seq_ != NULL;
}

std::string Peptide::getDecoyType() {
  string returnValue = "target";
  if (decoy_modified_seq_ != NULL) {
    returnValue = "decoy";
  }
  return(returnValue);
}

/**
 * \brief Get the modified peptide sequence
 *
 * If the peptide has no modifications, create a sequence of
 * MODIFIED_AA_T's in which none of them are actually modified.
 * \returns A newly allocated copy of the sequence of MODIFIED_AA_Ts.
 */
MODIFIED_AA_T* Peptide::getModifiedAASequence() {
  MODIFIED_AA_T* seq_copy = NULL;
  char* seq = getSequence();
  if (!varMods_.empty()) {
    string seqString(seq);
    seq_copy = Modification::ToSeq(seqString, varMods_);
  } else {
    convert_to_mod_aa_seq(seq, &seq_copy);
  }
  std::free(seq);
  return seq_copy;
}

string Peptide::unmodifySequence(const string& seq) {
  string outSeq;
  outSeq.reserve(seq.length());

  size_t i = 0, end = seq.length();
  bool flanks = seq.length() >= 5 && seq[1] == '.' && seq[seq.length() - 2] == '.';
  for (size_t i = 0; i < seq.length(); i++) {
    char c = seq[i];
    if ((flanks && (i < 2 || i >= seq.length() - 2)) || ('A' <= c && c <= 'Z')) {
      outSeq.push_back(c);
      continue;
    } else if (c == '[') {
      // Bracket mod
      size_t j = i;
      do {
        if (++j >= seq.length()) {
          throw runtime_error("Unclosed '[' in sequence '" + seq + "'");
        }
      } while (seq[j] != ']');
      i = j;
    } else {
      // Symbol mod
    }
  }
  outSeq.reserve(outSeq.length());
  return outSeq;
}

void Peptide::setUnmodifiedSequence(const string& sequence) {
  sequence_ = sequence;
  setLength(sequence.length());
}

/**
 * sets the modified sequence for the peptide
 */
void Peptide::setModifiedAASequence(
  MODIFIED_AA_T* mod_seq, ///< modified sequence to set
  bool decoy ///< is the peptide a decoy?
) {
  Modification::FromSeq(mod_seq, length_, &sequence_, &varMods_);
  if (decoy) {
    if (decoy_modified_seq_) {
      std::free(decoy_modified_seq_);
    }
    decoy_modified_seq_ = copy_mod_aa_seq(mod_seq);
  }
}

void Peptide::setDecoyModifiedSeq(MODIFIED_AA_T* decoy_modified_seq) {
  decoy_modified_seq_ = decoy_modified_seq;
}

/**
 * \brief Get the modified aa sequence in string form.
 *
 * If the peptide has no modifications, returns same string as
 * get_peptide_sequence.  If modified, adds the mod symbols to the string.
 * \returns The peptide sequence including any modifications.
 */
string Peptide::getModifiedSequenceWithSymbols() {
  if (decoy_modified_seq_) {
    char* seqTmp = 
      modified_aa_string_to_string_with_symbols(decoy_modified_seq_, length_);
    string seqString(seqTmp);
    free(seqTmp);
    return seqString;
  }

  char* seqC = getSequence();
  string seq(seqC);
  free(seqC);
  map<unsigned char, string> symbols;
  for (vector<Modification>::const_iterator i = varMods_.begin(); i != varMods_.end(); i++) {
    map<unsigned char, string>::const_iterator j = symbols.find(i->Index());
    if (j == symbols.end()) {
      symbols[i->Index()] = i->Symbol();
    } else {
      symbols[i->Index()] += i->Symbol();
    }
  }

  for (map<unsigned char, string>::const_reverse_iterator i = symbols.rbegin();
       i != symbols.rend();
       i++) {
    seq.insert(i->first + 1, i->second);
  }

  return seq;
}

/**
 * \brief Get the modified aa sequence in string form.
 *
 * If the peptide has no modifications, returns same string as
 * get_peptide_sequence.  If modified, adds in brackets the masses of
 * all modifications.  If merge_masses is true, prints the sum of all
 * modifications for a residue.  If false, prints all masses in a
 * comma separated list.
 * \returns The peptide sequence including any modifications.
 */
string Peptide::getModifiedSequenceWithMasses() {
  if (decoy_modified_seq_) {
    char* seqTmp = modified_aa_string_to_string_with_masses(
      decoy_modified_seq_, length_, GlobalParams::getModMassFormat());
    string seqString(seqTmp);
    free(seqTmp);
    return seqString;
  }

  char* seqC = getSequence();
  string seq(seqC);
  free(seqC);
  map< int, vector<double> > masses;
  for (vector<Modification>::const_iterator i = varMods_.begin(); i != varMods_.end(); i++) {
    map< int, vector<double> >::const_iterator j = masses.find(i->Index());
    if (j == masses.end()) {
      masses[i->Index()] = vector<double>(1, i->DeltaMass());
    } else {
      masses[i->Index()].push_back(i->DeltaMass());
    }
  }
  int precision = GlobalParams::getModPrecision();
  MASS_TYPE_T massType = GlobalParams::getIsotopicMass();
  for (map< int, vector<double> >::const_reverse_iterator i = masses.rbegin();
       i != masses.rend();
       i++) {
    char buffer[64];
    switch (GlobalParams::getModMassFormat()) {
      case MOD_MASS_ONLY: {
        double sum = accumulate(i->second.begin(), i->second.end(), 0.0);
        sprintf(buffer, "[%.*f]", precision, sum);
        break;
      }
      case AA_PLUS_MOD: {
        double sum = accumulate(i->second.begin(), i->second.end(),
                                get_mass_amino_acid(seq[i->first], massType));
        sprintf(buffer, "[%.*f]", precision, sum);
        break;
      }
      case MOD_MASSES_SEPARATE: {
        vector<string> massStrings;
        for (vector<double>::const_iterator j = i->second.begin(); j != i->second.end(); j++) {
          massStrings.push_back(StringUtils::ToString(*j, precision));
        }
        sprintf(buffer, "[%s]", StringUtils::Join(massStrings, ',').c_str());
        break;
      }
    }
    seq.insert(i->first + 1, buffer);
  }

  return seq;
}

/* getters requiring calculation */

/**
 * \brief Count the number of modified amino acids in the
 * peptide. This number is distnct from the number of aamods in the
 * peptide mod since one amino acid can have more than one
 * modification on it.
 * \returns The number of amino acids in the peptide that have at
 * least one modification.
 */
int Peptide::countModifiedAAs(){
  set<unsigned char> indices;
  for (vector<Modification>::const_iterator i = varMods_.begin(); i != varMods_.end(); i++) {
    indices.insert(i->Index());
  }
  return indices.size();
}

/**
 * \returns The mass of the given peptide as determined by the aa sequence.
 */
FLOAT_T Peptide::calcSequenceMass(
  const string& peptide, ///< the query peptide -in
  MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  ) {
  FLOAT_T peptide_mass = 0;
  for (string::const_iterator i = peptide.begin(); i != peptide.end(); i++) {
    peptide_mass += get_mass_amino_acid(*i, mass_type);
  }
  if(mass_type == AVERAGE){
    return peptide_mass + MASS_H2O_AVERAGE;
  }
  return peptide_mass + MASS_H2O_MONO;
}

/**
 * This appears to be the same as calc_sequence_mass??
 * \returns The mass of the given peptide.
 */
FLOAT_T Peptide::calcMass(MASS_TYPE_T mass_type) const {
  FLOAT_T mass = 0;
  char* seq = getSequence();
  for (char* i = seq; *i != '\0'; i++) {
    mass += AminoAcidUtil::GetMass(*i, mass_type != AVERAGE);
  }
  free(seq);

  return mass_type != AVERAGE ? mass + MASS_H2O_MONO : mass + MASS_H2O_AVERAGE;
}

FLOAT_T Peptide::calcModifiedMass(MASS_TYPE_T mass_type) const {
  FLOAT_T mass = calcMass(mass_type);
  for (vector<Modification>::const_iterator i = varMods_.begin(); i != varMods_.end(); i++) {
    mass += i->DeltaMass();
  }
  vector<Modification> staticMods = getStaticMods();
  for (vector<Modification>::const_iterator i = staticMods.begin(); i != staticMods.end(); i++) {
    mass += i->DeltaMass();
  }
  return mass;
}

FLOAT_T Peptide::calcModifiedMass() const {
  return calcModifiedMass(GlobalParams::getIsotopicMass());
}

/**
 * Examines the peptide sequence and counts how many tryptic missed
 * cleavage sites exist. 
 *\returns the number of missed cleavage sites in the peptide
 */
int Peptide::getMissedCleavageSites() {

  int missed_count = 0;
  int aa_idx = 0;
  char* sequence = getSequencePointer();

  // count the missed cleavage sites
  for(; aa_idx < length_-1; ++aa_idx){
    if(sequence[aa_idx] == 'K' ||
       sequence[aa_idx] == 'R'){
      
      // skip one that are followed by a P
      if(sequence[aa_idx+1] == 'P'){
        continue;
      }
      else{
        ++missed_count;
      }      
    } 
  }
  
  return missed_count;
}

int Peptide::getMissedCleavageSites(
  set<int> skip ///< skip these amino acid indices.
) {
  int missed_count = 0;
  char* sequence = getSequencePointer();

  // count the missed cleavage sites
  for (int aa_idx = 0; aa_idx < length_ - 1; ++aa_idx) {
    if (skip.find(aa_idx) != skip.end()) {
      continue;
    }

    bool cleavage_prevented = true;
    if (sequence[aa_idx] == 'K' || sequence[aa_idx] == 'R') {
      cleavage_prevented = false;
      // skip one that are followed by a P
      if (sequence[aa_idx + 1] == 'P') {
        cleavage_prevented = true;
        continue;
      }

      for (vector<Modification>::const_iterator modIter = varMods_.begin();
           modIter != varMods_.end();
           modIter++) {
        if (modIter->Index() == aa_idx && modIter->PreventsCleavage()) {
          cleavage_prevented = true;
          break;
        }
      }
    }

    if (!cleavage_prevented) {
      ++missed_count;
    }
  }

  return missed_count;
}


/**
 * \brief Find the distance from the c-terminus of the source protein
 * to the c-terminus of the peptide (seq[0]).  
 * In the case of multiple source proteins, return the smallest
 * distance.
 * \returns The distance from the protein c-terminus.
 */
int Peptide::getCDistance(){

  int min_index = MAX_PROTEIN_SEQ_LENGTH;
  for(vector<PeptideSrc*>::iterator iter= peptide_srcs_.begin();
       iter!= peptide_srcs_.end();
       ++iter
      ){
    PeptideSrc* cur_src = *iter;
    int index = cur_src->getStartIdx();
    if( index < min_index ){
      min_index = index;
    }
  }
  return min_index - 1;
}

/**
 * \brief Find the distance from the n-terminus of the source protein
 * to the n-terminus of the peptide.
 * In the case of multiple source proteins, return the smallest
 * distance.
 * \returns The distance from the protein c-terminus.
 */
int Peptide::getNDistance(){

  int min_index = MAX_PROTEIN_SEQ_LENGTH;
  int peptide_length = getLength();

  for(vector<PeptideSrc*>::iterator iter= peptide_srcs_.begin();
       iter!= peptide_srcs_.end();
       ++iter
      ){
    PeptideSrc* cur_src = *iter;
    int protein_length = cur_src->getParentProtein()->getLength();
    // get index of end
    int start_index = cur_src->getStartIdx();

    int cidx = protein_length - (start_index + peptide_length - 1);
    if( cidx < min_index){
      min_index = cidx;
  }
  }  
  return min_index;

}

/**
 * Change the given target peptide into a decoy by randomizing its sequence.
 * Uses settings in parameter.c to decide between shuffling and
 * reversing the sequence.  Any modifications that exist will be
 * maintained on the same amino acids whose position will move.  If
 * the peptide is already a decoy, replaces the existing decoy
 * sequence. 
 */
void Peptide::transformToDecoy() {
  bool reverse_seq = false; // For now, only generate shuffled decoy

  // delete any existing decoy sequence
  if (decoy_modified_seq_){  
    freeModSeq(decoy_modified_seq_); 
  }
  // if the peptide is already modified, shuffle the modified sequence
  if (!varMods_.empty()) {
    MODIFIED_AA_T* new_seq = NULL;
    if (reverse_seq) {
      new_seq = generateReversedModSequence();
    } else {
      new_seq = generateShuffledModSequence();
    }
    decoy_modified_seq_ = new_seq;
  } else {// shuffle the unmodified sequence
    char* new_seq = NULL;
    if (reverse_seq) {
      new_seq = generateReversedSequence();
    } else {
      new_seq = generateShuffledSequence();
    }
    convert_to_mod_aa_seq(new_seq, &(decoy_modified_seq_));
    std::free(new_seq);
  }
}

/**
 * \brief Return a randomly shuffled version of the given peptide's
 * sequence, leaving the terminal amino acids in place.  Ensures that
 * the shuffled version is not the same as the given peptide.
 * 
 * \returns A newly-allocated char array with the shuffled sequence.
 */
static const int MAX_SHUFFLES = 5; // Don't bother trying to shuffle more than this.
char* Peptide::generateShuffledSequence() {
  // Allocate a copy of the peptide.
  char* sequence = getSequence();
  int length = length_;

  // Shuffle from left to right, using the Knuth algorithm for shuffling.
  int num_shuffles = 0;
  do {

    // Don't move the n-term and c-term amino acids
    int start_idx = 1;
    int end_idx = length - 2;

    while(start_idx < end_idx){
      int switch_idx = get_random_number_interval(start_idx, end_idx);
      char temp_char = sequence[start_idx];
      sequence[start_idx] = sequence[switch_idx];
      sequence[switch_idx] = temp_char;
      ++start_idx;
    }
    num_shuffles++;
  } while (equal_peptides(sequence, this) && (num_shuffles < MAX_SHUFFLES));

  return sequence;
}

/**
 * \brief Return a reversed version of the given peptide's sequence as
 * an array of char (A-Z).  Leave the first and last residue
 * unchanged.  If the reversed sequence is identical to the target,
 * shuffle the sequence instead.
 *
 * \returns A newly-allocated char array of the reversed sequence.
 */
char* Peptide::generateReversedSequence() {
  char* sequence = getSequence();
  int length = length_;
  int start_idx = 1;       // leave first ...
  int end_idx = length -2; // ...and last residue in place
  char temp_char = 0;

  while(start_idx < end_idx){
    temp_char = sequence[end_idx];
    sequence[end_idx] = sequence[start_idx];
    sequence[start_idx] = temp_char;
    start_idx++;
    end_idx--;
  }

  // check to see if the reversed sequence is the same as original
  if( strncmp(sequence, getSequencePointer(), length_) == 0 ){
    carp(CARP_DETAILED_INFO, 
         "Peptide %s is a palindrome and will be shuffled instead of reversed.",
         sequence);
    std::free(sequence);
    sequence = generateShuffledSequence();
  }

  return sequence;
}

/**
 * \brief Return a randomly shuffled version of the given peptide's 
 * sequence as an array of MODIIFIED_AA_T.  Based on the peptide type,
 * will leave the end(s) unchanged to preserve the tryptic property.
 * 
 *\returns A newly-allcoated MODIFIED_AA_T array of the shuffled sequence.
 */
MODIFIED_AA_T* Peptide::generateShuffledModSequence() {
  //TODO (BF 6-Apr-09): should we warn if seq is len 3 and won't change?
  MODIFIED_AA_T* sequence = getModifiedAASequence();
  int length = length_;
  int start_idx = 1;
  int end_idx = length - 2;
  int switch_idx = 0;
  MODIFIED_AA_T temp_aa = 0;

  // shuffle from left to right, using the Knuth algorithm for shuffling.
  while (start_idx < end_idx) {
    switch_idx = get_random_number_interval(start_idx, end_idx);
    temp_aa = sequence[start_idx];
    sequence[start_idx] = sequence[switch_idx];
    sequence[switch_idx] = temp_aa;
    ++start_idx;
  }

  return sequence;
}

/**
 * \brief Return a reversed version of the given peptide's sequence as
 * an array of MODIFIED_AA_T.  Leave the first and last residue
 * unchanged.  If the reversed sequence is identical to the target,
 * shuffle the sequence instead.
 *
 * \returns A newly-allocated MODIFIED_AA_T array of the reversed sequence.
 */
MODIFIED_AA_T* Peptide::generateReversedModSequence() {

  MODIFIED_AA_T* sequence = getModifiedAASequence();
  int length = length_;
  int start_idx = 0;
  int end_idx = length - 1;
  MODIFIED_AA_T temp_aa = 0;

  // first check to see if it will yield a different seq when reversed
  if( modified_aa_seq_is_palindrome(sequence, length) == true){
    return generateShuffledModSequence();
  }

  // Do not move the first and last residue, regardless of enzyme
  ++start_idx;
  --end_idx;

  // reverse  
  while(start_idx < end_idx){
    temp_aa = sequence[start_idx];
    sequence[start_idx] = sequence[end_idx];
    sequence[end_idx] = temp_aa;
    ++start_idx;
    end_idx--;
  }

  return sequence;
}

string Peptide::getId() const {
  char* seqChars = getSequence();
  string seq(seqChars);
  free(seqChars);
  return getId(seq, varMods_);
}

string Peptide::getId(const string& unmodifiedSeq, const vector<Modification>& mods) {
  string id = unmodifiedSeq;
  vector<Modification> modsSorted(mods);
  std::sort(modsSorted.begin(), modsSorted.end(), Modification::SortFunction);
  for (vector<Modification>::const_iterator i = modsSorted.begin(); i != modsSorted.end(); i++) {
    char buf[64];
    sprintf(buf, " %x %d", i->Definition(), i->Index());
    id += buf;
  }
  return id;
}

/* Comparisons for sorting */

/**
 * Compare two peptide sequences.
 * \returns Zero (0) if the sequences are identical, -1 if the first
 * sequence is less than the first and 1 if the first sequence is
 * greater than teh first.
 */
int Peptide::triCompareSequence(
  Peptide* peptide_one,  ///< the peptide sequence to compare  -out
  Peptide* peptide_two  ///< the peptide sequence to compare  -out
  )
{
  // find the shorter peptide
  int short_len = 0;
  if( peptide_one->length_ < peptide_two->length_ ){
    short_len = peptide_one->length_;
  } else {
    short_len = peptide_two->length_;
  }

  char* seq_one = peptide_one->getPeptideSrc()->getSequencePointer();
  char* seq_two = peptide_two->getPeptideSrc()->getSequencePointer();
    
  // stop comparing as soon as they differ
  int pep_idx = 0;
  for(pep_idx = 0; pep_idx < short_len; pep_idx++ ){
      if(seq_one[pep_idx] != seq_two[pep_idx]){
        pep_idx++; // stop pointing one after the difference
        break;
      }
  }

  // move index back one to the last letter compared
  pep_idx--;

  // if the seqs are the same up to this point, then compare the lengths
  if( seq_one[pep_idx] == seq_two[pep_idx] ){
    if( peptide_one->length_ == peptide_two->length_ ){ // same seq
      return 0;
    } else if ( peptide_one->length_ < peptide_two->length_ ){ 
      return -1;
    } else {
      return 1;
    }
  }

  // else, the seqs are different
  if(seq_one[pep_idx] < seq_two[pep_idx]){
    return -1;
  } else {
    return 1;
  }

}
/**
 * Compare the sequence of two peptides and return true if the first
 * petpide sequence is less than (in a lexical sort) the second
 * peptide.  Return false if they are idential peptides.
 */
bool Peptide::lessThan(
  Peptide* peptide_one,
  Peptide* peptide_two
  ){
  // find the shorter peptide
  int short_len = 0;
  if( peptide_one->length_ < peptide_two->length_ ){
    short_len = peptide_one->length_;
  } else {
    short_len = peptide_two->length_;
  }

  char* seq_one = peptide_one->getPeptideSrc()->getSequencePointer();
  char* seq_two = peptide_two->getPeptideSrc()->getSequencePointer();
    
  // stop comparing as soon as they differ
  int pep_idx = 0;
  for(pep_idx = 0; pep_idx < short_len; pep_idx++ ){
      if(seq_one[pep_idx] != seq_two[pep_idx]){
        break;
      }
  }
  if( seq_one[pep_idx] == seq_two[pep_idx] ){
    return (peptide_one->length_ < peptide_two->length_);
  } 

  return (seq_one[pep_idx] < seq_two[pep_idx]);
}

/* Public functions--Printing / parsing */

/**
 * Fills the given vectors with the names and descriptions of all
 * proteins containing this peptide.  Makes the descriptions
 * xml-friendly by swapping double for single quotes and angled braces
 * for square. Returned in the same order as getFlankingAAs().  Clears
 * any existing values in the vectors.
 * Adapted from Match::get_information_of_proteins()
 * \returns The number of proteins.
 */
int Peptide::getProteinInfo(vector<string>& protein_ids,
                            vector<string>& protein_descriptions){

  protein_ids.clear();
  protein_descriptions.clear();

  for(PeptideSrcIterator iter = getPeptideSrcBegin(); 
      iter!=getPeptideSrcEnd();
      ++iter
   ){
    PeptideSrc* peptide_src = *iter; 
    Protein* protein = peptide_src->getParentProtein();
    protein_ids.push_back(protein->getIdPointer());

    string description = "";
    const string& annotation_pointer = protein->getAnnotationPointer();
    if (!annotation_pointer.empty()) {
      description = protein->getAnnotationPointer();
      // replace double quotes with single quotes
      replace(description.begin(), description.end(), '"', '\'');
      // remove any xml tags in the description by replacing <> with []
      replace(description.begin(), description.end(), '<', '[');
      replace(description.begin(), description.end(), '>', ']');
    }
    protein_descriptions.push_back(description);

  } 

  return protein_ids.size();
}

/* Public functions--Iterators */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(
  Peptide* peptide ///< peptide sequence to iterate -in
  )
{
  RESIDUE_ITERATOR_T* residue_iterator =
    (RESIDUE_ITERATOR_T*)mycalloc(1, sizeof(RESIDUE_ITERATOR_T));
  
  residue_iterator->peptide =  peptide;
  residue_iterator->residue_idx = 0;
  residue_iterator->sequence = peptide->getSequence();
  return residue_iterator;
}        

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(
  RESIDUE_ITERATOR_T* residue_iterator ///< free this object -in
  )
{
  free(residue_iterator->sequence);
  free(residue_iterator);
}

/**
 * The basic iterator functions.
 * \returns true if there are additional residues to iterate over, false if not.
 */
bool residue_iterator_has_next(
  RESIDUE_ITERATOR_T* residue_iterator ///< the query iterator -in
  )
{
  return (residue_iterator->residue_idx < residue_iterator->peptide->getLength());
}

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(
  RESIDUE_ITERATOR_T* residue_iterator  ///< the query iterator -in
  )
{
  ++residue_iterator->residue_idx;
  return residue_iterator->sequence[residue_iterator->residue_idx - 1];
}

/**
 * \brief Builds a comma delimited string listing the 
 * protein id(peptide start index) for the sources of 
 * a peptide
 *
 * \returns a string of the protein sources for this peptide
 */
string Peptide::getProteinIdsLocations() {

  set<string> protein_ids_locations;
  string protein_field_string;  

  std::ostringstream protein_field_stream;
  if (!peptide_srcs_.empty()) {
    for( PeptideSrcIterator iter = getPeptideSrcBegin();
      iter!=getPeptideSrcEnd();++iter){
      
      PeptideSrc* peptide_src =*iter;
      Protein* protein = peptide_src->getParentProtein();
      string& protein_id = protein->getIdPointer();
      std::ostringstream protein_loc_stream;
      protein_loc_stream << protein_id;

      if (!protein->isPostProcess()) {
       int peptide_loc = peptide_src->getStartIdx();      
        protein_loc_stream << "(" << peptide_loc << ")";
      } else if (peptide_src->getStartIdxOriginal() > 0) {
        protein_loc_stream << "(" << peptide_src->getStartIdxOriginal() << ")";
      }

      protein_ids_locations.insert(protein_loc_stream.str());
    }
  }

  set<string>::iterator result_iter = protein_ids_locations.begin();
  protein_field_string = *result_iter;

  while(++result_iter != protein_ids_locations.end()) {
    protein_field_string += "," + *result_iter;
  }

  return protein_field_string;
}

 
/**
 * \brief Builds a comma delimited string listing the protein ids
 * for the sources of a peptide.
 */
vector<string> Peptide::getProteinIds() {
  vector<string> ids;
  for (PeptideSrcIterator iter = getPeptideSrcBegin(); iter != getPeptideSrcEnd(); ++iter) {
    ids.push_back((*iter)->getParentProtein()->getIdPointer());
  }
  return ids;
}

/**
 * \brief Builds a comma delimited string listing the flanking amino acids
 * for the sources of a peptide.
 *
 * \returns a pointer to the string. Caller is responsible for freeing memeory.
 * If peptide has no sources returns NULL.
 */
char* Peptide::getFlankingAAs() {

  string flanking_string = "";

  // iterate over all PeptideSrc
  for (PeptideSrcIterator iter = getPeptideSrcBegin();
       iter != getPeptideSrcEnd();
       ++iter) {
    // add comma if there is already flanking AA(s) in the string
    if (!flanking_string.empty()) {
      flanking_string += ',';
    }
    PeptideSrc* src = *iter;
    Protein* protein = src->getParentProtein();
    int start_idx = src->getStartIdx();
    if (!protein->isPostProcess()) {
      // not post process, get flanking AAs from protein sequence
      int end_idx = start_idx + getLength() - 1;
      int protein_length = protein->getLength();
      char* protein_seq = protein->getSequencePointer();
      flanking_string += (start_idx > 1) ? protein_seq[start_idx - 2] : '-';
      flanking_string += (end_idx < protein_length) ? protein_seq[end_idx] : '-';
    } else {
      // post process, get flanking AAs from the protein object
      PostProcessProtein* post_process_protein = (PostProcessProtein*)protein;
      flanking_string += post_process_protein->getNTermFlankingAA(start_idx - 1);
      flanking_string += post_process_protein->getCTermFlankingAA(start_idx - 1);
    }
  }

  // convert to c string
  int length = flanking_string.length();
  char* flanking_field = (char*)malloc(sizeof(char)*(length + 1));
  memcpy(flanking_field, flanking_string.c_str(), length);
  flanking_field[length] = '\0';

  return flanking_field;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

