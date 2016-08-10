#include "Modification.h"
#include "io/carp.h"
#include "util/AminoAcidUtil.h"
#include "util/MathUtil.h"
#include "util/Params.h"
#include "util/StringUtils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace Crux;

ModificationDefinition::ModificationDefinition(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool preventsCleavage, bool preventsXLink, char symbol)
  : deltaMass_(deltaMass), position_(position),
    symbol_(symbol), preventsCleavage_(preventsCleavage), preventsXLink_(preventsXLink) {
  AddAminoAcids(aminoAcids);
}

ModificationDefinition::~ModificationDefinition() {
}

const ModificationDefinition* ModificationDefinition::New(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool isStatic, bool preventsCleavage, bool preventsXLink) {
  return isStatic
    ? NewStaticMod(aminoAcids, deltaMass, position, preventsCleavage, preventsXLink)
    : NewVarMod(aminoAcids, deltaMass, position, preventsCleavage, preventsXLink, '\0');
}

const ModificationDefinition* ModificationDefinition::NewStaticMod(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool preventsCleavage, bool preventsXLink) {
  // Look for existing
  vector<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  for (vector<ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    if (MathUtil::AlmostEqual((*i)->deltaMass_, deltaMass, Params::GetInt("mod-precision")) &&
        (position == UNKNOWN || (*i)->position_ == UNKNOWN || (*i)->position_ == position) &&
        (*i)->preventsCleavage_ == preventsCleavage &&
        (*i)->preventsXLink_ == preventsXLink) {
      string added = (*i)->AddAminoAcids(aminoAcids);
      for (string::const_iterator j = added.begin(); j != added.end(); j++) {
        modContainer_.staticMods_[*j].push_back(*i);
      }
      if (position != UNKNOWN && (*i)->position_ == UNKNOWN) {
        (*i)->position_ = position;
      }
      return *i;
    }
  }
  // Does not exist yet
  ModificationDefinition* mod = new ModificationDefinition(
    aminoAcids, deltaMass, position, '\0', preventsCleavage, preventsXLink);
  modContainer_.Add(mod);
  return mod;
}

const ModificationDefinition* ModificationDefinition::NewVarMod(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool preventsCleavage, bool preventsXLink, char symbol) {
  // Look for existing
  for (vector<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    if (MathUtil::AlmostEqual((*i)->deltaMass_, deltaMass, Params::GetInt("mod-precision")) &&
        (position == UNKNOWN || (*i)->position_ == UNKNOWN || (*i)->position_ == position) &&
        (*i)->preventsCleavage_ == preventsCleavage &&
        (*i)->preventsXLink_ == preventsXLink &&
        (symbol == '\0' || (*i)->symbol_ == symbol)) {
      (*i)->AddAminoAcids(aminoAcids);
      if (position != UNKNOWN && (*i)->position_ == UNKNOWN) {
        (*i)->position_ = position;
      }
      return *i;
    }
  }
  // Does not exist yet
  if (symbol == '\0') {
    symbol = modContainer_.NextSymbol();
  } else {
    modContainer_.ConsumeSymbol(symbol);
  }
  ModificationDefinition* mod = new ModificationDefinition(
    aminoAcids, deltaMass, position, preventsCleavage, preventsXLink, symbol);
  modContainer_.Add(mod);
  return mod;
}

string ModificationDefinition::String() const {
  // For debugging
  stringstream ss;
  ss << '[' << this << ']';
  if (!aminoAcids_.empty()) {
    ss << '[' << StringUtils::Join(aminoAcids_) << ']';
  } else {
    ss << "[?]";
  }
  ss << '[' << DeltaMass() << ']';
  if (Static()) {
    ss << "[static]";
  } else {
    ss << "[variable " << Symbol() << ']';
  }
  ss << '[';
  switch (Position()) {
    default:        ss << "unknown"; break;
    case ANY:       ss << "any"; break;
    case PEPTIDE_N: ss << "peptide N"; break;
    case PEPTIDE_C: ss << "peptide C"; break;
    case PROTEIN_N: ss << "protein N"; break;
    case PROTEIN_C: ss << "protein C"; break;
  }
  ss << ']';
  return ss.str();
}

void ModificationDefinition::ListAll() {
  // For debugging
  ListStaticMods();
  ListVarMods();
}

void ModificationDefinition::ListStaticMods() {
  // For debugging
  carp(CARP_INFO, "Listing static modifications");
  vector<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  for (vector<ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    carp(CARP_INFO, "%s", (*i)->String().c_str());
  }
}

void ModificationDefinition::ListVarMods() {
  // For debugging
  carp(CARP_INFO, "Listing variable modifications");
  for (vector<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    carp(CARP_INFO, "%s", (*i)->String().c_str());
  }
}

void ModificationDefinition::ClearAll() {
  modContainer_ = ModificationDefinitionContainer();
}

void ModificationDefinition::ClearStaticMods() {
  vector<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  for (vector<ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    delete *i;
  }
  modContainer_.staticMods_.clear();
}

void ModificationDefinition::ClearVarMods() {
  for (vector<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    delete *i;
  }
  modContainer_.varMods_.clear();
  modContainer_.InitSymbolPool();
}

vector<const ModificationDefinition*> ModificationDefinition::AllMods() {
  vector<const ModificationDefinition*> mods = StaticMods();
  vector<const ModificationDefinition*> varMods = VarMods();
  mods.insert(mods.end(), varMods.begin(), varMods.end());
  return mods;
}

vector<const ModificationDefinition*> ModificationDefinition::StaticMods(char c) {
  vector<ModificationDefinition*> mods = modContainer_.StaticMods(c);
  vector<const ModificationDefinition*> modsConst;
  for (vector<ModificationDefinition*>::const_iterator i = mods.begin();
       i != mods.end();
       i++) {
    modsConst.push_back(*i);
  }
  return modsConst;
}

vector<const ModificationDefinition*> ModificationDefinition::VarMods() {
  vector<const ModificationDefinition*> mods;
  for (vector<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    mods.push_back(*i);
  }
  return mods;
}

double ModificationDefinition::DeltaMass(char c, ModPosition position) {
  double mass = 0.0;
  vector<const ModificationDefinition*> staticMods = StaticMods(c);
  for (vector<const ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    if ((*i)->position_ == UNKNOWN || (*i)->position_ == ANY || position == (*i)->position_) {
      mass += (*i)->DeltaMass();
    }
  }
  return mass;
}

const set<char>& ModificationDefinition::AminoAcids() const {
  return aminoAcids_;
}

double ModificationDefinition::DeltaMass() const {
  return deltaMass_;
}

bool ModificationDefinition::Static() const {
  return symbol_ == '\0';
}

ModPosition ModificationDefinition::Position() const {
  return position_;
}

char ModificationDefinition::Symbol() const {
  return symbol_;
}

bool ModificationDefinition::PreventsCleavage() const {
  return preventsCleavage_;
}

bool ModificationDefinition::PreventsXLink() const {
  return preventsXLink_;
}

const ModificationDefinition* ModificationDefinition::Find(char symbol) {
  for (vector<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    if ((*i)->symbol_ == symbol) {
      return *i;
    }
  }
  return NULL;
}

const ModificationDefinition* ModificationDefinition::Find(
  double deltaMass, bool isStatic, ModPosition position) {
  vector<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  vector<ModificationDefinition*>* mods = isStatic ? &staticMods : &modContainer_.varMods_;
  for (vector<ModificationDefinition*>::const_iterator i = mods->begin();
       i != mods->end();
       i++) {
    if ((position == UNKNOWN || position == (*i)->Position()) &&
        MathUtil::AlmostEqual((*i)->deltaMass_, deltaMass, Params::GetInt("mod-precision"))) {
      return *i;
    }
  }
  return NULL;
}

string ModificationDefinition::AddAminoAcids(const string& aminoAcids) {
  string added;
  for (string::const_iterator i = aminoAcids.begin(); i != aminoAcids.end(); i++) {
    if (*i == 'X') {
      for (char j = 'A'; j <= 'Z'; j++) {
        if (aminoAcids_.insert(j).second) {
          added.push_back(j);
        }
      }
      break;
    }
    if (aminoAcids_.insert(*i).second) {
      added.push_back(*i);
    }
  }
  return added;
}

void swap(Modification& x, Modification& y) {
  using std::swap;
  swap(x.index_, y.index_);
  swap(x.mod_, y.mod_);
}

ostream& operator<<(ostream& stream, const Modification& mod) {
  stream << mod.String();
  return stream;
}

Modification::Modification(const ModificationDefinition* mod, unsigned char index)
  : index_(index), mod_(mod) {
}

Modification::Modification(const Modification& other)
  : index_(other.index_), mod_(other.mod_) {
}

Modification::~Modification() {
}

Modification& Modification::operator=(Modification other) {
  swap(*this, other);
  return *this;
}

bool Modification::operator==(const Modification& other) const {
  return index_ == other.index_ && mod_ == other.mod_;
}

bool Modification::operator!=(const Modification& other) const {
  return !(*this == other);
}

string Modification::String() const {
  char buffer[32];
  string positionStr;
  switch (mod_->Position()) {
    case PEPTIDE_N: positionStr = "_n"; break;
    case PEPTIDE_C: positionStr = "_c"; break;
    case PROTEIN_N: positionStr = "_N"; break;
    case PROTEIN_C: positionStr = "_C"; break;
  }
  sprintf(buffer, "%d_%c_%.*f%s",
          index_ + 1,
          mod_->Static() ? 'S' : 'V',
          Params::GetInt("mod-precision"), mod_->DeltaMass(),
          positionStr.c_str());
  return string(buffer);
}

Modification Modification::Parse(const string& modString, Peptide* peptide) {
  vector<string> pieces = StringUtils::Split(modString, '_');
  if (3 > pieces.size() || pieces.size() > 4) {
    throw runtime_error("Could not parse modification string '" + modString + "'");
  }
  unsigned char index = (unsigned char)StringUtils::FromString<unsigned int>(pieces[0]) - 1;
  bool isStatic = false;
  if (StringUtils::IEquals(pieces[1], "S")) {
    isStatic = true;
  } else if (!StringUtils::IEquals(pieces[1], "V")) {
    throw runtime_error("Could not parse modification string '" + modString + "'");
  }
  double deltaMass = StringUtils::FromString<double>(pieces[2]);
  ModPosition position = ANY;
  if (pieces.size() == 4) {
    if (pieces[3] == "n") {
      position = PEPTIDE_N;
    } else if (pieces[3] == "c") {
      position = PEPTIDE_C;
    } else if (pieces[3] == "N") {
      position = PROTEIN_N;
    } else if (pieces[3] == "C") {
      position = PROTEIN_C;
    } else {
      throw runtime_error("Could not parse modification string '" + modString + "'");
    }
  }
  const ModificationDefinition* definition =
    ModificationDefinition::Find(deltaMass, isStatic, position);
  string aa;
  if (peptide) {
    char* seq = peptide->getSequence();
    aa.push_back(seq[index]);
    free(seq);
  }
  if (definition == NULL) {
    // Create new modification
    definition = ModificationDefinition::New(
      !aa.empty() ? aa : "X", deltaMass, position, isStatic);
  } else if (!aa.empty() &&
             definition->AminoAcids().find(aa[0]) == definition->AminoAcids().end()) {
    definition = ModificationDefinition::New(
      aa, definition->DeltaMass(), definition->Position(), definition->Static(),
      definition->PreventsCleavage(), definition->PreventsXLink());
  }

  return Modification(definition, index);
}

unsigned char Modification::Index() const {
  return index_;
}

const ModificationDefinition* Modification::Definition() const {
  return mod_;
}

double Modification::DeltaMass() const {
  return mod_->DeltaMass();
}

bool Modification::Static() const {
  return mod_->Static();
}

ModPosition Modification::Position() const {
  return mod_->Position();
}

char Modification::Symbol() const {
  return mod_->Symbol();
}

bool Modification::PreventsCleavage() const {
  return mod_->PreventsCleavage();
}

bool Modification::PreventsXLink() const {
  return mod_->PreventsXLink();
}

// Turns a MODIFIED_AA_T* and its length into a sequence and vector of Modification
// TODO Remove this
void Modification::FromSeq(MODIFIED_AA_T* seq, int length,
                           string* outSeq, vector<Modification>* outMods) {
  *outSeq = string(length, 'X');
  if (outMods != NULL) {
    *outMods = vector<Modification>();
  }
  if (seq != NULL && length != 0) {
    AA_MOD_T** allMods = NULL;
    int modCount = get_all_aa_mod_list(&allMods);
    for (int i = 0; i < length; i++) {
      (*outSeq)[i] = modified_aa_to_char(seq[i]);
      if (outMods == NULL) {
        continue;
      }
      for (int j = 0; j < modCount; j++) {
        if (is_aa_modified(seq[i], allMods[j])) {
          ModPosition position = UNKNOWN;
          switch (aa_mod_get_position(allMods[j])) {
          case N_TERM: position = PEPTIDE_N; break;
          case C_TERM: position = PEPTIDE_C; break;
          }
          string aaList = aa_mod_get_aa_list_string(allMods[j]);
          bool preventsCleavage = aa_mod_get_prevents_cleavage(allMods[j]);
          bool preventsXLink = aa_mod_get_prevents_xlink(allMods[j]);
          double massChange = aa_mod_get_mass_change(allMods[j]);
          const ModificationDefinition* def =
            ModificationDefinition::Find(massChange, false, position);
          if (!def) {
            char symbol = aa_mod_get_symbol(allMods[j]);
            def = ModificationDefinition::NewVarMod(
              aaList, massChange, position, preventsCleavage, preventsXLink, symbol);
          }
          outMods->push_back(Modification(def, (unsigned char)i));
        }
      }
    }
  }
}

// Turns an unmodified sequence and vector of mods into a MODIFIED_AA_T*
// TODO Remove this
MODIFIED_AA_T* Modification::ToSeq(const string& seq, const vector<Modification>& mods) {
  MODIFIED_AA_T* modSeq = (MODIFIED_AA_T*)mycalloc(seq.length() + 1, sizeof(MODIFIED_AA_T));

  for (unsigned i = 0; i < seq.length(); i++) {
    modSeq[i] = char_aa_to_modified(seq[i]);
  }
  modSeq[seq.length()] = MOD_SEQ_NULL;

  for (vector<Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
    if (!i->Static()) {
      modify_aa(&modSeq[i->Index()], get_aa_mod_from_mass((FLOAT_T)i->DeltaMass()));
    }
  }

  return modSeq;
}

ModificationDefinitionContainer::ModificationDefinitionContainer() {
  InitSymbolPool();
}

ModificationDefinitionContainer::~ModificationDefinitionContainer() {
  for (vector<ModificationDefinition*>::const_iterator i = varMods_.begin();
       i != varMods_.end();
       i++) {
    delete *i;
  }
  vector<ModificationDefinition*> staticMods = StaticMods();
  for (vector<ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    delete *i;
  }
}

void ModificationDefinitionContainer::InitSymbolPool() {
  symbolPool_.clear();
  symbolPool_.push_back('*');
  symbolPool_.push_back('#');
  symbolPool_.push_back('@');
  symbolPool_.push_back('^');
  symbolPool_.push_back('~');
  symbolPool_.push_back('%');
  symbolPool_.push_back('$');
  symbolPool_.push_back('&');
  symbolPool_.push_back('!');
  symbolPool_.push_back('?');
  symbolPool_.push_back('+');
}

vector<ModificationDefinition*> ModificationDefinitionContainer::StaticMods(char c) {
  vector<ModificationDefinition*> mods;
  for (map< char, vector<ModificationDefinition*> >::const_iterator i = staticMods_.begin();
       i != staticMods_.end();
       i++) {
    if (c != '\0' && c != i->first) {
      continue;
    }
    for (vector<ModificationDefinition*>::const_iterator j = i->second.begin();
         j != i->second.end();
         j++) {
      if (find(mods.begin(), mods.end(), *j) == mods.end()) {
        mods.push_back(*j);
      }
    }
  }
  return mods;
}

void ModificationDefinitionContainer::Add(ModificationDefinition* def) {
  if (def->Static()) {
    for (set<char>::const_iterator i = def->AminoAcids().begin();
         i != def->AminoAcids().end();
         i++) {
      map< char, vector<ModificationDefinition*> >::const_iterator j = staticMods_.find(*i);
      if (j == staticMods_.end()) {
        staticMods_[*i] = vector<ModificationDefinition*>(1, def);
      } else {
        staticMods_[*i].push_back(def);
      }
    }
  } else {
    varMods_.push_back(def);
  }
  carp(CARP_DEBUG, "Added new modification [%x]: aa[%s] dM[%f] static[%d] symbol[%c]",
       def, StringUtils::Join(def->AminoAcids()).c_str(), def->DeltaMass(),
       def->Static(), def->Symbol());
}

char ModificationDefinitionContainer::NextSymbol() {
  if (!symbolPool_.empty()) {
    char symbol = symbolPool_.front();
    symbolPool_.pop_front();
    return symbol;
  }
  carp(CARP_WARNING, "No more symbols for variable modifications available");
  return '+';
}

void ModificationDefinitionContainer::ConsumeSymbol(char c) {
  for (deque<char>::iterator i = symbolPool_.begin(); i != symbolPool_.end(); i++) {
    if (*i == c) {
      symbolPool_.erase(i);
      return;
    }
  }
}

