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

ModificationDefinition::ModificationDefinition()
  : deltaMass_(numeric_limits<double>::quiet_NaN()), position_(UNKNOWN),
    symbol_('\0'), preventsCleavage_(false), preventsXLink_(false), monoLink_(false) {
}

ModificationDefinition::ModificationDefinition(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool preventsCleavage, bool preventsXLink, bool monoLink, char symbol)
  : deltaMass_(deltaMass), position_(position),
    symbol_(symbol), preventsCleavage_(preventsCleavage), preventsXLink_(preventsXLink),
    monoLink_(monoLink) {
  AddAminoAcids(aminoAcids);
}

ModificationDefinition::~ModificationDefinition() {
}

const ModificationDefinition* ModificationDefinition::New(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool isStatic, bool preventsCleavage, bool preventsXLink, bool monoLink) {
  return isStatic
    ? NewStaticMod(aminoAcids, deltaMass, position, preventsCleavage, preventsXLink)
    : NewVarMod(aminoAcids, deltaMass, position, preventsCleavage, preventsXLink, monoLink, '\0');
}

const ModificationDefinition* ModificationDefinition::NewStaticMod(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool preventsCleavage, bool preventsXLink) {
  // Look for existing
  set<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  for (set<ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    if (MathUtil::AlmostEqual((*i)->deltaMass_, deltaMass, Params::GetInt("mod-precision")) &&
        (position == UNKNOWN || (*i)->position_ == UNKNOWN || (*i)->position_ == position) &&
        (*i)->preventsCleavage_ == preventsCleavage &&
        (*i)->preventsXLink_ == preventsXLink) {
      string added = (*i)->AddAminoAcids(aminoAcids);
      for (string::const_iterator j = added.begin(); j != added.end(); j++) {
        modContainer_.staticMods_[*j].insert(*i);
      }
      if (position != UNKNOWN && (*i)->position_ == UNKNOWN) {
        (*i)->position_ = position;
      }
      return *i;
    }
  }
  // Does not exist yet
  ModificationDefinition* mod = new ModificationDefinition(
    aminoAcids, deltaMass, position, '\0', preventsCleavage, preventsXLink, false);
  modContainer_.Add(mod);
  return mod;
}

const ModificationDefinition* ModificationDefinition::NewVarMod(
  const string& aminoAcids, double deltaMass, ModPosition position,
  bool preventsCleavage, bool preventsXLink, bool monoLink, char symbol) {
  // Look for existing
  for (set<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    if (MathUtil::AlmostEqual((*i)->deltaMass_, deltaMass, Params::GetInt("mod-precision")) &&
        (position == UNKNOWN || (*i)->position_ == UNKNOWN || (*i)->position_ == position) &&
        (*i)->preventsCleavage_ == preventsCleavage &&
        (*i)->preventsXLink_ == preventsXLink &&
        (*i)->monoLink_ == monoLink &&
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
    aminoAcids, deltaMass, position, preventsCleavage, preventsXLink, monoLink, symbol);
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
  set<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  for (set<ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    carp(CARP_INFO, "%s", (*i)->String().c_str());
  }
}

void ModificationDefinition::ListVarMods() {
  // For debugging
  carp(CARP_INFO, "Listing variable modifications");
  for (set<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    carp(CARP_INFO, "%s", (*i)->String().c_str());
  }
}

void ModificationDefinition::Remove(const ModificationDefinition* mod) {
  if (mod == NULL) {
    return;
  } else if (mod->symbol_ == '\0') {
    // Static
    for (set<char>::const_iterator i = mod->aminoAcids_.begin(); i != mod->aminoAcids_.end(); i++) {
      modContainer_.staticMods_[*i].erase((ModificationDefinition*)mod);
    }
  } else {
    // Variable
    modContainer_.varMods_.erase((ModificationDefinition*)mod);
    modContainer_.symbolPool_.push_front(mod->symbol_);
  }
  delete mod;
}

void ModificationDefinition::ClearAll() {
  modContainer_ = ModificationDefinitionContainer();
}

void ModificationDefinition::ClearStaticMods() {
  set<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  for (set<ModificationDefinition*>::const_iterator i = staticMods.begin();
       i != staticMods.end();
       i++) {
    delete *i;
  }
  modContainer_.staticMods_.clear();
}

void ModificationDefinition::ClearVarMods() {
  for (set<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
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
  std::copy(varMods.begin(), varMods.end(), std::back_inserter(mods));
  return mods;
}

vector<const ModificationDefinition*> ModificationDefinition::StaticMods(char c) {
  set<ModificationDefinition*> mods = modContainer_.StaticMods(c);
  vector<const ModificationDefinition*> modsConst;
  for (set<ModificationDefinition*>::const_iterator i = mods.begin();
       i != mods.end();
       i++) {
    modsConst.push_back(*i);
  }
  std::sort(modsConst.begin(), modsConst.end(), ModificationDefinition::SortFunction);
  return modsConst;
}

vector<const ModificationDefinition*> ModificationDefinition::VarMods() {
  vector<const ModificationDefinition*> mods;
  for (set<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    mods.push_back(*i);
  }
  std::sort(mods.begin(), mods.end(), ModificationDefinition::SortFunction);
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

string ModificationDefinition::Title() const {
  return "";
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

bool ModificationDefinition::MonoLink() const {
  return monoLink_;
}

const ModificationDefinition* ModificationDefinition::Find(char symbol) {
  for (set<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
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
  set<ModificationDefinition*> staticMods = modContainer_.StaticMods();
  set<ModificationDefinition*>* mods = isStatic ? &staticMods : &modContainer_.varMods_;
  for (set<ModificationDefinition*>::const_iterator i = mods->begin();
       i != mods->end();
       i++) {
    if ((position == UNKNOWN || position == (*i)->Position()) &&
        MathUtil::AlmostEqual((*i)->deltaMass_, deltaMass, Params::GetInt("mod-precision"))) {
      return *i;
    }
  }
  return NULL;
}

bool ModificationDefinition::SortFunction(const ModificationDefinition* x, const ModificationDefinition* y) {
  if (x == y) {
    return false;
  } else if (x->Static() != y->Static()) {
    return x->Static();
  } else if (x->position_ != y->position_) {
    map<ModPosition, int> posValues;
    posValues[PROTEIN_N] = 0;
    posValues[PEPTIDE_N] = 1;
    posValues[ANY] = 2;
    posValues[PEPTIDE_C] = 3;
    posValues[PROTEIN_C] = 4;
    posValues[UNKNOWN] = 5;
    return posValues[x->position_] - posValues[y->position_] < 0;
  } else if (x->deltaMass_ != y->deltaMass_) {
    return x->deltaMass_ > y->deltaMass_;
  }
  const set<char>& xAminos = x->aminoAcids_;
  const set<char>& yAminos = y->aminoAcids_;
  if (xAminos.size() != yAminos.size()) {
    return xAminos.size() > yAminos.size();
  }
  set<char>::const_iterator i = xAminos.begin();
  set<char>::const_iterator j = yAminos.begin();
  for ( ; i != xAminos.end(); i++, j++) {
    if (*i != *j) {
      return *i < *j;
    }
  }
  if (x->preventsCleavage_ != y->preventsCleavage_) {
    return x->preventsCleavage_;
  } else if (x->preventsXLink_ != y->preventsXLink_) {
    return x->preventsXLink_;
  } else if (x->monoLink_ != y->monoLink_) {
    return x->monoLink_;
  }
  return false;
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

const UnimodDefinition* UnimodDefinition::Get(int unimodId) {
  Unimod::Modification unimod = Unimod::Get(unimodId);
  if (unimod.getTitle().empty()) {
    throw runtime_error("Unknown Unimod ID " + StringUtils::ToString(unimodId));
  }
  // Look for existing
  for (set<ModificationDefinition*>::const_iterator i = modContainer_.varMods_.begin();
       i != modContainer_.varMods_.end();
       i++) {
    if ((*i)->Title() == unimod.getTitle()) {
      return (const UnimodDefinition*)*i;
    }
  }
  // Does not exist yet
  UnimodDefinition* mod = new UnimodDefinition();
  // Initialize modification
  set< pair<char, Unimod::Position> > specificities = unimod.getSpecificities();
  set<Unimod::Position> allPositions;
  for (set< pair<char, Unimod::Position> >::const_iterator i = specificities.begin();
       i != specificities.end();
       i++) {
    if (i->first != 'n' && i->first != 'c') {
      mod->aminoAcids_.insert(i->first);
    }
    allPositions.insert(i->second);
  }
  if (mod->aminoAcids_.empty()) {
    mod->aminoAcids_.insert('X');
  }
  switch (allPositions.size()) {
    case 0: mod->position_ = UNKNOWN; break;
    case 1:
      switch (*allPositions.begin()) {
        case Unimod::ANYWHERE: mod->position_ = ANY; break;
        case Unimod::ANY_N: mod->position_ = PEPTIDE_N; break;
        case Unimod::ANY_C: mod->position_ = PEPTIDE_C; break;
        case Unimod::PROTEIN_N: mod->position_ = PROTEIN_N; break;
        case Unimod::PROTEIN_C: mod->position_ = PROTEIN_C; break;
      }
    default: mod->position_ = ANY; break;
  }
  mod->deltaMass_ = unimod.getMonoMass();
  mod->symbol_ = modContainer_.NextSymbol();
  mod->title_ = unimod.getTitle();
  modContainer_.Add(mod);
  return mod;
}

string UnimodDefinition::Title() const {
  return title_;
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

vector<Modification> Modification::Parse(const string& modString, const string* unmodifiedSequence) {
  vector<Modification> mods;
  vector<string> all = StringUtils::Split(modString, ',');
  for (vector<string>::const_iterator i = all.begin(); i != all.end(); i++) {
    try {
      mods.push_back(ParseOne(StringUtils::Trim(*i), unmodifiedSequence));
    } catch (runtime_error& e) {
      carp(CARP_ERROR, "Error parsing modification string: %s", e.what());
    }
  }
  return mods;
}

Modification Modification::ParseOne(const string& modString, const string* unmodifiedSequence) {
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
  double mass = StringUtils::FromString<double>(pieces[2]);
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
  const ModificationDefinition* mod = ModificationDefinition::Find(mass, isStatic, position);
  string aa;
  if (unmodifiedSequence != NULL) {
    if (unmodifiedSequence->length() > index) {
      aa = string(1, unmodifiedSequence->at(index));
    } else {
      throw runtime_error("Index " + StringUtils::ToString(index) +
                          " is out of bounds for sequence '" + *unmodifiedSequence + "'");
    }
  }
  if (mod == NULL) {
    // Create new modification
    mod = ModificationDefinition::New(!aa.empty() ? aa : "X", mass, position, isStatic);
  } else if (!aa.empty() && mod->AminoAcids().find(aa[0]) == mod->AminoAcids().end()) {
    mod = ModificationDefinition::New(
      aa, mod->DeltaMass(), mod->Position(), mod->Static(),
      mod->PreventsCleavage(), mod->PreventsXLink());
  }
  return Modification(mod, index);
}

unsigned char Modification::Index() const {
  return index_;
}

const ModificationDefinition* Modification::Definition() const {
  return mod_;
}

std::string Modification::Title() const {
  return mod_->Title();
}

const set<char>& Modification::AminoAcids() const {
  return mod_->AminoAcids();
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

bool Modification::MonoLink() const {
  return mod_->MonoLink();
}

void Modification::FromSeq(const string& seq,
                           string* outSeq, vector<Modification>* outMods) {
  size_t i = 0, end = seq.length();
  if (seq.length() >= 5 && seq[1] == '.' && seq[seq.length() - 2] == '.') {
    i += 2;
    end -= 2;
  }
  if (outSeq != NULL) {
    *outSeq = "";
    outSeq->reserve(end - i);
  }
  if (outMods != NULL) {
    *outMods = vector<Modification>();
  }
  int aaCount = 0;
  for (; i < end; i++) {
    char c = seq[i];
    if ('A' <= c && c <= 'Z') {
      ++aaCount;
      if (outSeq != NULL) {
        *outSeq += c;
      }
      continue;
    }
    const ModificationDefinition* mod;
    string aa = i > 0 ? string(1, seq[i - 1]) : "";
    bool handledMod = false;
    switch (c) {
    case '[': {
      // Bracket mod
      size_t j = i;
      do {
        if (++j >= seq.length()) {
          throw runtime_error("Unclosed '[' in sequence '" + seq + "'");
        }
      } while (seq[j] != ']');
      string contents = seq.substr(i + 1, j - (i + 1));
      double mass;
      if (StringUtils::TryFromString(contents, &mass)) {
        if ((mod = ModificationDefinition::Find(mass, false)) == NULL) {
          // Maybe it is static?
          if ((mod = ModificationDefinition::Find(mass, true)) == NULL) {
            mod = ModificationDefinition::NewVarMod(aa, mass, ANY); // TODO N/C mods?
          }
        }
      } else if (StringUtils::StartsWith(contents, "UNIMOD:")) {
          string unimodIdStr = contents.substr(7);
          int unimodId;
          if (!StringUtils::TryFromString(unimodIdStr, &unimodId)) {
            throw runtime_error("Invalid Unimod ID '" + unimodIdStr + "' in sequence '" + seq + "'");
          }
          mod = UnimodDefinition::Get(unimodId);
      } else {
        throw runtime_error("Could not convert '" + contents + "' to double in sequence '" + seq + "'");
      }
      i = j;
      handledMod = true;
      break;
    }
    case 'n':
      if (aaCount == 0) {
        continue;
      }
      break;
    case 'c':
      bool aaRemaining = false;
      for (size_t j = i + 1; j < end; j++) {
        if ('A' <= seq[j] && seq[j] <= 'Z') {
          aaRemaining = true;
          break;
        }
      }
      if (!aaRemaining) {
        continue;
      }
      break;
    }
    if (!handledMod) {
      // Symbol mod
      if ((mod = ModificationDefinition::Find(c)) == NULL) {
        throw runtime_error("Invalid character '" + string(1, c) + "' in sequence '" + seq + "'");
      }
    }
    // Add AA
    if (!aa.empty() && mod->AminoAcids().find(aa[0]) == mod->AminoAcids().end()) {
      mod = ModificationDefinition::New(
        aa, mod->DeltaMass(), mod->Position(), mod->Static(),
        mod->PreventsCleavage(), mod->PreventsXLink());
    }
    if (outMods != NULL && !mod->Static()) {
      outMods->push_back(Modification(mod, aaCount > 0 ? aaCount - 1 : 0));
    }
  }
  if (outSeq != NULL) {
    outSeq->reserve(outSeq->length());
  }
}

// Turns a MODIFIED_AA_T* and its length into a sequence and vector of Modification
// TODO Remove this
void Modification::FromSeq(MODIFIED_AA_T* seq, int length,
                           string* outSeq, vector<Modification>* outMods) {
  if (outSeq != NULL) {
    *outSeq = string(length, 'X');
  }
  if (outMods != NULL) {
    *outMods = vector<Modification>();
  }
  if (seq != NULL && length != 0) {
    AA_MOD_T** allMods = NULL;
    int modCount = get_all_aa_mod_list(&allMods);
    for (int i = 0; i < length; i++) {
      if (outSeq != NULL) {
        (*outSeq)[i] = modified_aa_to_char(seq[i]);
      }
      if (outMods == NULL) {
        continue;
      }
      for (int j = 0; j < modCount; j++) {
        if (allMods[j]->isModified(seq[i])) {
          ModPosition position = UNKNOWN;
          switch (allMods[j]->getPosition()) {
          case N_TERM: position = PEPTIDE_N; break;
          case C_TERM: position = PEPTIDE_C; break;
          }
          string aaList = allMods[j]->getAAListString();
          bool preventsCleavage = allMods[j]->getPreventsCleavage();
          bool preventsXLink = allMods[j]->getPreventsXLink();
          bool monoLink = allMods[j]->getMonoLink();
          double massChange = allMods[j]->getMassChange();
          const ModificationDefinition* def =
            ModificationDefinition::Find(massChange, false, position);
          if (!def) {
            char symbol = allMods[j]->getSymbol();
            def = ModificationDefinition::NewVarMod(
              aaList, massChange, position, preventsCleavage, preventsXLink, monoLink, symbol);
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
  MODIFIED_AA_T* modSeq = newModSeq();
  for (unsigned i = 0; i < seq.length(); i++) {
    modSeq[i] = char_aa_to_modified(seq[i]);
  }
  modSeq[seq.length()] = MOD_SEQ_NULL;

  for (vector<Modification>::const_iterator i = mods.begin(); i != mods.end(); i++) {
    if (!i->Static()) {
      get_aa_mod_from_mass((FLOAT_T)i->DeltaMass())->modify(&modSeq[i->Index()]);
    }
  }

  return modSeq;
}

bool Modification::SortFunction(const Modification& x, const Modification& y) {
  if (x.mod_ == y.mod_) {
    return x.index_ < y.index_;
  } else if (x.Static() != y.Static()) {
    return x.Static();
  } else if (x.index_ != y.index_) {
    return x.index_ < y.index_;
  }
  return ModificationDefinition::SortFunction(x.mod_, y.mod_);
}

ModificationDefinitionContainer::ModificationDefinitionContainer() {
  InitSymbolPool();
}

ModificationDefinitionContainer::~ModificationDefinitionContainer() {
  for (set<ModificationDefinition*>::const_iterator i = varMods_.begin();
       i != varMods_.end();
       i++) {
    delete *i;
  }
  set<ModificationDefinition*> staticMods = StaticMods();
  for (set<ModificationDefinition*>::const_iterator i = staticMods.begin();
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

set<ModificationDefinition*> ModificationDefinitionContainer::StaticMods(char c) {
  set<ModificationDefinition*> mods;
  for (map< char, set<ModificationDefinition*> >::const_iterator i = staticMods_.begin();
       i != staticMods_.end();
       i++) {
    if (c != '\0' && c != i->first) {
      continue;
    }
    for (set<ModificationDefinition*>::const_iterator j = i->second.begin();
         j != i->second.end();
         j++) {
      if (find(mods.begin(), mods.end(), *j) == mods.end()) {
        mods.insert(*j);
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
      map< char, set<ModificationDefinition*> >::const_iterator j = staticMods_.find(*i);
      if (j == staticMods_.end()) {
        staticMods_[*i] = set<ModificationDefinition*>();
      }
      staticMods_[*i].insert(def);
    }
  } else {
    varMods_.insert(def);
  }
  carp(CARP_DEBUG, "Added new modification [%x]: aa[%s] dM[%f] static[%d] symbol[%c]",
       def, StringUtils::Join(def->AminoAcids()).c_str(), def->DeltaMass(),
       def->Static(), !def->Static() ? def->Symbol() : ' ');
}

char ModificationDefinitionContainer::NextSymbol() {
  if (!symbolPool_.empty()) {
    char symbol = symbolPool_.front();
    symbolPool_.pop_front();
    return symbol;
  }
  //carp_once(CARP_WARNING, "No more symbols for variable modifications available");
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

