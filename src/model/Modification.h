#ifndef CRUX_MODIFICATION_H
#define CRUX_MODIFICATION_H

#include "model/Peptide.h"
#include "util/modifications.h"

#include <deque>
#include <ostream>
#include <set>
#include <string>
#include <vector>

enum ModPosition {
  UNKNOWN, ANY,
  PEPTIDE_N, PEPTIDE_C,
  PROTEIN_N, PROTEIN_C
};

class ModificationDefinition {
public:
  ModificationDefinition(
    const std::string& aminoAcids, double deltaMass, ModPosition position,
    bool preventsCleavage, bool preventsXLink, char symbol);
  virtual ~ModificationDefinition();

  static const ModificationDefinition* New(
    const std::string& aminoAcids, double deltaMass, ModPosition position,
    bool isStatic, bool preventsCleavage = false, bool preventsXLink = false);
  static const ModificationDefinition* NewStaticMod(
    const std::string& aminoAcids, double deltaMass, ModPosition position,
    bool preventsCleavage = false, bool preventsXLink = false);
  static const ModificationDefinition* NewVarMod(
    const std::string& aminoAcids, double deltaMass, ModPosition position,
    bool preventsCleavage = false, bool preventsXLink = false, char symbol = '\0');

  std::string String() const;
  static void ListAll();
  static void ListStaticMods();
  static void ListVarMods();

  static void ClearAll();
  static void ClearStaticMods();
  static void ClearVarMods();

  static std::vector<const ModificationDefinition*> AllMods();
  static std::vector<const ModificationDefinition*> StaticMods(char c = '\0');
  static std::vector<const ModificationDefinition*> VarMods();
  static double DeltaMass(char c, ModPosition position);

  const std::set<char>& AminoAcids() const;
  double DeltaMass() const;
  bool Static() const;
  ModPosition Position() const;
  char Symbol() const;
  bool PreventsCleavage() const;
  bool PreventsXLink() const;

  static const ModificationDefinition* Find(char symbol);
  static const ModificationDefinition* Find(double deltaMass,
    bool isStatic, ModPosition position = UNKNOWN);
protected:
  std::string AddAminoAcids(const std::string& aminoAcids);

  std::set<char> aminoAcids_;
  double deltaMass_;
  ModPosition position_;
  char symbol_;
  bool preventsCleavage_;
  bool preventsXLink_;
};

namespace Crux { class Modification; }
void swap(Crux::Modification& x, Crux::Modification& y);

namespace Crux {
class Modification {
public:
  friend void ::swap(Modification& x, Modification& y);
  friend std::ostream& operator<<(std::ostream& stream, const Modification& mod);

  Modification(const ModificationDefinition* mod, unsigned char index);
  Modification(const Modification& other);
  virtual ~Modification();

  Modification& operator=(Modification other);
  bool operator==(const Modification& other) const;
  bool operator!=(const Modification& other) const;

  std::string String() const;
  static Modification Parse(const std::string& modString, Peptide* peptide);

  unsigned char Index() const;
  const ModificationDefinition* Definition() const;
  double DeltaMass() const;
  bool Static() const;
  ModPosition Position() const;
  char Symbol() const;
  bool PreventsCleavage() const;
  bool PreventsXLink() const;

  static void FromSeq(MODIFIED_AA_T* seq, int length,
                      std::string* outSeq, std::vector<Modification>* outMods);
  static MODIFIED_AA_T* ToSeq(const std::string& seq, const std::vector<Modification>& mods);
protected:
  unsigned char index_; // 0 based position
  const ModificationDefinition* mod_;
};
}

class ModificationDefinitionContainer {
public:
  friend class ModificationDefinition;

  ModificationDefinitionContainer();
  virtual ~ModificationDefinitionContainer();
protected:
  void InitSymbolPool();
  std::vector<ModificationDefinition*> StaticMods(char c = '\0');
  void Add(ModificationDefinition* def);
  char NextSymbol();
  void ConsumeSymbol(char c);

  std::vector<ModificationDefinition*> varMods_;
  std::map< char, std::vector<ModificationDefinition*> > staticMods_;
  std::deque<char> symbolPool_;
};

static ModificationDefinitionContainer modContainer_;

#endif

