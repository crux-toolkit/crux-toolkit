// Benjamin Diament
//
// Generate peptides (enzymatic or non-enzymatic) from a set of proteins.
// 
// Example command-line:
// make_peptides [ options ] --proteins=<raw_proteins.proto input file> \
//                           --peptides=<peptides.proto output file> 
// See other command-line options below.
//
// Implementation overview
//
// We do not wish to write an unsorted list of all the Peptides and then sort
// them separately, as we expect this would be slow. Instead, we use a heap to
// write the peptides in sorted order as we generate them. We do this as
// follows.
//
// All proteins sequences are stored in memory (This isn't scalable, but the
// same work could be done in chunks that fit in memory and then merged
// later. Such merging is not implemented. TODO 246)
//
// To each amino acid in all the protein sequences we assign a Peptide (see
// class definition below) beginning at that position. This Peptide's length
// will be incremented over the course of the run until either the maximum
// peptide mass or length is exceeded. Barring exceeding the mass and length
// limits, therefore, there are always as many Peptides in memory as there are
// amino acids.
//
// We always maintain a heap of all these Peptides, arranged by mass.
// We think of identical Peptides (the same amino acid sequence) as a group.
// We ensure that a group of identical Peptides will come off the heap 
// consecutively, ordered by protein_id and position within the protein.
// To guarantee this we must use a multi-critereon comparator for the heap.
// (See Peptide::Compare()).
//
// Until there are no Peptides remaining, we repeatedly remove the lightest
// group of Peptides from the heap, create a protocol buffer representing that
// peptide and all its locations, and store it to disk as a record. Then we 
// grow each Peptide in the group by one amino acid and insert it back into the
// heap (assuming it hasn't exceeded the mass and length limits).

#include <stdio.h>
#include <iostream>
#include <string>
#include <climits>
#include <vector>
#include <algorithm>
#include <bitset>
#include <gflags/gflags.h>
#include "records.h"
#include "records_to_vector-inl.h"
#include "header.pb.h"
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "mass_constants.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)
typedef pb::Header_PeptidesHeader Settings;

///////////////////////////////////////////////////////////////////////////////
// Represents a peptide begining at a fixed position. One instance appears on
// the heap for each position. The length varies over the course of the run.
// (See comments above.)
//
// After construction, call Init(), which will return true if the peptide
// should enter the heap. (It might not if it can't be extended to meet
// the minimum mass and length threshold.)
//
// Subclasses may override Advance().
class Peptide {
 public:
  Peptide(const pb::Protein& protein, int pos, int max_len = INT_MAX)
    : protein_id_(protein.id()),
      pos_(pos),
      residues_(protein.residues().data() + pos),
      protein_end_(min(max_len, int(protein.residues().length() - pos))),
      length_(0),
      mass_(monoisotopic_precursor_ ? MassConstants::fixp_mono_h2o
	    : MassConstants::fixp_avg_h2o) {
    assert(protein_end_ >= 0); 
  }

  int Compare(const Peptide& other, bool allow_equal) const {
    // To ensure that identical Peptides come off the heap consecutively and
    // ordered by protein and position.
    // Return -1 if this < other, return 1 if this > other.
    // If allow_equal, return 0 when Peptides are identical.

    // mass first
    if (mass_ > other.mass_)
      return 1;
    if (mass_ < other.mass_)
      return -1;
    // then length
    if (length_ < other.length_)
      return -1;
    if (length_ > other.length_)
      return 1;
    // if all shortcuts fail, have to check the chars
    const char* u = residues_;
    const char* v = other.residues_;
    for (int i = 0; i < length_; ++i, ++u, ++v) {
      int diff = int(*u) - int(*v);
      if (diff != 0)
        return diff;
    }
    if (allow_equal)
      return 0;

    if (protein_id_ < other.protein_id_)
      return -1;
    if (protein_id_ > other.protein_id_)
      return 1;
    return pos_ - other.pos_;
  }

  virtual bool Advance() {
    // Increment length. Return false if limits exceeded.
    if (protein_end_ == 0 || length_ >= max_length_)
      return false;
    --protein_end_;
    mass_ += AAMass(residues_[length_++]);
    return mass_ <= max_mass_;
  }

  int StartPos() { return pos_; }
  int EndPos() { return pos_ + length_; }

  bool Init() {
    while (length_ < min_length_ || mass_ < min_mass_)
      if (!Advance())
        return false;
    return true;
  }

  void Print() {
    // Debugging
    cerr << mass_ << "\t  len:" << length_ << "\t  prot:" 
      << protein_id_ << "\t  pos:" << pos_ << "\n";
  }

  void StartPB(pb::Peptide* pb_peptide) {
    // Initialize a protocol buffer with the first Peptide in a group.
    pb_peptide->set_id(peptide_count_++);
    pb_peptide->set_mass(MassConstants::ToDouble(mass_));
    pb_peptide->set_length(length_);
    pb_peptide->mutable_first_location()->set_protein_id(protein_id_);
    pb_peptide->mutable_first_location()->set_pos(pos_);
  }

  void BuildAuxLocationsPB(pb::AuxLocation* aux_location) {
    // Add another location to a protocol buffer representing a group of
    // locations for the current peptide.
    pb::Location* location = aux_location->add_location();
    location->set_protein_id(protein_id_);
    location->set_pos(pos_);
  }

  static bool SetMinMaxMassAndLength(const Settings& settings) {
    // Returns true if successful (settings are valid).
    if (!settings.has_min_length() || 
        !settings.has_max_length() ||
        !settings.has_min_mass() ||
        !settings.has_max_mass()) 
      return false;
    min_length_ = settings.min_length();
    max_length_ = settings.max_length();
    min_mass_ = MassConstants::ToFixPt(settings.min_mass());
    max_mass_ = MassConstants::ToFixPt(settings.max_mass());
    monoisotopic_precursor_ = (settings.has_monoisotopic_precursor() &&
			       settings.monoisotopic_precursor());
    return true;
  }

 protected:
  static FixPt AAMass(char aa) {
    return monoisotopic_precursor_ ? MassConstants::fixp_mono_table[aa]
      : MassConstants::fixp_avg_table[aa];
  }

  static int peptide_count_;
  static int  min_length_, max_length_;
  static FixPt min_mass_, max_mass_;
  static bool monoisotopic_precursor_;

  int protein_id_;
  unsigned short int pos_;
  const char* residues_; // points into protein sequence
  unsigned short int protein_end_; // remaining residues in protein
  unsigned short int length_;
  FixPt mass_;
};

int Peptide::peptide_count_ = 0;
int Peptide::min_length_;
int Peptide::max_length_;
FixPt Peptide::min_mass_;
FixPt Peptide::max_mass_;
bool Peptide::monoisotopic_precursor_;

struct greater_peptide : public binary_function<Peptide*, Peptide*, bool> {
  bool operator()(Peptide* x, Peptide* y) {
    if (x == y)
      return false;
    return x->Compare(*y, false) > 0;
  }
};

class ReversePeptide : public Peptide {
 public:
  // pos should be given at the *beginning* of the peptide. max_len must be
  // provided accurately (not e.g. INT_MAX). As for Peptide, pos + max_len
  // points just beyond where the peptide ends.
  ReversePeptide(const pb::Protein& protein, int pos, int max_len)
    : Peptide(protein, pos + max_len, 0) {
    CHECK(protein.residues().length() - pos >= max_len);
    protein_end_ = max_len;
  }

  virtual bool Advance() {
    // Same as parent, but move backwards
    if (protein_end_ == 0 || length_ >= max_length_)
      return false;
    --protein_end_;
    mass_ += AAMass(*--residues_);
    ++length_;
    --pos_;
    return mass_ <= max_mass_;
  }
};

// Make only one peptide (if within limits) of length len (or up to protein end)
class WholePeptide : public Peptide {
 public:
  WholePeptide(const pb::Protein& protein, int pos, int len = INT_MAX)
    : Peptide(protein, pos, len) {
  }

  virtual bool Advance() {
    while (Peptide::Advance())
      if (protein_end_ == 0)
        return true;
    return false;
  }
};

class SelectivePeptide : public Peptide {
 public:
  SelectivePeptide(const pb::Protein& protein, const vector<bool>& mask,
                   int pos, int len = INT_MAX)
    : Peptide(protein, pos, len), mask_(mask) {
  }

  virtual bool Advance() {
    while (true) {
      if (!Peptide::Advance())
        return false;
      if (mask_[EndPos()])
        return true;
    }
  }

 private:
  const vector<bool>& mask_;
};

///////////////////////////////////////////////////////////////////////////////
// Cleaves a protein according to the specification of an enzyme as described in
// the construction parameters. Returns a vector<int> of the start positions of
// the resulting peptides.
class Enzyme {
 public:
  // 'left' and 'right' are sets of amino acids. A cleavage will occur at a
  // given site if, to the left is an amino acid in 'left' and to the right is
  // an amino acid in 'right'. If invert_XXX is true, then the cleavage occurs
  // *unless* there is an amino acid in the corresponding set.
  Enzyme(const char* left, bool invert_left,
         const char* right, bool invert_right) {
    Init(left, invert_left, &left_);
    Init(right, invert_right, &right_);
  }

  bool Left(char c) const { return left_[Index(c)]; }
  bool Right(char c) const { return right_[Index(c)]; }

  // Digest protein to find cleavage points.
  void Cleave(const pb::Protein& protein, vector<int>* result) const {
    // marks FIRST position of peptides resulting from digestion but excludes
    // (implicit) position 0 from list. A final entry is inserted marking the
    // end of the last peptide, just beyond its last amino acid.
    result->clear();
    int len = protein.residues().length();
    for (int i = 0; i < len - 1; ++i)
      if (Left(protein.residues()[i]) && Right(protein.residues()[i+1]))
        result->push_back(i+1);
    result->push_back(len);
  }

  static vector<bool>* AsBitVec(vector<int>* cleavages) {
    // caller will own result (and must delete it).
    int len = cleavages->back() + 1;
    vector<bool>* bitvec = new vector<bool>(len, false);
    bitvec->clear();
    vector<int>::iterator i = cleavages->begin();
    for (; i != cleavages->end(); ++i)
      (*bitvec)[*i] = true;
    return bitvec;
  }

  // Map enzyme names to variables
  static const Enzyme* FromName(const string& name);

 private:
  void Init(const char* chars, bool invert, bitset<26>* bits) {
    bits->reset();
    if (invert)
      bits->set();
    bool val = !invert;
    for (; *chars; ++chars)
      (*bits)[Index(*chars)] = val;
  }

  static int Index(char c) {
    CHECK(c >= 'A' && c <= 'Z');
    return c - 'A';
  }

  bitset<26> left_;
  bitset<26> right_;
  
  // We predefine the following enzymes:
  static Enzyme trypsin;
  static Enzyme chymotrypsin;
  static Enzyme elastase;
  static Enzyme clostripain;
  static Enzyme cyanogen_bromide;
  static Enzyme idosobenzoate;
  static Enzyme proline_endopeptidase;
  static Enzyme staph_protease;
  static Enzyme modified_chymotrypsin;
  static Enzyme elastase_trypsin_chymotrypsin;
  static Enzyme aspn;
};

Enzyme Enzyme::trypsin("KR", false, "P", true);
Enzyme Enzyme::chymotrypsin("FWY", false, "P", true);
Enzyme Enzyme::elastase("ALIV", false, "P", true);
Enzyme Enzyme::clostripain("R", false, "", true);
Enzyme Enzyme::cyanogen_bromide("M", false, "", true);
Enzyme Enzyme::idosobenzoate("W", false, "", true);
Enzyme Enzyme::proline_endopeptidase("P", false, "", true);
Enzyme Enzyme::staph_protease("E", false, "", true);
Enzyme Enzyme::modified_chymotrypsin("FWYL", false, "P", true);
Enzyme Enzyme::elastase_trypsin_chymotrypsin("ALIVKRWFY", false, "P", true);
Enzyme Enzyme::aspn("", true, "D", false);

// Map enzyme names to variables
const Enzyme* Enzyme::FromName(const string& name) {
  if (name == "trypsin")
    return &trypsin;
  if (name == "chymotrypsin") 
    return &chymotrypsin;
  if (name == "elastase") 
    return &elastase;
  if (name == "clostripain") 
    return &clostripain;
  if (name == "cyanogen-bromide") 
    return &cyanogen_bromide;
  if (name == "idosobenzoate") 
    return &idosobenzoate;
  if (name == "proline-endopeptidase") 
    return &proline_endopeptidase;
  if (name == "staph-protease") 
    return &staph_protease;
  if (name == "modified-chymotrypsin") 
    return &modified_chymotrypsin;
  if (name == "elastase-trypsin-chymotrypsin") 
    return &elastase_trypsin_chymotrypsin;
  if (name == "aspn") 
    return &aspn;
  return NULL; // No match
}

///////////////////////////////////////////////////////////////////////////////
// Create and manage the heap of peptides; write out the peptide file.
// Subclasses should override AddProtein() to initialize the heap by calling
// AddPeptide() for each peptide in the initial heap.
class PeptideHeap {
 public:
  virtual ~PeptideHeap() {}

  void MakePeptides(const vector<const pb::Protein*>& proteins,
                    HeadedRecordWriter* peptide_writer, 
                    HeadedRecordWriter* aux_loc_writer);

 protected:
  bool AddPeptide(Peptide* new_peptide) {
    if (new_peptide->Init()) {
      peptides_.push_back(new_peptide);
      return true;
    }
    delete new_peptide;
    return false;
  }

 private:
  typedef vector<Peptide*> PepVec;
  typedef PepVec::iterator PepIter;

  virtual void AddProtein(const pb::Protein& protein) = 0;

  void GetGroup(pb::Peptide* pb_peptide, pb::AuxLocation* aux_location);
  void PrintHeap(const string& name, PepIter iter, PepIter end);

  PepVec peptides_;
};

void PeptideHeap::PrintHeap(const string& name, PepIter iter, PepIter end) {
  // Debugging
  cerr << "=========  " << name << "  =========\n";
  for (; iter != end; ++iter) {
    cerr << "     ";
    (*iter)->Print();
  }
}

void PeptideHeap::GetGroup(pb::Peptide* pb_peptide, pb::AuxLocation* aux_location) {
  pb_peptide->Clear();
  aux_location->Clear();
  
  // PrintHeap("Start", peptides_.begin(), peptides_.end());
  
  // Pop a group off the heap. 
  // The group will sit at the end of the peptides_ vector in reverse order.
  // (STL's pop_heap() operation puts the smallest heap element at the end of 
  // the heap range so we call pop_heap() once for each peptide in the group). 
  Peptide* peptide = peptides_.front(); // smallest heap element
  vector<Peptide*>::iterator group_begin = peptides_.end();
  pop_heap(peptides_.begin(), group_begin--, greater_peptide());
  
  while ((group_begin != peptides_.begin()) &&
         (peptide->Compare(*peptides_.front(), true) == 0))
    pop_heap(peptides_.begin(), group_begin--, greater_peptide());
  
  
  // PrintHeap("Heap Portion", peptides_.begin(), group_begin);
  // PrintHeap("Group", group_begin, peptides_.end());
  
  // Fill in pb_peptide with group
  vector<Peptide*>::iterator iter = peptides_.end() - 1;
  (*iter)->StartPB(pb_peptide);

  // Add the auxiliary locations
  for (--iter; iter != group_begin - 1; --iter)    
    (*iter)->BuildAuxLocationsPB(aux_location);
  
  // Advance all group members, deleting any that are done.
  iter = group_begin;
  vector<Peptide*>::iterator end = peptides_.end();
  while (iter != end) {
    if ((*iter)->Advance()) {
      // re-push still-valid peptides and advance iterator
      push_heap(peptides_.begin(), ++iter, greater_peptide());
    } else { // Peptide has exceeded mass or length bounds.
      delete *iter;
      *iter = *(--end);
      // don't advance iterator
    }
  }
  
  // PrintHeap("After Advance", peptides_.begin(), end);
  
  peptides_.resize(end - peptides_.begin()); // remove deleted peptides
}

void PeptideHeap::MakePeptides(const vector<const pb::Protein*>& proteins,
                               HeadedRecordWriter* peptide_writer, 
                               HeadedRecordWriter* aux_loc_writer) {
  CHECK(peptide_writer->OK());
  CHECK(aux_loc_writer->OK());

  cerr << "Creating Peptide Array" << endl;
  vector<const pb::Protein*>::const_iterator i = proteins.begin();
  for (; i != proteins.end(); ++i) {
    if (!(*i)->residues().empty())
      AddProtein(**i);
  }

  cerr << "heapifying" << endl;
  make_heap(peptides_.begin(), peptides_.end(), greater_peptide());

  pb::Peptide pb_peptide;
  pb::AuxLocation aux_location;
  int aux_loc_idx = 0;
  for (int count = 1; !peptides_.empty(); ++count) {
    GetGroup(&pb_peptide, &aux_location);
    
    // Not all peptides have aux locations associated with them. Check to see
    // if GetGroup added any locations to aux_location. If yes, only then
    // assign the corresponding array index to the peptide and write it out.
    if (aux_location.location_size() > 0) {
      pb_peptide.set_aux_locations_index(aux_loc_idx++);
      aux_loc_writer->Write(&aux_location);
    }

    // Write the peptide AFTER the aux_locations check, in case we added an
    // aux_locations_index to the peptide.
    peptide_writer->Write(&pb_peptide);

    if (count % 100000 == 0)
      cerr << "Wrote " << count << " peptides\n";
  }
}

///////////////////////////////////////////////////////////////////////////////
// Subclasses of PeptideHeap for various peptide restrictions

class NonEnzymaticPeptideHeap : public PeptideHeap {
 private:
  virtual void AddProtein(const pb::Protein& protein) {
    for (int pos = 0; AddPeptide(new Peptide(protein, pos)); ++pos);
  }
};

class NonEnzymaticPeptideHeapReverse : public PeptideHeap {
  // testing only; should give exact same results as NonEnzymaticPeptideHeap
 private:
  virtual void AddProtein(const pb::Protein& protein) {
    for (int len = protein.residues().length();
         AddPeptide(new ReversePeptide(protein, 0, len)); --len);
  }
};

class FullyEnzymaticPeptideHeap : public PeptideHeap {
 public:
  FullyEnzymaticPeptideHeap(const Enzyme* enzyme) : enzyme_(enzyme) {}

  virtual ~FullyEnzymaticPeptideHeap() {}

 private:
  virtual void AddProtein(const pb::Protein& protein) {
    vector<int> cleavages;
    enzyme_->Cleave(protein, &cleavages);
    int begin_pos = 0;
    vector<int>::iterator end_pos = cleavages.begin();
    while (end_pos != cleavages.end()) {
      AddPeptide(new WholePeptide(protein, begin_pos, *end_pos - begin_pos));
      begin_pos = *end_pos++; // get next peptide
    }
  }

  const Enzyme* enzyme_;
};

class PartiallyEnzymaticPeptideHeap : public PeptideHeap {
 public:
  PartiallyEnzymaticPeptideHeap(const Enzyme* enzyme) : enzyme_(enzyme) {}

  virtual ~PartiallyEnzymaticPeptideHeap() {}

 private:
  virtual void AddProtein(const pb::Protein& protein) {
    vector<int> cleavages;
    enzyme_->Cleave(protein, &cleavages);
    int begin_pos = 0;
    vector<int>::iterator end_pos = cleavages.begin();
    while (end_pos != cleavages.end()) {
      int max_len = *end_pos - begin_pos;
      AddPeptide(new ReversePeptide(protein, begin_pos, max_len));
      AddPeptide(new Peptide(protein, begin_pos, max_len - 1));
      begin_pos = *end_pos++; // get next peptide
    }
  }

  const Enzyme* enzyme_;
};

class AllowMissedCleavagesPeptideHeap : public PeptideHeap {
 public:
  AllowMissedCleavagesPeptideHeap(const Enzyme* enzyme, int max_cleavages)
    : enzyme_(enzyme), max_cleavages_(max_cleavages) {
  }

  virtual ~AllowMissedCleavagesPeptideHeap() {
    vector<vector<bool>*>::iterator i = masks_.begin();
    for (; i != masks_.end(); ++i)
      delete *i;
  }

 protected:
  vector<int> cleavages_; // current set of cleavages
  vector<bool>* mask_; // current mask

 private:
  virtual void FlipMaskIfNeeded() {}
  virtual void AddFwdPeptide(const pb::Protein& protein, int pos, int len) = 0;
  virtual void AddRevPeptide(const pb::Protein& protein, int pos, int len) {}

  virtual void AddProtein(const pb::Protein& protein) {
    GetCleavages(protein);
    FlipMaskIfNeeded();
    vector<int>::iterator next_pos = cleavages_.begin();
    vector<int>::iterator last_pos = cleavages_.end() - 1;
    vector<int>::iterator end_pos = next_pos + max_cleavages_;
    if (max_cleavages_ >= cleavages_.size())
      end_pos = last_pos;
    for (int pos = 0; next_pos != cleavages_.end(); pos = *next_pos++) {
      AddFwdPeptide(protein, pos, *end_pos - pos);
      if (end_pos != last_pos)
	end_pos++;
    }
    next_pos = cleavages_.begin();
    end_pos = cleavages_.begin();
    int pos = 0;
    for (; end_pos != cleavages_.end(); ++end_pos) {
      AddRevPeptide(protein, pos, *end_pos - pos);
      if (end_pos - cleavages_.begin() >= max_cleavages_)
	pos = *next_pos++;
    }
  }

  virtual void GetCleavages(const pb::Protein& protein) {
    enzyme_->Cleave(protein, &cleavages_);
    mask_ = Enzyme::AsBitVec(&cleavages_);
    masks_.push_back(mask_);
  }

  vector<vector<bool>*> masks_;
  const Enzyme* enzyme_;
  int max_cleavages_;
};

class FullyEnzymaticAllowMissedCleavagesPeptideHeap
  : public AllowMissedCleavagesPeptideHeap {
 public:
  FullyEnzymaticAllowMissedCleavagesPeptideHeap(const Enzyme* enzyme,
						int max_cleavages)
    : AllowMissedCleavagesPeptideHeap(enzyme, max_cleavages) {
  }

  virtual ~FullyEnzymaticAllowMissedCleavagesPeptideHeap() {}

 private:
  virtual void AddFwdPeptide(const pb::Protein& protein, int pos, int len) {
    AddPeptide(new SelectivePeptide(protein, *mask_, pos, len));
  }
};

class PartiallyEnzymaticAllowMissedCleavagesPeptideHeap
  : public AllowMissedCleavagesPeptideHeap {
 public:
  PartiallyEnzymaticAllowMissedCleavagesPeptideHeap(const Enzyme* enzyme,
						    int max_cleavages)
    : AllowMissedCleavagesPeptideHeap(enzyme, max_cleavages) {
  }

  virtual ~PartiallyEnzymaticAllowMissedCleavagesPeptideHeap() {}

 private:
  virtual void FlipMaskIfNeeded() { mask_->flip(); }

  virtual void AddFwdPeptide(const pb::Protein& protein, int pos, int len) {
    AddPeptide(new SelectivePeptide(protein, *mask_, pos, len));
  }

  virtual void AddRevPeptide(const pb::Protein& protein, int pos, int len) {
    AddPeptide(new ReversePeptide(protein, pos, len));
  }
};

#define CHECK_SETTINGS(x) do { \
  if (!(x)) { \
    fprintf(stderr, "SETTINGS CHECK FAILED: %s\n", #x); \
    return false; \
  } \
} while (false)

bool MakePeptides(pb::Header* header, const string& peptides_file, 
                  const string& aux_locs_file) {

  // header is passed in as the settings desiderata.
  // Returns true if settings specified in header are valid.

  // First, collect the source file of raw proteins, read it and check validity
  CHECK_SETTINGS(header->source_size() == 1);
  pb::Header_Source& source = *header->mutable_source(0);
  CHECK_SETTINGS(source.has_filename());
  CHECK_SETTINGS(!source.has_filetype()); // source.header() is given instead
  string proteins_filename = source.filename();
  vector<const pb::Protein*> proteins;
  pb::Header prot_header;
  cerr << "Reading Proteins" << endl;
  CHECK_SETTINGS(ReadRecordsToVector<pb::Protein>(&proteins, proteins_filename,
                                                  &prot_header));
  CHECK_SETTINGS(prot_header.file_type() == pb::Header::RAW_PROTEINS);
  // TODO: reinstate this. CHECK_SETTINGS(prot_header.has_raw_proteins_header());
  // The raw proteins file is read in. It's a valid source file;
  // remember it as such:
  source.mutable_header()->CopyFrom(prot_header);

  // Now check other desired settings
  CHECK_SETTINGS(header->has_peptides_header());
  const Settings& settings = header->peptides_header();
  CHECK_SETTINGS(Peptide::SetMinMaxMassAndLength(settings));

  CHECK_SETTINGS(settings.has_enzyme() && !settings.enzyme().empty());

//  CHECK_SETTINGS(MassConstants::Init(&settings.mods()));

  // Determine what type of PeptideHeap to create from settings indicated 
  // in header->
  PeptideHeap* heap = NULL;
  if (settings.enzyme() == "none") {
    CHECK_SETTINGS(!settings.has_full_digestion() && 
                   !settings.has_max_missed_cleavages());
    heap = new NonEnzymaticPeptideHeap;
  } else {
    CHECK_SETTINGS(settings.has_full_digestion() &&
                   settings.has_max_missed_cleavages());
    int max_missed_cleavages = settings.max_missed_cleavages();
    if (max_missed_cleavages < 0)
      max_missed_cleavages = INT_MAX;

    const Enzyme* enzyme = Enzyme::FromName(settings.enzyme());
    CHECK_SETTINGS(enzyme != NULL); // Unrecognized enzyme name

    if (settings.full_digestion()) {
      if (max_missed_cleavages != 0) {
        heap = new FullyEnzymaticAllowMissedCleavagesPeptideHeap(
	  enzyme, max_missed_cleavages);
      } else {
	heap = new FullyEnzymaticPeptideHeap(enzyme);
      }
    } else { // Partial digest
      if (max_missed_cleavages != 0) {
        heap = new PartiallyEnzymaticAllowMissedCleavagesPeptideHeap(
          enzyme, max_missed_cleavages);
      } else {
	heap = new PartiallyEnzymaticPeptideHeap(enzyme);
      }
    }
  }
  assert(heap != NULL);

  header->set_file_type(pb::Header::PEPTIDES);
  header->mutable_peptides_header()->set_has_peaks(false);
  HeadedRecordWriter peptide_writer(peptides_file, *header); // put header in outfile

  // Create the auxiliary locations header and writer
  pb::Header aux_locs_header;
  aux_locs_header.set_file_type(pb::Header::AUX_LOCATIONS);
  pb::Header_Source* aux_locs_source = aux_locs_header.add_source();
  aux_locs_source->set_filename(peptides_file);
  aux_locs_source->mutable_header()->CopyFrom(*header);
  HeadedRecordWriter aux_loc_writer(aux_locs_file, aux_locs_header);
  
  heap->MakePeptides(proteins, &peptide_writer, &aux_loc_writer);
  delete heap;

  for (int i = 0; i < proteins.size(); ++i)
    delete proteins[i];

  return true;
}
