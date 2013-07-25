// Benjamin Diament
//

#include <stdio.h>
#include <string>
#include <deque>
#include <vector>
#include <algorithm>
#include <gflags/gflags.h>
#include "abspath.h"
#include "records.h"
#include "records_to_vector-inl.h"
#include "header.pb.h"
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "mass_constants.h"
#include "modifications.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(proteins, "", "File of proteins corresponding to peptides, as raw_proteins.proto");
DEFINE_string(in_peptides, "", "File of unfragmented peptides, as peptides.proto");
DEFINE_string(out_peptides, "", "File to create of fragmented peptides, as peptides.proto");

typedef deque<pb::Peptide*> PepQ;
typedef PepQ::iterator PepIter;

class PeptideStream {
 public:
  PeptideStream(PepQ* pep_q, double added_mass)
    : pep_q_(pep_q),
      is_heaviest_(false),
      added_mass_(added_mass),
      iter_(pep_q->begin()) {
  }

  virtual ~PeptideStream() {}

  void MarkAsHeaviest() { is_heaviest_ = true; }

  bool Done() { return iter_ == pep_q_->end(); }

  pb::Peptide* Current() const { return *iter_; }

  virtual void Advance() {
    if (is_heaviest_) {
      assert(iter_ == pep_q_->begin());
      printf("Deleting mass %f\n", (*iter_)->mass());
      delete *iter_;
      pep_q_->pop_front();
    }
    ++iter_;    
  }

  double Mass() const { return (*iter_)->mass() + added_mass_; }

 protected:
  PepIter iter_;
  PepQ* pep_q_;

 private:
  bool is_heaviest_;
  double added_mass_;
};

class LightestPeptideStream : public PeptideStream {
 public:
  LightestPeptideStream(PepQ* pep_q, HeadedRecordReader* reader)
    : PeptideStream(pep_q, 0), reader_(reader) {
    Read();
  }

  virtual void Advance() {
    PeptideStream::Advance();
    assert(iter_ == pep_q_->end());
    Read();
  }

  virtual ~LightestPeptideStream() {}

 private:
  void Read() {
    if (reader_->Done())
      return;
    pb::Peptide* p = new pb::Peptide;
    reader_->Read(p);
    pep_q_->push_back(p);
    assert(*iter_ == p);
  }

  HeadedRecordReader* reader_;
};

struct greater_ps
  : public binary_function<PeptideStream*, PeptideStream*, bool> {
  bool operator()(PeptideStream* x, PeptideStream* y) {
    return x->Mass() > y->Mass();
  }
};

class StreamHeap {
 public:
  StreamHeap(HeadedRecordReader* reader) {
    Init(new LightestPeptideStream(&pep_q_, reader));
  }

  void AddStream(double added_mass) {
    CHECK(added_mass > 0); // also, no two added masses should be the same!
    Init(new PeptideStream(&pep_q_, added_mass));
  }
  
  void DoneAddingStreams() {
    heap_end_ = heap_.end();
    if (!Done()) {
      HIter heaviest = min_element(heap_.begin(), heap_end_, greater_ps());
      (*heaviest)->MarkAsHeaviest();
      make_heap(heap_.begin(), heap_end_, greater_ps());
      pop_heap(heap_.begin(), heap_end_, greater_ps());      
    }
  }

  bool Done() const { return heap_.begin() == heap_end_; }

  PeptideStream* Current() const { return *(heap_end_ - 1); }

  void Advance() {
    PeptideStream* ps = Current();
    ps->Advance();
    if (ps->Done()) {
      delete ps;
      --heap_end_;
    } else {
      push_heap(heap_.begin(), heap_end_, greater_ps());
    }
    if (!Done())
      pop_heap(heap_.begin(), heap_end_, greater_ps());
  }

 private:
  void Init(PeptideStream* str) {
    if (!str->Done()) {
      heap_.push_back(str);
    } else {
      delete str;
    }
  }

  PepQ pep_q_;
  vector<PeptideStream*> heap_;
  typedef vector<PeptideStream*>::iterator HIter;
  HIter heap_end_;
};

void Test(HeadedRecordReader* reader, HeadedRecordWriter* writer) {
  StreamHeap str(reader);
  str.AddStream(10);
  str.AddStream(20);

  for (str.DoneAddingStreams(); !str.Done(); str.Advance()) {
    printf("Mass: %f\n", str.Current()->Mass());
    writer->Write(str.Current()->Current());
  }
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  vector<const pb::Protein*> proteins;
  CHECK(ReadRecordsToVector<pb::Protein>(&proteins, FLAGS_proteins));

  pb::Header orig_header, new_header;
  HeadedRecordReader reader(FLAGS_in_peptides, &orig_header);
  CHECK(orig_header.file_type() == pb::Header::PEPTIDES);
  CHECK(orig_header.has_peptides_header());

  pb::ModTable* mod_table = orig_header.mutable_peptides_header()->mutable_mods();
  CHECK(MassConstants::Init(mod_table));
  VariableModTable var_mod_table;
  var_mod_table.Init(*mod_table);
  var_mod_table.SerializeUniqueDeltas(mod_table);

  new_header.set_file_type(pb::Header::PEPTIDES);
  pb::Header_PeptidesHeader* subheader = new_header.mutable_peptides_header();
  subheader->CopyFrom(orig_header.peptides_header());
  // subheader->set_has_mods(true);
  pb::Header_Source* source = new_header.add_source();
  source->mutable_header()->CopyFrom(orig_header);
  source->set_filename(AbsPath(FLAGS_in_peptides));
  HeadedRecordWriter writer(FLAGS_out_peptides, new_header);
  CHECK(reader.OK());
  CHECK(writer.OK());

  Test(&reader, &writer);

  return 0;
}
