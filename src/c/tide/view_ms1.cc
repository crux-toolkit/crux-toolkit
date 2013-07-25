#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <gflags/gflags.h>
#include "view_ms1_cmds.pb.h"
#include "spectrum_collection.h"
#include "records.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(request_file, "/tmp/view_ms1_request_pipe", "sequence of requests");
DEFINE_string(response_file, "/tmp/view_ms1_response_pipe", "sequence of responses");

class Handler {
 public:
  Handler() : coll_(NULL) {}

  bool Handle(const pb::ViewMS1Command& cmd) {
    bool ok = true;
    if (ok && cmd.has_load())
      ok = Load(cmd.load());
    if (ok && cmd.has_bucket())
      Bucket(cmd.bucket());
    if (ok && cmd.has_min_mz())
      SetMinMZ(cmd.min_mz());
    if (ok && cmd.has_max_mz())
      SetMaxMZ(cmd.max_mz());
    if (ok && cmd.has_min_rtime())
      SetMinRTime(cmd.min_rtime());
    if (ok && cmd.has_max_rtime())
      SetMaxRTime(cmd.max_rtime());
    if (ok && cmd.has_write_matrix())
      ok = WriteMatrix(cmd.write_matrix());
    return ok;
  }

 private:
  bool Load(const string& filename) {
    delete coll_;
    pb::Header header;
    coll_ = new SpectrumCollection;
    bool ok = coll_->ReadSpectrumRecords(filename, &header);
    if (ok)
      return true;
    delete coll_;
    coll_ = NULL;
    return false;
  }

  void Bucket(double binsize) { binsize_ = binsize; }
  void SetMinMZ(double min_mz) { min_mz_ = min_mz; }
  void SetMaxMZ(double max_mz) { max_mz_ = max_mz; }
  void SetMinRTime(double min_rtime) { min_rtime_ = min_rtime; }
  void SetMaxRTime(double max_rtime) { max_rtime_ = max_rtime; }

  bool WriteMatrix(const string& filename) {
    if (filename[0] != '/') {
      cerr << "filename needs to be absolute\n";
      return false;
    }
    ofstream out(filename.c_str());
    Specs* spectra = coll_->Spectra();
    int num_bins = int((max_mz_ - min_mz_)/binsize_ + 0.99999);
    double bins[num_bins];
    for (int i = 0; i < num_bins; ++i)
      bins[i] = 0;
    for (Specs::iterator i = spectra->begin(); i!=spectra->end(); ++i) {      
      const Spectrum& s = **i;
      for (int j = 0; j < s.Size(); ++j) {
	int bucket = int((s.M_Z(j) - min_mz_)/binsize_);
	bins[bucket] += s.Intensity(j);
      }
      out << s.RTime();
      for (int i = 0; i < num_bins; ++i)
	out << '\t' << bins[i];
      out << '\n';
    }
    return true;
  }

  typedef vector<Spectrum*> Specs;
  SpectrumCollection* coll_;
  double binsize_;
  double min_mz_;
  double max_mz_;
  double min_rtime_;
  double max_rtime_;
};

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_request_file.empty());
  CHECK(!FLAGS_response_file.empty());

  Handler handler;
  RecordReader reader(FLAGS_request_file);
  ofstream out(FLAGS_response_file.c_str());
  while (!reader.Done()) {
    pb::ViewMS1Command cmd;
    CHECK(reader.Read(&cmd));
    bool ok = handler.Handle(cmd);    
    out << (ok ? "OK" : "COMMAND FAILED") << "\n";
    out.flush();
  }
  CHECK(reader.OK());
  out << "Done\n";

  return 0;
}
