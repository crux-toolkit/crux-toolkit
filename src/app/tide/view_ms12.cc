#include <stdio.h>
#include <algorithm>
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

DEFINE_string(ms1_file, "", "original spectra");
DEFINE_string(ms2_file, "", "original ms/ms spectra file");
DEFINE_string(additional_ms2, "", "additional corresponding ms/ms spectra files");
DEFINE_string(matrix_file, "", "output matrix file");
DEFINE_double(mz_bin_size, 5, "size of buckets on m/z axis");
DEFINE_double(min_mz, 0, "lowest m/z");
DEFINE_double(max_mz, 0, "highest m/z");
DEFINE_double(min_rtime, 0, "lowest retention time");
DEFINE_double(max_rtime, 0, "highest retention time");

DEFINE_double(trim, 0.001, "remove this fraction of outliers.");

class MZLine {
 public:
  MZLine() : rtime_(0) {
    assert(num_bins_ > 0);
    arr_ = new double[num_bins_];
    for (int i = 0; i < num_bins_; ++i)
      arr_[i] = 0;
  }

  ~MZLine() { delete[] arr_; }

  static void CheckIndex(int index) {
    assert(index >= 0);
    assert(index < num_bins_);
  }

  double operator[](int index) const { CheckIndex(index); return arr_[index]; }
  double& operator[](int index) { CheckIndex(index); return arr_[index]; }
  double RTime() const { return rtime_; }
  void SetRTime(double rtime) { rtime_ = rtime; }

  static void SetNumBins(int num_bins) { num_bins_ = num_bins; }
  static int NumBins() { return num_bins_; }

  bool operator<(const MZLine& other) const { return rtime_ < other.rtime_; }

 private:
  static int num_bins_;
  double* arr_;
  double rtime_;
};

struct greater_line : public binary_function<const MZLine*, const MZLine*, bool> {
  bool operator()(const MZLine* x, const MZLine* y) { return (*y) < (*x); }
};

struct less_line : public binary_function<const MZLine*, const MZLine*, bool> {
  bool operator()(const MZLine* x, const MZLine* y) { return (*x) < (*y); }
};

int MZLine::num_bins_ = -1;

void Split(vector<string>* v, string s, const string& delimiters) {
  string::size_type last = 0;
  string::size_type pos = s.find_first_of(delimiters, last);

  while (string::npos != pos || string::npos != last) {
    v->push_back(s.substr(last, pos - last));
    last = s.find_first_not_of(delimiters, pos);
    pos = s.find_first_of(delimiters, last);
  }
}

class Handler {
 public:
  Handler() : coll_(NULL) {
    binsize_ = FLAGS_mz_bin_size;
    min_mz_ = FLAGS_min_mz;
    max_mz_ = FLAGS_max_mz;
    min_rtime_ = FLAGS_min_rtime;
    max_rtime_ = FLAGS_max_rtime;
  }

  bool Run() {
    vector<string> additional_ms2;
    Split(&additional_ms2, FLAGS_additional_ms2, ", \t");
    cerr << "Loading ms1 file " << FLAGS_ms1_file << endl;
    if (!Load(FLAGS_ms1_file))
      return false;
    cerr << "Bucketing.\n";
    Bucket();
    cerr << "Discretizing.\n";
    Discretize();
    cerr << "Adding MS2 file " << FLAGS_ms2_file << endl;
    if (!AddMS2(FLAGS_ms2_file, -1))
      return false;
    for (vector<string>::iterator i = additional_ms2.begin();
	 i != additional_ms2.end(); ++i) {
      cerr << "Adding MS2 file " << *i << endl;
      if (!AddMS2(*i, -2))
	return false;
    }
    cerr << "Writing matrix file " << FLAGS_matrix_file << endl;
    if (!WriteMatrix(FLAGS_matrix_file))
      return false;
    cerr << "Done\n";
    return true;
  }

 private:
  typedef vector<Spectrum*> Specs;
  typedef vector<MZLine*> LineVec;

  static SpectrumCollection* LoadColl(const string& filename) {
    pb::Header header;
    SpectrumCollection* coll = new SpectrumCollection;
    bool ok = coll->ReadSpectrumRecords(filename, &header);
    if (!ok) {
      delete coll;
      return NULL;
    }
    return coll;
  }

  bool Load(const string& filename) {
    delete coll_;
    coll_ = LoadColl(filename);
    return (coll_ != NULL);
  }

  int Bin(double mz) { return int((mz - min_mz_) / binsize_); }

  bool IsSorted(const LineVec& lines) {
    return (adjacent_find(lines.begin(), lines.end(), greater_line()) ==
	    lines.end());
  }

  bool AddMS2(const string& filename, int label) {
    SpectrumCollection* coll = LoadColl(filename);
    if (coll == NULL)
      return false;
    Specs* spectra = coll->Spectra();
    LineVec lines;
    lines.reserve(spectra->size());
    for (int i = 0; i < spectra->size(); ++i) {
      const Spectrum& s = *(*spectra)[i];
      double mz = s.PrecursorMZ();
      if (mz >= min_mz_ && mz <= max_mz_) {
	MZLine* line = new MZLine();
	line->SetRTime(s.RTime());
	int bucket = Bin(mz);
	for (int j = bucket - 2; j <= bucket+2; ++j)
	  if (j >= 0 && j < MZLine::NumBins())
	    (*line)[j] = label;
	lines.push_back(line);
      }
    }    
    CHECK(IsSorted(lines));
    LineVec result;
    result.reserve(bins_.size() + lines.size());
    merge(bins_.begin(), bins_.end(), lines.begin(), lines.end(),
	  back_insert_iterator<LineVec>(result), less_line());
    bins_.swap(result);
    return true;
  }

  void Bucket() {
    Specs* spectra = coll_->Spectra();
    MZLine::SetNumBins(int((max_mz_ - min_mz_)/binsize_ + 1));
    bins_.reserve(spectra->size());

    for (int i = 0; i < spectra->size(); ++i) {
      bins_.push_back(new MZLine());
      const Spectrum& s = *(*spectra)[i];
      bins_[i]->SetRTime(s.RTime());
      for (int j = 0; j < s.Size() && s.M_Z(j) < max_mz_; ++j) {
	if (s.M_Z(j) < min_mz_)
	  continue;
	int bucket = Bin(s.M_Z(j));
	(*bins_[i])[bucket] += s.Intensity(j);
      }
    }
    CHECK(IsSorted(bins_));
  }

  void Discretize() {
    const int kColors = 128;
    // Trim
    vector<double> vals;
    for (int i = 0; i < bins_.size(); ++i) {
      for (int j = 0; j < MZLine::NumBins(); ++j) {
	double val = (*bins_[i])[j];
	if (val > 0)
	  vals.push_back(val);
      }
    }
    int nvals = vals.size();
    int trim_point = int(nvals * FLAGS_trim);
    nth_element(vals.begin(), vals.end() - trim_point, vals.end());
    double max_val = vals[vals.size() - trim_point - 1];
    cerr << "MAX VAL: " << max_val << endl;
    double window = max_val/kColors;
    for (int i = 0; i < bins_.size(); ++i) {
      for (int j = 0; j < MZLine::NumBins(); ++j) {
	double& entry = (*bins_[i])[j];
	if (entry >= max_val) {
	  entry = kColors-1;
	} else {
	  entry = int(entry / window);
	  if (entry > kColors-1)
	    entry = kColors-1;
	}
      }
    }
  }

  void Discretize2() {
    const int kColors = 128;

    vector<double> vals;
    for (int i = 0; i < bins_.size(); ++i) {
      for (int j = 0; j < MZLine::NumBins(); ++j) {
	double val = (*bins_[i])[j];
	if (val > 0)
	  vals.push_back(val);
      }
    }
    sort(vals.begin(), vals.end());

    double cutoffs[kColors-1];
    for (int i = 1; i < kColors; ++i)
      cutoffs[i-1] = vals[i*vals.size()/kColors];

    for (int i = 0; i < bins_.size(); ++i) {
      for (int j = 0; j < MZLine::NumBins(); ++j) {
	double& entry = (*bins_[i])[j];
	entry = lower_bound(cutoffs, cutoffs + kColors - 1, entry) - cutoffs;
      }
    }
  }

  bool WriteMatrix(const string& filename) {
    ofstream out(filename.c_str());
    if (!out.is_open()) {
      cerr << "ERROR: failed to open file " << filename << endl;
      return false;
    }
    for (int i = 0; i < bins_.size(); ++i) {
      out << bins_[i]->RTime();
      for (int j = 0; j < MZLine::NumBins(); ++j)
	out << '\t' << (*bins_[i])[j];
      out << '\n';
    }
    out.close();
    return true;
  }

  SpectrumCollection* coll_;
  double binsize_;
  double min_mz_;
  double max_mz_;
  double min_rtime_;
  double max_rtime_;
  LineVec bins_;
};

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  CHECK(!FLAGS_ms1_file.empty());
  CHECK(!FLAGS_matrix_file.empty());

  Handler handler;
  return handler.Run() ? 0 : 1;
}

// for --dmap option:
// for i in ["%d\t%d:0:0\tlabel%d" % (i,i,i) for i in range(256)]: print >>f, i
