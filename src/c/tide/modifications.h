// Benjamin Diament

#include<climits>
#include<algorithm>
#include<iostream>
#include<stdio.h>
#include<tr1/unordered_map>
#include "mod_coder.h"

using namespace std;

#define SHOW(x) { cout << (#x) << ": " << (x) << "\n"; }
#if 0
#define SHOW_ALL(x) {							  \
  cout << (#x) << ":";                                                    \
  for (typeof(x)::const_iterator it = (x).begin(); it != (x).end(); ++it) \
    cout << " " << (*it);                                                 \
  cout << "\n";                                                           \
 }
#endif
#define SHOW_PAIR(x) ((x).first) << ", " << ((x).second)

template<typename T> inline void ShowAll(const T& x) {
  for (typename T::const_iterator it = x.begin(); it != x.end(); ++it)
    cout << " " << (*it);
  cout << "\n";
}

#define SHOW_ALL(x) { cout << (#x) << ":"; ShowAll((x)); }


struct first_eq : public binary_function<pair<int,int>, pair<int,int>, bool> {
  bool operator()(pair<int, int> x, pair<int, int> y) {
    return x.first == y.first;
  }
};

static bool IsAA(char c) {
  const char* AA = "ACDEFGHIKLMNPQRSTVWY";
  for (; (*AA != '\0') && (*AA != c); ++AA);
  return (*AA == c);
}

class VariableModTable {
 public:
  VariableModTable() {}

  bool Init(const pb::ModTable& pb_mod_table) {
    if (pb_mod_table.variable_mod_size() == 0)
      return true;
    vector<double>& UD = unique_delta_;
    for (int i = 0; i < pb_mod_table.variable_mod_size(); ++i) {
      UD.push_back(pb_mod_table.variable_mod(i).delta());
      max_counts_.push_back(pb_mod_table.variable_mod(i).max_count());
    }
    original_deltas_.resize(UD.size());
    copy(UD.begin(), UD.end(), original_deltas_.begin());
    sort(UD.begin(), UD.end());
    UD.resize(unique(UD.begin(), UD.end()) - UD.begin());
    if (UD[0] <= 0) {
      cerr << "ERROR: Non-positive amino acid modification " << UD[0]
	   << " appears in modifications table.\n";
      return false;
    }
    coder_.Init(UD.size());
    tr1::unordered_map<double, int> deltas;
    for (int i = 0; i < UD.size(); ++i)
      deltas[UD[i]] = i;
    for (int i = 0; i < pb_mod_table.variable_mod_size(); ++i) {
      const pb::Modification& mod = pb_mod_table.variable_mod(i);
      const string& aa = mod.amino_acids();
      for (int j = 0; j < aa.size(); ++j)
	possibles_[aa[j]].push_back(make_pair(deltas[mod.delta()], i));
    }
    // Check possibles lists
    for (int i = 0; i < 256; ++i) {
      vector<pair<int, int> >& p = possibles_[i];
      if (!IsAA(char(i)) && p.size() > 0)
	return false;
      sort(p.begin(), p.end());
      vector<pair<int, int> >::iterator dup = adjacent_find(p.begin(), p.end(),
							    first_eq());
      if (dup != p.end()) {
	cerr << "ERROR: Amino acid modification " << char(i) << "+"
	     << unique_delta_[dup->first] << " appears more than once in "
	     << "modifications table.\n";
	return false;
      }
    }
    return true;
  }

  bool Parse(const char* spec_text) {
    pb_mod_table_.Clear();
    int pos = 0;
    while (true) {
      char c;
      int next_pos = -1;
      sscanf(spec_text + pos, " %c%n", &c, &next_pos);
      if (next_pos == -1)
	return Error(spec_text, pos, "Expected modification specification.");
      pos += next_pos - 1;
      unsigned int limit = 0;
      if (c >= '1' && c <= '9') {
	sscanf(spec_text + pos, "%u %n", &limit, &next_pos);
	if (limit == UINT_MAX)
	  return Error(spec_text, pos, "Limit too big.");
	pos += next_pos;
      }
      int aa_len = -1, plus_pos = -1, delta_pos = -1, end_pos = -1;
      sscanf(spec_text + pos, "%*[ACDEFGHIKLMNPQRSTVWY]%n %n+ %n%*[0-9.] %n",
	     &aa_len, &plus_pos, &delta_pos, &end_pos);
      if (aa_len == -1)
	return Error(spec_text, pos, "Expected amino acid symbol.");
      assert(plus_pos != -1);
      if (delta_pos == -1)
	return Error(spec_text, pos + plus_pos, "Expected '+' and "
		     "modification amount.");
      if (end_pos == -1)
	return Error(spec_text, pos + delta_pos, "Expected modification "
		     "amount.");
      if ((limit == 0) && (aa_len != 1))
	return Error(spec_text, pos, "Static modifications must be specified "
		     "for one amino acid at a time.");
      int confirm_end_pos = -1;
      double delta;
      sscanf(spec_text + pos + delta_pos, "%lg %n", &delta, &confirm_end_pos);
      if (delta_pos + confirm_end_pos != end_pos)
	return Error(spec_text, pos + delta_pos, "Cannot parse modification "
		     "amount.");
      if (delta <= 0)
	return Error(spec_text, pos + delta_pos, "Modification amount must be "
		     "positive.");

      pb::Modification* mod;
      if (limit == 0) {
	mod = pb_mod_table_.add_static_mod();
      } else {
	mod = pb_mod_table_.add_variable_mod();
	mod->set_max_count(limit);
      }
      mod->set_amino_acids(string(spec_text + pos, aa_len));
      mod->set_delta(delta);

      pos += end_pos;
      if (spec_text[pos] == '\0')
	break;
      if (spec_text[pos] == ',')
	++pos;
    }
    if (!Init(pb_mod_table_))
      return false;
    return SerializeUniqueDeltas(&pb_mod_table_);
  }

  int NumPoss(char aa) const { return possibles_[aa].size(); }

  int PossMaxCt(char aa, int index) const {
    return possibles_[aa][index].second;
  }

  int PossDeltIx(char aa, int index) const {
    return possibles_[aa][index].first;
  }

  double PossDelta(char aa, int index) const {
    return unique_delta_[PossDeltIx(aa, index)];
  }

#if 0
  const IntPairVec* Possibles(char aa) const {
    return &possibles_[aa];
  }
#endif

  bool SerializeUniqueDeltas(pb::ModTable* pb_mod_table) {
    if (pb_mod_table->unique_deltas_size() == 0) {
      vector<double>::iterator iter = unique_delta_.begin();
      for (; iter != unique_delta_.end(); ++iter)
	pb_mod_table->add_unique_deltas(*iter);
      return true;
    }
    // if unique deltas already specified, just confirm equality
    if (pb_mod_table->unique_deltas_size() != unique_delta_.size())
      return false;
    for (int i = 0; i < unique_delta_.size(); ++i)
      if (unique_delta_[i] != pb_mod_table->unique_deltas(i))
	return false;
    return true;
  }

  const pb::ModTable* ParsedModTable() const { return &pb_mod_table_; }

  const vector<int>* MaxCounts() const { return &max_counts_; }
  const vector<double>* OriginalDeltas() const { return &original_deltas_; }

  int EncodeMod(int aa_index, int unique_delta_index) {
    return coder_.EncodeMod(aa_index, unique_delta_index);
  }

  void Show() {
    SHOW(unique_delta_.size());
    SHOW_ALL(unique_delta_);
    SHOW(max_counts_.size());
    SHOW_ALL(max_counts_);
    const char* aa = "ACDEFGHIKLMNPQRSTVWY";
    for (const char* c = aa; *c; ++c) {
      cout << "possibles_[" << (*c) << "] = ";
      for (IntPairVec::iterator i = possibles_[*c].begin();
	   i != possibles_[*c].end(); ++i) {
	cout << "delta: " << i->first << "(" << unique_delta_[i->first] << ")"
	     << "  max_count: " << i->second << "(" << max_counts_[i->second]
	     << ")  ";
      }
      cout << "\n";
    }
  }

 private:
  bool Error(const char* spec_text, int err_pos, const char* msg) {
    cerr << "Error: couldn't parse modification specification:\n" << spec_text
	 << "\n" << string(err_pos, ' ') << "^\n" << msg << "\n";
    pb_mod_table_.Clear();
    return false;
  }

  typedef vector<pair<int, int> > IntPairVec;
  IntPairVec possibles_[256]; // unique_delta_, max_count_
  vector<double> unique_delta_, original_deltas_;
  vector<int> max_counts_;
  ModCoder coder_;
  pb::ModTable pb_mod_table_;
};

