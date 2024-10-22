// Benjamin Diament

/*
* In order to update the code for the unimod parser (parse_unimod)
follow the instructions below:
  1. Navigate to crux-toolkit's  bin folder this should be found at the root level of the crux-toolkit folder
  2. Download the latest unimod xml files from http://www.unimod.org/xml/unimod.xml
  3. Run the commad: "python unimod_parser.py unimod.xml > unimod.h"
  4. Copy the unimods.h header to src/model/unimod.h
  5. Modify parse_unimod function as needed
*/

#ifndef TIDE_MODIFICATIONS_H
#define TIDE_MODIFICATIONS_H

#include <stdio.h>

#include <algorithm>
#include <climits>
#include <iostream>
#if defined(_MSC_VER) || defined(DARWIN)
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif
#include <regex>
#include <vector>

#include "header.pb.h"
#include "io/carp.h"
#include "mod_coder.h"
#include "model/Unimod.h"
#include "util/Params.h"
#include <sstream>

using namespace std;

#define SHOW_ONE(x) \
    { cout << (#x) << ": " << (x) << "\n"; }
#define SHOW_PAIR(x) ((x).first) << ", " << ((x).second)

enum mods_spec_type {
    MOD_SPEC,  // table for regular amino acid modifications
    NTPEP,     // type for n-terminal peptide modifications
    CTPEP,     // type for c-terminal peptide modifications
    NTPRO,     // type for n-terminal protein modifications
    CTPRO      // type for c-terminal protein modifications
};

typedef enum mods_spec_type MODS_SPEC_TYPE_T;

template <typename T>
inline void ShowAll(const T& x) {
    for (typename T::const_iterator it = x.begin(); it != x.end(); ++it)
        cout << " " << (*it);
    cout << "\n";
}

#define SHOW_ALL(x)          \
    {                        \
        cout << (#x) << ":"; \
        ShowAll((x));        \
    }

struct first_eq : public binary_function<pair<int, int>, pair<int, int>, bool> {
    bool operator()(pair<int, int> x, pair<int, int> y) {
      return x.first == y.first;
    }
};

static bool IsAA(char c) {
    const char* AA = "ACDEFGHIKJLMNOPQRSTUVWYX";
    for (; (*AA != '\0') && (*AA != c); ++AA)
        ;
    return (*AA == c);
}

class VariableModTable {
   public:
    VariableModTable() { offset_ = 0; }

    /**
     * @brief Loads and validates possibles vector
     *
     * @param pb_mod_table - the type of the vector to load - regular or one of terminals
     * @return true if completed successfully
     * @return false otherwise
     */
    bool Init(const pb::ModTable& pb_mod_table) {
        IntPairVec* possibles;
        if (&pb_mod_table == &pb_mod_table_)
            possibles = possibles_;
        if (&pb_mod_table == &pb_ntpep_mod_table_)
            possibles = possibles_ntpe_;
        if (&pb_mod_table == &pb_ctpep_mod_table_)
            possibles = possibles_ctpe_;
        if (&pb_mod_table == &pb_ntpro_mod_table_)
            possibles = possibles_ntpr_;
        if (&pb_mod_table == &pb_ctpro_mod_table_)
            possibles = possibles_ctpr_;

        if (pb_mod_table.variable_mod_size() == 0)
            return true;
        vector<double>& UD = unique_delta_;

#if defined _MSC_VER || defined DARWIN
        unordered_map<double, int> deltas;
#else
        tr1::unordered_map<double, int> deltas;
#endif
        for (int i = 0; i < UD.size(); ++i)
            deltas[UD[i]] = i;  // lookup table delta value to its position in the _unique_deltas vector

        for (int i = 0; i < pb_mod_table.variable_mod_size(); ++i) {
            const pb::Modification& mod = pb_mod_table.variable_mod(i);
            const string& aa = mod.amino_acids();
            for (int j = 0; j < aa.size(); ++j) {
                // for each AA make a list of (position in _unique_delta, global mod index) pairs
                // global mod index is the number of the given mod in the list of mods specified on the command line.
                possibles[aa[j]].push_back(make_pair(deltas[mod.delta()], i + offset_));
            }
        }
        offset_ += pb_mod_table.variable_mod_size();

        // Check possibles lists
        for (int i = 0; i < 256; ++i) {
            vector<pair<int, int> >& p = possibles[i];
            if (!IsAA(char(i)) && p.size() > 0) {
                return false;
            }
            sort(p.begin(), p.end());
            vector<pair<int, int> >::iterator dup = adjacent_find(p.begin(), p.end(), first_eq());
            if (dup != p.end()) {
                carp(CARP_FATAL,
                     "Amino acid modification %c+%g appears more than once in modifications table.",
                     char(i), unique_delta_[dup->first]);
                return false;
            }
        }
        return true;
    }

    struct ParserResponse {
        std::string actual_data;
        std::vector<std::string> pb_data;
    };

    ParserResponse parse_unimod(const char* spec_text) {
        ParserResponse response;
        std::regex regexPattern("\\[Unimod:[0-9]+\\]");
        std::string mod_specs = spec_text;

        std::string mod_spec, actual_data, token;
        std::vector<std::string> pb_data;
        std::smatch all_match, num_match;
        std::regex numberRegex("\\d+");
        Unimod::Modification mod = Unimod::Modification();
        string isotopic_mass = Params::GetString("isotopic-mass");
        int unimod_id;
        std::vector<std::string> tokens;
        std::stringstream ss(mod_specs);
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        for (std::vector<std::string>::iterator it = tokens.begin(); it != tokens.end(); ++it) {
            mod_spec = *it;

            if (regex_search(mod_spec, all_match, regexPattern)) {
                std::string unimod = all_match.str();
                if (regex_search(unimod, num_match, numberRegex)) {
                    unimod_id = stoi(num_match.str());
                    mod = Unimod::Get(unimod_id);
                    double mass;
                    if (isotopic_mass == "average") {
                        mass = mod.getAvgMass();
                    } else if (isotopic_mass == "mono") {
                        mass = mod.getMonoMass();
                    }
                    std::string new_spec_text = "";

                    if (mass > 0) {
                        new_spec_text = mod_spec.substr(0, mod_spec.length() - unimod.length()) + "+" + std::to_string(mass);
                    } else {
                        new_spec_text = mod_spec.substr(0, mod_spec.length() - unimod.length()) + std::to_string(mass);
                    }

                    if (it == tokens.end() - 1) {
                        actual_data = actual_data + new_spec_text;
                    } else {
                        actual_data = actual_data + new_spec_text + ",";
                    }
                    std::string unimod_title = mod.getTitle();
                    std::string _pb_data = "UNIMOD, UNIMOD:" + std::to_string(unimod_id) + ", " + unimod_title  + ", " + std::to_string(mass);
                    pb_data.push_back(_pb_data);
                }
            } else {
                if (it == tokens.end() - 1) {
                    actual_data = actual_data + mod_spec;
                } else {
                    actual_data = actual_data + mod_spec + ",";
                }
                size_t pos = mod_spec.find_first_of("+-");
                if (pos != std::string::npos) {
                    mod_spec = mod_spec.substr(pos);
                    if (!mod_spec.empty()) {
                        char firstChar = mod_spec[0];
                        if(firstChar == '+') {
                            mod_spec.erase(0, 1);  // Remove the plus sign
                        }
                    }
                }
                std::string _pb_data = "CHEMMOD, CHEMMOD:" + mod_spec + ", unknown modification, ";
                pb_data.push_back(_pb_data);
            }
        }
        response.actual_data = actual_data;
        response.pb_data = pb_data;
        return response;
    }

    /**
     * @brief Parses the modspec string from the command line and loads them into
     * Protobuf collections. Loads all mod deltas into _unique_delta vector _without deduplication_.
     *
     * @param spec_text - comma-separated modspecs from the command line
     * @param mod_table - type of the table to parse - regular or one of terminal mods
     * @return true - if completed successfully
     * @return false otherwise
     */
    bool Parse(const char* input_spec_text, MODS_SPEC_TYPE_T mod_table = MOD_SPEC) {
        ParserResponse parse_response = parse_unimod(input_spec_text);
        string tmp = parse_response.actual_data;
        const char* spec_text = tmp.c_str();

        pb::ModTable* pb_mod_table_ptr = NULL;
        switch (mod_table) {
            case MOD_SPEC:
                pb_mod_table_ptr = &pb_mod_table_;
                break;
            case NTPEP:
                pb_mod_table_ptr = &pb_ntpep_mod_table_;
                break;
            case CTPEP:
                pb_mod_table_ptr = &pb_ctpep_mod_table_;
                break;
            case NTPRO:
                pb_mod_table_ptr = &pb_ntpro_mod_table_;
                break;
            case CTPRO:
                pb_mod_table_ptr = &pb_ctpro_mod_table_;
                break;
        }
        if (pb_mod_table_ptr == NULL)
            return false;

        int pos = 0;
        int pb_data_index = 0;
        while (true) {
            char c;
            int next_pos = -1;
            sscanf(spec_text + pos, " %c%n", &c, &next_pos);
            if (next_pos == -1)
                return Error(spec_text, pos, "Expected modification specification.");
            pos += next_pos - 1;
            unsigned int limit = 0;
            if (c >= '1' && c <= '9') {
                sscanf(spec_text + pos, "%u%n", &limit, &next_pos);
                if (limit == UINT_MAX) {
                    return Error(spec_text, pos, "Limit too big.");
                }
                pos += next_pos;
            }
            int aa_len = -1, plus_pos = -1, delta_pos = -1, end_pos = -1;

            if (mod_table == MOD_SPEC) {
                sscanf(spec_text + pos, "%*[ACDEFGHIJKLMNOPQRSTUVWY]%n%n%*[+-]%n%*[0-9.]%n",
                       &aa_len, &plus_pos, &delta_pos, &end_pos);
            } else {
                sscanf(spec_text + pos, "%*[ACDEFGHIJKLMNOPQRSTUVWYX]%n%n%*[+-]%n%*[0-9.]%n",
                       &aa_len, &plus_pos, &delta_pos, &end_pos);
            }

            if (aa_len == -1)
                return Error(spec_text, pos, "Expected amino acid symbol.");
            assert(plus_pos != -1);
            if (delta_pos == -1)
                return Error(spec_text, pos + plus_pos, "Expected '+' 'or' - and modification amount.");
            if (end_pos == -1)
                return Error(spec_text, pos + delta_pos, "Expected modification amount.");
            if ((limit == 0) && (aa_len != 1))
                return Error(spec_text, pos,
                             "Static modifications must be specified "
                             "for one amino acid at a time.");
            int confirm_end_pos = -1;
            double delta;
            sscanf(spec_text + pos + delta_pos, "%lg%n", &delta, &confirm_end_pos);
            if (delta_pos + confirm_end_pos != end_pos)
                return Error(spec_text, pos + delta_pos,
                             "Cannot parse modification "
                             "amount.");
            if (*(spec_text + pos + plus_pos) == '-')
                delta *= -1;

            pb::Modification* mod;

            if (limit > 1 && pb_mod_table_ptr != &pb_mod_table_)
                limit = 1;

            if (limit == 0) {
                mod = pb_mod_table_ptr->add_static_mod();
            } else {
                mod = pb_mod_table_ptr->add_variable_mod();
                mod->set_max_count(limit);
                unique_delta_.push_back(delta);
                max_counts_.push_back(limit);
            }
            mod->set_amino_acids(string(spec_text + pos, aa_len));
            mod->set_delta(delta);

            mod->set_name(parse_response.pb_data.at(pb_data_index));
            pos += end_pos;
            if (spec_text[pos] == '\0')
                break;

            if (spec_text[pos] == ',')
                ++pos;
            ++pb_data_index;
        }
        return true;
    }
    void ClearTables() {
        pb_mod_table_.Clear();
        pb_ntpep_mod_table_.Clear();
        pb_ctpep_mod_table_.Clear();
        pb_ntpro_mod_table_.Clear();
        pb_ctpro_mod_table_.Clear();
    }

    int NumPoss(char aa, MODS_SPEC_TYPE_T mod_table = MOD_SPEC) const {
        switch (mod_table) {
            case MOD_SPEC:
                return possibles_[aa].size();
            case NTPEP:
                return possibles_ntpe_[aa].size();
            case CTPEP:
                return possibles_ctpe_[aa].size();
            case NTPRO:
                return possibles_ntpr_[aa].size();
            case CTPRO:
                return possibles_ctpr_[aa].size();
        }
        return 0;
    }

    int PossMaxCt(char aa, int index, MODS_SPEC_TYPE_T mod_table = MOD_SPEC) const {
        switch (mod_table) {
            case MOD_SPEC:
                return possibles_[aa][index].second;
            case NTPEP:
                return possibles_ntpe_[aa][index].second;
            case CTPEP:
                return possibles_ctpe_[aa][index].second;
            case NTPRO:
                return possibles_ntpr_[aa][index].second;
            case CTPRO:
                return possibles_ctpr_[aa][index].second;
        }
        return 0;
    }

    int PossDeltIx(char aa, int index, MODS_SPEC_TYPE_T mod_table = MOD_SPEC) const {
        switch (mod_table) {
            case MOD_SPEC:
                return possibles_[aa][index].first;
            case NTPEP:
                return possibles_ntpe_[aa][index].first;
            case CTPEP:
                return possibles_ctpe_[aa][index].first;
            case NTPRO:
                return possibles_ntpr_[aa][index].first;
            case CTPRO:
                return possibles_ctpr_[aa][index].first;
        }
        return 0;
    }

    double PossDelta(int possDeltIx) const {
        return unique_delta_[possDeltIx];
    }

    double PossDelta(char aa, int index) const {
        return unique_delta_[PossDeltIx(aa, index)];
    }

#if 0
  const IntPairVec* Possibles(char aa) const {
    return &possibles_[aa];
  }
#endif

    /**
     * @brief Deduplicates the unique_deltas_ collection after reserving a copy
     * in the original_deltas_. Loads possibles and unique_deltas for each mod type.
     *
     * @return true
     * @return false
     */
    bool SerializeUniqueDeltas() {
        if (unique_delta_.size() == 0)
            return (0);
        original_deltas_.resize(unique_delta_.size());
        copy(unique_delta_.begin(), unique_delta_.end(), original_deltas_.begin());
        sort(unique_delta_.begin(), unique_delta_.end());
        unique_delta_.resize(unique(unique_delta_.begin(), unique_delta_.end()) - unique_delta_.begin());
        coder_.Init(unique_delta_.size());

        // The following these need to be in exactly the same order as parsed in the TideIndexApplication.cpp
        Init(pb_mod_table_);
        Init(pb_ctpep_mod_table_);
        Init(pb_ntpep_mod_table_);
        Init(pb_ctpro_mod_table_);
        Init(pb_ntpro_mod_table_);
        SerializeUniqueDeltas(&pb_mod_table_);
        SerializeUniqueDeltas(&pb_ctpep_mod_table_);
        SerializeUniqueDeltas(&pb_ntpep_mod_table_);
        SerializeUniqueDeltas(&pb_ctpro_mod_table_);
        SerializeUniqueDeltas(&pb_ntpro_mod_table_);
        return (1);
    }

    /**
     * @brief Copies unique_deltas_ values into the protobuf mods collection. As a result each protobuf collection
     * has all the deltas of the transforms specified on the command line.
     *
     * @param pb_mod_table - protobuf table to load
     * @return true if the load is successful
     * @return false if the protobuf collection already contains unique deltas and they are different from that we have globally.
     */
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
    const pb::ModTable* ParsedNtpepModTable() const { return &pb_ntpep_mod_table_; }
    const pb::ModTable* ParsedCtpepModTable() const { return &pb_ctpep_mod_table_; }
    const pb::ModTable* ParsedNtproModTable() const { return &pb_ntpro_mod_table_; }
    const pb::ModTable* ParsedCtproModTable() const { return &pb_ctpro_mod_table_; }
    const vector<int>* MaxCounts() const { return &max_counts_; }
    const vector<double>* OriginalDeltas() const { return &original_deltas_; }

    int EncodeMod(int aa_index, int unique_delta_index) const {
        return coder_.EncodeMod(aa_index, unique_delta_index);
    }

    void DecodeMod(int code, int* aa_index, int* unique_delta_index) const {
        coder_.DecodeMod(code, aa_index, unique_delta_index);
    }

    void Show() {
        SHOW_ONE(unique_delta_.size());
        SHOW_ALL(unique_delta_);
        SHOW_ONE(max_counts_.size());
        SHOW_ALL(max_counts_);
        const char* aa = "ACDEFGHIJKLMNOPQRSTUVWYX";
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
    int Unique_delta_size() { return unique_delta_.size(); }

   private:
    bool Error(const char* spec_text, int err_pos, const char* msg) {
        cerr << "Error: couldn't parse modification specification:\n"
             << spec_text
             << "\n"
             << string(err_pos, ' ') << "^\n"
             << msg << "\n";
        pb_mod_table_.Clear();
        return false;
    }

    /**
     * @brief possibles_ is an array of vectors where for each AA (indexed by the AA character as array position)
     * there is a list of (position in _unique_delta, global mod index) pairs.
     * global mod index is the number of the given mod in the list of mods specified on the command line.
     */
    typedef vector<pair<int, int> > IntPairVec;
    IntPairVec possibles_[256];       // unique_delta_, max_count_
    IntPairVec possibles_ctpe_[256];  // unique_delta_, max_count_  cterminal peptide
    IntPairVec possibles_ntpe_[256];  // unique_delta_, max_count_  nterminal peptide
    IntPairVec possibles_ctpr_[256];  // unique_delta_, max_count_  cterminal protein
    IntPairVec possibles_ntpr_[256];  // unique_delta_, max_count_  nterminal protein
    vector<double> unique_delta_, original_deltas_;
    vector<int> max_counts_;
    int offset_;
    ModCoder coder_;
    pb::ModTable pb_mod_table_;        // modification table for regular amino acid modifications
    pb::ModTable pb_ntpep_mod_table_;  // modification table for n-terminal peptide modifications
    pb::ModTable pb_ctpep_mod_table_;  // modification table for c-terminal peptide modifications
    pb::ModTable pb_ntpro_mod_table_;  // modification table for n-terminal protein modifications
    pb::ModTable pb_ctpro_mod_table_;  // modification table for c-terminal protein modifications
};

#endif
