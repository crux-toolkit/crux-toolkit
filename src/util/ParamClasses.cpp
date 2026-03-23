#include "Params.h"
#include "StringUtils.h"

#include <algorithm>

using namespace std;

// ***** Parameter classes ***** //
//
// Param (base class)
//
Param::Param(const string& name,
             const string& usage,
             const string& fileNotes,
             bool visible)
  : name_(name), usage_(usage), fileNotes_(fileNotes), visible_(visible) {}
Param::~Param() {}
string Param::GetName() const { return name_; }
string Param::GetUsage() const { return usage_; }
string Param::GetFileNotes() const { return fileNotes_; }
bool Param::IsVisible() const { return visible_; }
bool Param::IsArgument() const { return false; }
void Param::ThrowIfInvalid() const {}
string Param::GetParamFileString(bool defaultValue) const {
  vector<string> lines =
    StringUtils::Split(Params::ProcessHtmlDocTags(usage_), '\n');
  vector<string> noteLines =
    StringUtils::Split(Params::ProcessHtmlDocTags(fileNotes_), '\n');
  lines.insert(lines.end(), noteLines.begin(), noteLines.end());
  stringstream ss;
  for (vector<string>::const_iterator i = lines.begin(); i != lines.end(); i++) {
    vector<string> formatted = StringUtils::Split(StringUtils::LineFormat(*i, 78), '\n');
    for (vector<string>::const_iterator j = formatted.begin(); j != formatted.end(); j++) {
      ss << "# " << *j << endl;
    }
  }
  ss << name_ << '=' << (defaultValue ? GetStringDefault() : GetString()) << endl;
  return ss.str();
}
void Param::Set(bool value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from bool");
}
void Param::Set(int value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from int");
}
void Param::Set(double value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from double");
}
void Param::Set(const char* value) {
  Set(string(value));
}
void Param::Set(const string& value) {
  throw runtime_error("Cannot set value of '" + name_ + "' from string");
}
//
// BoolParam
//
BoolParam::BoolParam(const string& name,
                     const string& usage,
                     const string& fileNotes,
                     bool visible,
                     bool value)
  : Param(name, usage, fileNotes, visible), value_(value), original_(value) {}
string BoolParam::GetAcceptedValues() const { return "T|F"; }
bool BoolParam::IsDefault() const { return value_ == original_; }
bool BoolParam::GetBool() const { return value_; }
int BoolParam::GetInt() const { return IntParam::From(value_); }
double BoolParam::GetDouble() const { return DoubleParam::From(value_); }
string BoolParam::GetString() const { return StringParam::From(value_); }
bool BoolParam::GetBoolDefault() const { return original_; }
int BoolParam::GetIntDefault() const { return IntParam::From(original_); }
double BoolParam::GetDoubleDefault() const { return DoubleParam::From(original_); }
string BoolParam::GetStringDefault() const { return StringParam::From(original_); }
void BoolParam::Set(bool value) { value_ = value; }
void BoolParam::Set(int value) { value_ = From(value); }
void BoolParam::Set(double value) { value_ = From(value); }
void BoolParam::Set(const string& value) {
  try {
    value_ = From(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected boolean)");
  }
}
bool BoolParam::From(int i) { return i != 0; }
bool BoolParam::From(double d) { return d != 0; }
bool BoolParam::From(string s) {
  s = StringUtils::ToLower(s);
  if (s == "t" || s == "true") {
    return true;
  } else if (s == "f" || s == "false") {
    return false;
  }
  throw runtime_error("Cannot convert '" + s + "' to boolean");
}
//
// IntParam
//
IntParam::IntParam(const string& name,
                   const string& usage,
                   const string& fileNotes,
                   bool visible,
                   int value,
                   int min,
                   int max)
  : Param(name, usage, fileNotes, visible),
    value_(value), original_(value), min_(min), max_(max) {}
void IntParam::ThrowIfInvalid() const {
  if (value_ < min_ || value_ > max_) {
    throw runtime_error("Value of '" + name_ + "' must be between " +
                        StringUtils::ToString(min_) + " and " +
                        StringUtils::ToString(max_));
  }
}
string IntParam::GetAcceptedValues() const { return "<integer>"; }
bool IntParam::IsDefault() const { return value_ == original_; }
bool IntParam::GetBool() const { return BoolParam::From(value_); }
int IntParam::GetInt() const { return value_; }
double IntParam::GetDouble() const { return DoubleParam::From(value_); }
string IntParam::GetString() const { return StringParam::From(value_); }
bool IntParam::GetBoolDefault() const { return BoolParam::From(original_); }
int IntParam::GetIntDefault() const { return original_; }
double IntParam::GetDoubleDefault() const { return DoubleParam::From(original_); }
string IntParam::GetStringDefault() const { return StringParam::From(original_); }
void IntParam::Set(bool value) { value_ = From(value); }
void IntParam::Set(int value) { value_ = value; }
void IntParam::Set(double value) { value_ = From(value); }
void IntParam::Set(const string& value) {
  try {
    value_ = From(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected int)");
  }
}
int IntParam::From(bool b) { return b ? 1 : 0; }
int IntParam::From(double d) { return (int)d; }
int IntParam::From(const string& s) { return StringUtils::FromString<int>(s); }
//
// DoubleParam
//
DoubleParam::DoubleParam(const string& name,
                         const string& usage,
                         const string& fileNotes,
                         bool visible,
                         double value,
                         double min,
                         double max)
  : Param(name, usage, fileNotes, visible),
    value_(value), original_(value), min_(min), max_(max) {}
void DoubleParam::ThrowIfInvalid() const {
  if (value_ < min_ || value_ > max_) {
    throw runtime_error("Value of '" + name_ + "' must be between " +
                        StringUtils::ToString(min_) + " and " +
                        StringUtils::ToString(max_));
  }
}
string DoubleParam::GetAcceptedValues() const { return "<float>"; }
bool DoubleParam::IsDefault() const { return value_ == original_; }
bool DoubleParam::GetBool() const { return BoolParam::From(value_); }
int DoubleParam::GetInt() const { return IntParam::From(value_); }
double DoubleParam::GetDouble() const { return value_; }
string DoubleParam::GetString() const { return StringParam::From(value_); }
bool DoubleParam::GetBoolDefault() const { return BoolParam::From(original_); }
int DoubleParam::GetIntDefault() const { return IntParam::From(original_); }
double DoubleParam::GetDoubleDefault() const { return original_; }
string DoubleParam::GetStringDefault() const { return StringParam::From(original_); }
void DoubleParam::Set(bool value) { value_ = From(value); }
void DoubleParam::Set(int value) { value_ = From(value); }
void DoubleParam::Set(double value) { value_ = value; }
void DoubleParam::Set(const string& value) {
  try {
    value_ = From(value);
  } catch (...) {
    throw runtime_error("Invalid value for '" + name_ + "': " + "'" + value + "' "
                        "(expected float)");
  }
}
double DoubleParam::From(bool b) { return b ? 1 : 0; }
double DoubleParam::From(int i) { return (double)i; }
double DoubleParam::From(const string& s) { return StringUtils::FromString<double>(s); }
//
// StringParam
//
StringParam::StringParam(const string& name,
                         const string& usage,
                         const string& fileNotes,
                         bool visible,
                         const string& value,
                         const vector<string>& validValues)
  : Param(name, usage, fileNotes, visible),
    original_(value), validValues_(validValues) {
  Set(value);
}
void StringParam::ThrowIfInvalid() const {
  if (!validValues_.empty() &&
      find(validValues_.begin(), validValues_.end(), value_) == validValues_.end()) {
    throw runtime_error("Invalid value for '" + name_ + "'; must be one of <" +
                        StringUtils::Join(validValues_, '|') + ">");
  }
}
string StringParam::GetAcceptedValues() const {
  return validValues_.empty() ? "<string>" : StringUtils::Join(validValues_, '|');
}
bool StringParam::IsDefault() const { return value_ == original_; }
bool StringParam::GetBool() const { return BoolParam::From(value_); }
int StringParam::GetInt() const { return IntParam::From(value_); }
double StringParam::GetDouble() const { return DoubleParam::From(value_); }
string StringParam::GetString() const { return value_; }
bool StringParam::GetBoolDefault() const { return BoolParam::From(original_); }
int StringParam::GetIntDefault() const { return IntParam::From(original_); }
double StringParam::GetDoubleDefault() const { return DoubleParam::From(original_); }
string StringParam::GetStringDefault() const { return original_; }
void StringParam::Set(bool value) { value_ = From(value); }
void StringParam::Set(int value) { value_ = From(value); }
void StringParam::Set(double value) { value_ = From(value); }
void StringParam::Set(const string& value) {
  value_ = value != "__NULL_STR" ? value : "";
}
string StringParam::From(bool b) { return b ? "true" : "false"; }
string StringParam::From(int i) { return StringUtils::ToString(i); }
string StringParam::From(double d) { return StringUtils::ToString(d); }
//
// ArgParam
//
ArgParam::ArgParam(const string& name, const string& usage)
  : Param(name, usage, "", false), values_(vector<string>()) {}
string ArgParam::GetAcceptedValues() const { return "<string>"; }
bool ArgParam::IsArgument() const { return true; }
bool ArgParam::IsDefault() const { return false; }
bool ArgParam::GetBool() const { return BoolParam::From(GetString()); }
int ArgParam::GetInt() const { return StringUtils::FromString<int>(GetString()); }
double ArgParam::GetDouble() const { return StringUtils::FromString<double>(GetString()); }
string ArgParam::GetString() const {
  if (values_.empty()) {
    throw runtime_error("No value for argument '" + name_ + "'");
  }
  return values_.front();
}
bool ArgParam::GetBoolDefault() const { return false; }
int ArgParam::GetIntDefault() const { return 0; }
double ArgParam::GetDoubleDefault() const { return 0; }
string ArgParam::GetStringDefault() const { return ""; }
const vector<string>& ArgParam::GetStrings() const { return values_; }
void ArgParam::Set(bool value) { values_ = vector<string>(1, StringParam::From(value)); }
void ArgParam::Set(int value) { values_ = vector<string>(1, StringParam::From(value)); }
void ArgParam::Set(double value) { values_ = vector<string>(1, StringParam::From(value)); }
void ArgParam::Set(const string& value) { values_ = vector<string>(1, value); }
void ArgParam::AddValue(const string& value) { values_.push_back(value); }

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

