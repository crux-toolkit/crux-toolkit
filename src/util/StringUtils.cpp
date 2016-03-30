#include "StringUtils.h"

#include "boost/algorithm/string.hpp"

using namespace std;

const char* StringUtils::WHITESPACE_CHARS = " \t\n\v\f\r";

vector<string> StringUtils::Split(const string& s, char delimiter) {
  return Split<string>(s, delimiter);
}

vector<string> StringUtils::Split(const string& s, const string& delimiter) {
  return Split<string>(s, delimiter);
}

vector<string> StringUtils::Fields(const string& s) {
  return Fields<string>(s);
}

string StringUtils::ToLower(string s) {
  boost::to_lower(s);
  return s;
}

string StringUtils::ToUpper(string s) {
  boost::to_upper(s);
  return s;
}

bool StringUtils::IEquals(const string& x, const string& y) {
  return boost::iequals(x, y);
}

string StringUtils::Replace(string s, const string& oldSubstring, const string& newSubstring) {
  size_t find = 0;
  while ((find = s.find(oldSubstring, find)) != string::npos) {
    s.replace(find, oldSubstring.length(), newSubstring);
    find += newSubstring.length();
  }
  return s;
}

bool StringUtils::StartsWith(const string& s, const string& substring) {
  return substring.length() <= s.length()
    ? s.compare(0, substring.length(), substring) == 0
    : false;
}

bool StringUtils::IStartsWith(const string& s, const string& substring) {
  return StartsWith(ToLower(s), ToLower(substring));
}

bool StringUtils::EndsWith(const string& s, const string& substring) {
  return substring.length() <= s.length()
    ? s.compare(s.length() - substring.length(), substring.length(), substring) == 0
    : false;
}

bool StringUtils::IEndsWith(const string& s, const string& substring) {
  return EndsWith(ToLower(s), ToLower(substring));
}

string StringUtils::Trim(string s) {
  boost::trim(s);
  return s;
}

string StringUtils::LTrim(string s) {
  boost::trim_left(s);
  return s;
}

string StringUtils::RTrim(string s) {
  boost::trim_right(s);
  return s;
}

bool StringUtils::IsNumeric(const string& s, bool allowNegative, bool allowDecimal) {
  if (s.empty()) {
    return false;
  }

  string::const_iterator i = s.begin();
  if (allowNegative && *i == '-') {
    ++i;
  }

  bool foundDecimal = false;

  for (; i != s.end(); ++i) {
    if (allowDecimal && *i == '.') {
      if (foundDecimal) {
        // Found two decimals
        return false;
      }
      foundDecimal = true;
    } else if (!isdigit(*i)) {
      // Found non-digit character
      return false;
    }
  }

  // Decimal can't be last character
  return !foundDecimal || *(--i) == '.';
}

string StringUtils::LineFormat(string s, unsigned limit, unsigned indentSize) {
  if (indentSize >= limit) {
    throw runtime_error("Limit must be greater than indent size");
  }
  stringstream lines;
  string indent(indentSize, ' ');
  unsigned contentLimit = limit - indentSize;

  while (s.length() > contentLimit) {
    size_t lineEnd;
    if (isspace(s[contentLimit])) {
      lineEnd = s.find_last_not_of(WHITESPACE_CHARS, contentLimit);
    } else {
      lineEnd = s.find_last_of(WHITESPACE_CHARS, contentLimit);
      lineEnd = (lineEnd == string::npos)
        ? lineEnd = contentLimit - 1
        : lineEnd = s.find_last_not_of(WHITESPACE_CHARS, lineEnd);
    }
    if (lineEnd == string::npos) {
      s = LTrim(s);
      continue;
    }
    lines << indent << s.substr(0, lineEnd + 1) << endl;
    s = LTrim(s.substr(lineEnd + 2));
  }
  lines << indent << s;
  return lines.str();
}

StringUtils::StringUtils() {}
StringUtils::~StringUtils() {}

