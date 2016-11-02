#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include <iomanip>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

class StringUtils {
 public:
  // Convert from a string
  template<typename T>
  static T FromString(const std::string& s) {
    T out;
    if (!TryFromString(s, &out)) {
      throw std::runtime_error("Could not convert string '" + s + "'");
    }
    return out;
  }

  //Convert from a string, returning false if conversion fails
  template<typename T>
  static bool TryFromString(const std::string& s, T* out) {
    std::stringstream ss(s);
    ss >> *out;
    return ss.eof() && !ss.fail();
  }
  static bool TryFromString(const std::string& s, std::string* out) {
    *out = s;
    return true;
  }

  // Convert to a string
  // Description added by Andy Lin
  // If fixedFloat is true, then the enum (most likely defined in objects.h)
  // will be printed out with the precision stored in variable 'decimals'
  // If fixedFloat is false, then enum (obj) will be printed out in scientic notation 
  // If decimals is not given, then enum (obj) will be printed with a precision of 8
  template<typename T>
  static std::string ToString(const T& obj, int decimals = -1, bool fixedFloat = true) {
    std::stringstream converter;
    if (decimals >= 0) {
      converter << std::fixed << std::setprecision(decimals);
      if (fixedFloat) {
        converter << std::fixed;
      } else {
        converter.unsetf(std::ios_base::floatfield);
      }
    } else {
      converter << std::setprecision(8);
    }
    converter << obj;
    return converter.str();
  }

  // Joins a vector of strings into a single string separated by a delimiter
  template<typename T>
  static std::string Join(const T values, const char delimiter ='\0') {
    std::stringstream ss;
    for (typename T::const_iterator i = values.begin(); i != values.end(); i++) {
      if (i != values.begin() && delimiter != '\0') {
        ss << delimiter;
      }
      ss << *i;
    }
    return ss.str();
  }

  // Split a string on a delimiter
  static std::vector<std::string> Split(const std::string& s, char delimiter);
  static std::vector<std::string> Split(const std::string& s, const std::string& delimiter);

  template<typename T>
  static std::vector<T> Split(const std::string& s, char delimiter) {
    std::vector<T> tokens;
    std::string::const_iterator from = s.begin();
    for (std::string::const_iterator i = from; i != s.end(); ++i) {
      if (*i == delimiter) {
        tokens.push_back(FromString<T>(std::string(from, i)));
        from = i + 1;
      }
    }
    tokens.push_back(FromString<T>(std::string(from, s.end())));
    return tokens;
  }

  template<typename T>
  static std::vector<T> Split(const std::string& s, const std::string& delimiter) {
    std::vector<T> tokens;
    size_t from = 0;
    size_t find;
    while ((find = s.find(delimiter, from)) != std::string::npos) {
      tokens.push_back(FromString<T>(s.substr(from, find - from)));
      from = find + delimiter.length();
    }
    tokens.push_back(FromString<T>(s.substr(from)));
    return tokens;
  }

  static std::vector<std::string> Fields(const std::string& s);

  template<typename T>
  static std::vector<T> Fields(const std::string& s) {
    std::vector<T> fields;
    std::string current;
    for (std::string::const_iterator i = s.begin(); i != s.end(); i++) {
      if (*i == ' ' || *i == '\t') {
        if (!current.empty()) {
          fields.push_back(FromString<T>(current));
          current.clear();
        }
      } else {
        current.push_back(*i);
      }
    }
    if (!current.empty()) {
      fields.push_back(FromString<T>(current));
    }
    return fields;
  }

  // Convert a string to lowercase
  static std::string ToLower(std::string s);

  // Convert a string to uppercase
  static std::string ToUpper(std::string s);

  // Compare two strings without case sensitivity
  static bool IEquals(const std::string& x, const std::string& y);

  // Replace all instances of a substring with a new substring
  static std::string Replace(
    std::string s, const std::string& oldSubstring, const std::string& newSubstring);

  // Return whether a string begins with a substring
  static bool StartsWith(const std::string& s, const std::string& substring);

  // Return whether a string starts with a substring (case insensitive)
  static bool IStartsWith(const std::string& s, const std::string& substring);

  // Return whether a string ends with a substring
  static bool EndsWith(const std::string& s, const std::string& substring);

  // Return whether a string ends with a substring (case insensitive)
  static bool IEndsWith(const std::string& s, const std::string& substring);

  // Trim whitespace from the beginning and end of a string
  static std::string Trim(std::string s);

  // Trim whitespace from the beginning of a string
  static std::string LTrim(std::string s);

  // Trim whitespace from the end of a string
  static std::string RTrim(std::string s);

  // Return whether a string is numeric
  static bool IsNumeric(
    const std::string& s, bool allowNegative = true, bool allowDecimal = true);

  // Break a string into lines limited by length
  static std::string LineFormat(std::string s, unsigned limit, unsigned indentSize = 0);

 private:
  static const char* WHITESPACE_CHARS;

  StringUtils();
  ~StringUtils();
};

#endif

