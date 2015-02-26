#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include "boost/lexical_cast.hpp"

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
    try {
      return boost::lexical_cast<T>(s);
    } catch (...) {
      throw std::runtime_error("Could not convert string '" + s + "'");
    }
  }

  // Convert to a string
  template<typename T>
  static std::string ToString(const T& obj, int decimals = -1) {
    std::stringstream converter;
    if (decimals >= 0) {
      converter << std::fixed << std::setprecision(decimals);
    }
    converter << obj;
    return converter.str();
  }

  // Joins a vector of strings into a single string separated by a delimiter
  template<typename T>
  static std::string Join(const std::vector<T>& v, const char delimiter = '\0') {
    std::stringstream ss;
    for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); i++) {
      if (i != v.begin() && delimiter != '\0') {
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

  // Return whether a string ends with a substring
  static bool EndsWith(const std::string& s, const std::string& substring);

  // Trim whitespace from the beginning and end of a string
  static std::string Trim(const std::string& s);

  // Trim whitespace from the beginning of a string
  static std::string LTrim(const std::string& s);

  // Trim whitespace from the end of a string
  static std::string RTrim(const std::string& s);

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

