#ifndef AMINOACIDUTIL_H
#define AMINOACIDUTIL_H

#include <string>

const double CYSTEINE_DEFAULT = 57.021464;

class AminoAcidUtil {
 public:
  static std::string GetName(char c);
 private:
  AminoAcidUtil();
  ~AminoAcidUtil();
};

#endif

