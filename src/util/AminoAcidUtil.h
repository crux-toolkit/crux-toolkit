#ifndef AMINOACIDUTIL_H
#define AMINOACIDUTIL_H

#include <string>

const double CYSTEINE_DEFAULT = 57.021464;

class AminoAcidUtil {
 public:
  static std::string GetName(char c);
  static double GetMass(char c, bool monoisotopic = true);
 private:
  AminoAcidUtil();
  ~AminoAcidUtil();
};

#endif

