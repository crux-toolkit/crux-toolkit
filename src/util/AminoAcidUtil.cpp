#include "AminoAcidUtil.h"

#include <limits>
#include <stdexcept>

using namespace std;

string AminoAcidUtil::GetName(char c) {
  switch (c) {
    case 'A': return "alanine";
    case 'B': return "";
    case 'C': return "cysteine";
    case 'D': return "aspartic acid";
    case 'E': return "glutamic acid";
    case 'F': return "phenylalanine";
    case 'G': return "glycine";
    case 'H': return "histidine";
    case 'I': return "isoleucine";
    case 'J': return "";
    case 'K': return "lysine";
    case 'L': return "leucine";
    case 'M': return "methionine";
    case 'N': return "asparagine";
    case 'O': return "ornithine";
    case 'P': return "proline";
    case 'Q': return "glutamine";
    case 'R': return "arginine";
    case 'S': return "serine";
    case 'T': return "threonine";
    case 'U': return "selenocysteine";
    case 'V': return "valine";
    case 'W': return "tryptophan";
    case 'X': return "";
    case 'Y': return "tyrosine";
    case 'Z': return "";
  }
  throw runtime_error("'" + string(1, c) + "' is not a valid amino acid");
}

double AminoAcidUtil::GetMass(char c, bool monoisotopic) {
  if (monoisotopic) {
    // Monoisotopic masses
    switch (c) {
      case 'A': return  71.03711;
      case 'B': return numeric_limits<double>::quiet_NaN();
      case 'C': return 103.00919;
      case 'D': return 115.02694;
      case 'E': return 129.04259;
      case 'F': return 147.06841;
      case 'G': return  57.02146;
      case 'H': return 137.05891;
      case 'I': return 113.08406;
      case 'J': return 113.08406;
      case 'K': return 128.09496;
      case 'L': return 113.08406;
      case 'M': return 131.04049;
      case 'N': return 114.04293;
      case 'O': return 114.07931;
      case 'P': return  97.05276;
      case 'Q': return 128.05858;
      case 'R': return 156.10111;
      case 'S': return  87.03203;
      case 'T': return 101.04768;
      case 'U': return 150.04344;
      case 'V': return  99.06841;
      case 'W': return 186.07931;
      case 'X': return numeric_limits<double>::quiet_NaN();
      case 'Y': return 163.06333;
      case 'Z': return numeric_limits<double>::quiet_NaN();
    }
  } else {
    // Average masses
    switch (c) {
      case 'A': return  71.0788;
      case 'B': return numeric_limits<double>::quiet_NaN();
      case 'C': return 103.1388;
      case 'D': return 115.0886;
      case 'E': return 129.1155;
      case 'F': return 147.1766;
      case 'G': return  57.0519;
      case 'H': return 137.1411;
      case 'I': return 113.1594;
      case 'J': return 113.1594;
      case 'K': return 128.1741;
      case 'L': return 113.1594;
      case 'M': return 131.1926;
      case 'N': return 114.1038;
      case 'O': return 114.1472;
      case 'P': return  97.1167;
      case 'Q': return 128.1307;
      case 'R': return 156.1875;
      case 'S': return  87.0782;
      case 'T': return 101.1051;
      case 'U': return 150.0388;
      case 'V': return  99.1326;
      case 'W': return 186.2132;
      case 'X': return numeric_limits<double>::quiet_NaN();
      case 'Y': return 163.1760;
      case 'Z': return numeric_limits<double>::quiet_NaN();
    }
  }
  throw runtime_error("'" + string(1, c) + "' is not a valid amino acid");
}

AminoAcidUtil::AminoAcidUtil() {}
AminoAcidUtil::~AminoAcidUtil() {}

