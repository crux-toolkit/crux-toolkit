#include "AminoAcidUtil.h"

#include <stdexcept>

using namespace std;

string AminoAcidUtil::GetName(char c) {
  if (c < 'A' || c > 'Z') {
    throw runtime_error("'" + string(1, c) + "' is not a valid amino acid");
  }
  switch (c) {
    case 'A': return "alanine";
    case 'C': return "cysteine";
    case 'D': return "aspartic acid";
    case 'E': return "glutamic acid";
    case 'F': return "phenylalanine";
    case 'G': return "glycine";
    case 'H': return "histidine";
    case 'I': return "isoleucine";
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
    case 'Y': return "tyrosine";
    default : return "";
  }
}

AminoAcidUtil::AminoAcidUtil() {}
AminoAcidUtil::~AminoAcidUtil() {}

