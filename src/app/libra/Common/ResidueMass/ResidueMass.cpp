#include "ResidueMass.h"


//
// stuff to load user-specified residue masses
// Copyright (C) Insilicos LLC 2006-2010 All Rights Reserved
//
#include "Common/util.h"
#include "Common/tpp_hashmap.h" // hashmap isn't in std space for some compilers
static void load_table(); // declared below
typedef struct {double m_monomass;double m_mass;bool m_nterm;bool m_cterm;} resinfo_t;
class cResidue {
public:
   cResidue()
   {};
   cResidue(const cResidue &rhs) :
      m_residue(rhs.m_residue),
      m_masses(rhs.m_masses)
   {};
   ~cResidue()
   {};

   std::string m_residue;
   std::vector<resinfo_t> m_masses;
};

// store as a hashmap of title:residues
typedef TPP_STDSTRING_HASH_MULTIMAP(cResidue) rmap;
static rmap residues; // list of residues

//
// end declarations for stuff for user defined mods (see below for implementation)
//

double ResidueMass::monoisotopic_nterm_mass_ = (double) 1.0078250321;
double ResidueMass::average_nterm_mass_ = (double) 1.0078250321;
double ResidueMass::monoisotopic_cterm_mass_ = (double) 17.0027396542;
double ResidueMass::average_cterm_mass_ = (double) 17.0027396542;

double ResidueMass::average_masses_[] = {
  (double) 71.0788, // A
  (double) 114.5962, // B
  (double) 103.1388, // C
  (double) 115.0886, // D
  (double) 129.1155, // E
  (double) 147.1766, // F
  (double) 57.0519, // G
  (double) 137.1411, // H
  (double) 113.1594, // I
  (double) 0.0000, // J
  (double) 128.1741, // K
  (double) 113.1594, // L
  (double) 131.1926, // M
  (double) 114.1038, // N
  (double) 114.1472, // O
  (double) 97.1167, // P
  (double) 128.1307, // Q
  (double) 156.1875, // R
  (double) 87.0782, // S
  (double) 101.1051, // T
  (double) 150.0379, // U
  (double) 99.1326, // V
  (double) 186.2132, // W
  (double) 113.1594, // X
  (double) 163.1760, //  Y
  (double) 128.6231  // Z
};

double ResidueMass::monoisotopic_masses_[] = {
  (double) 71.03711, // A
  (double) 114.53494, // B
  (double) 103.00919, // C
  (double) 115.02694, // D
  (double) 129.04259, // E
  (double) 147.06841, // F
  (double) 57.02146, // G
  (double) 137.05891, // H
  (double) 113.08406, // I
  (double) 0.0000, // J
  (double) 128.09496, // K
  (double) 113.08406, // L
  (double) 131.04049, // M
  (double) 114.04293, // N
  (double) 114.07931, // O
  (double) 97.05276, // P
  (double) 128.05858, // Q
  (double) 156.10111, // R
  (double) 87.03203, // S
  (double) 101.04768, // T
  (double) 150.95363, // U
  (double) 99.06841, // V
  (double) 186.07931, // W
  (double) 113.08406, // X
  (double) 163.06333, //  Y
  (double) 128.55059 // Z
};
/*
add_G_Glycine = 0.0000                 ; added to G - avg.  57.0519, mono.  57.02146
add_A_Alanine = 0.0000                 ; added to A - avg.  71.0788, mono.  71.03711
add_S_Serine = 0.0000                  ; added to S - avg.  87.0782, mono.  87.02303
add_P_Proline = 0.0000                 ; added to P - avg.  97.1167, mono.  97.05276
add_V_Valine = 0.0000                  ; added to V - avg.  99.1326, mono.  99.06841
add_T_Threonine = 0.0000               ; added to T - avg. 101.1051, mono. 101.04768
add_C_Cysteine = 57.05                 ; added to C - avg. 103.1388, mono. 103.00919
add_L_Leucine = 0.0000                 ; added to L - avg. 113.1594, mono. 113.08406
add_I_Isoleucine = 0.0000              ; added to I - avg. 113.1594, mono. 113.08406
add_X_LorI = 0.0000                    ; added to X - avg. 113.1594, mono. 113.08406
add_N_Asparagine = 0.0000              ; added to N - avg. 114.1038, mono. 114.04293
add_O_Ornithine = 0.0000               ; added to O - avg. 114.1472, mono  114.07931
add_B_avg_NandD = 0.0000               ; added to B - avg. 114.5962, mono. 114.53494
add_D_Aspartic_Acid = 0.0000           ; added to D - avg. 115.0886, mono. 115.02694
add_Q_Glutamine = 0.0000               ; added to Q - avg. 128.1307, mono. 128.05858
add_K_Lysine = 0.0000                  ; added to K - avg. 128.1741, mono. 128.09496
add_Z_avg_QandE = 0.0000               ; added to Z - avg. 128.6231, mono. 128.55059
add_E_Glutamic_Acid = 0.0000           ; added to E - avg. 129.1155, mono. 129.04259
add_M_Methionine = 0.0000              ; added to M - avg. 131.1926, mono. 131.04049
add_H_Histidine = 0.0000               ; added to H - avg. 137.1411, mono. 137.05891
add_F_Phenyalanine = 0.0000            ; added to F - avg. 147.1766, mono. 147.06841
add_R_Arginine = 0.0000                ; added to R - avg. 156.1875, mono. 156.10111
add_Y_Tyrosine = 0.0000                ; added to Y - avg. 163.1760, mono. 163.06333
add_W_Tryptophan = 0.0000              ; added to W - avg. 186.2132, mono. 186.07931
*/


double ResidueMass::getMass(char res, Boolean monoisotopic) {
  // special cases for termini
  if(res == 'n' || res == '1') {
    if(monoisotopic)
      return monoisotopic_nterm_mass_;
    return average_nterm_mass_;
  }
  else if(res == 'c' || res == '2') {
    if(monoisotopic)
      return monoisotopic_cterm_mass_;
    return average_cterm_mass_;
  }
  else if(res < 'A' || res > 'Z') { // error
	  std::cerr << "WARNING: Trying to compute mass of non-residue: " << res << std::endl;
    return 0;
  }
  if(monoisotopic)
    return getMonoisotopicMass(res);
  return getAverageMass(res);
}


double ResidueMass::getProteinMass(const char* prot, Boolean monoisotopic) {

  double total = getMass('n', monoisotopic) + getMass('c', monoisotopic);
  if(prot != NULL) {
    for(int k = 0; prot[k]; k++)
      total += getMass(prot[k], monoisotopic);
  }
  return total;

}

double ResidueMass::getAverageMass(char res) {
  return average_masses_[res-'A'];
}

double ResidueMass::getMonoisotopicMass(char res) {
  return monoisotopic_masses_[res-'A'];
}

double ResidueMass::getStdModMass(const char* mod, Boolean monoisotopic, char label) {
//
// first check for custom mods
//
load_table(); // in case we didn't already
std::pair<rmap::const_iterator, rmap::const_iterator> p =
  residues.equal_range(std::string(mod)); //  set up iterator
for (rmap::const_iterator i = p.first; i != p.second; ++i) {
  const cResidue &r = (*i).second;
  for (int i=0;r.m_residue[i];i++) {
     if (label==r.m_residue[i]) {
       return monoisotopic?r.m_masses[i].m_monomass:r.m_masses[i].m_mass;
     }
  }
}

// now the usual suspects
if(! strcmp(mod, "+N-formyl-met (Protein)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 161.051;
  else
     return (double) 161.222;
}
if(! strcmp(mod, "13C6-15N2 (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
    return (double) 136.109;
  else
    return (double) 136.115;
}
if(! strcmp(mod, "13C6-15N4 (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
    return (double) 166.109;
  else
    return (double) 166.115;
}
if(! strcmp(mod, "2-amino-3-oxo-butanoic_acid (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 99.032;
  else
     return (double) 99.088;
}
if(! strcmp(mod, "AB_old_ICATd0 (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 545.234;
  else
     return (double) 545.716;
}
if(! strcmp(mod, "AB_old_ICATd8 (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 553.284;
  else
     return (double) 553.765;
}
if(! strcmp(mod, "Acetyl (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 170.106;
  else
     return (double) 170.209;
}
if(! strcmp(mod, "Acetyl (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 43.0184;
  else
     return (double) 43.0446;
}
if(! strcmp(mod, "Acetyl_heavy (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 173.124;
  else
     return (double) 173.227;
}
if(! strcmp(mod, "Acetyl_heavy (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 46.0372;
  else
     return (double) 46.0631;
}
if(! strcmp(mod, "Acetyl_light (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 170.106;
  else
     return (double) 170.209;
}
if(! strcmp(mod, "Acetyl_light (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 43.0184;
  else
     return (double) 43.0446;
}
if(! strcmp(mod, "Acrylamide_heavy (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 177.065;
  else
     return (double) 177.239;
}
if(! strcmp(mod, "Acrylamide D0 (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 174.04576;
  else
     return (double) 174.1765;
}
if(! strcmp(mod, "Acrylamide D3 (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 177.05582;
  else
     return (double) 173.1865;
}
if(! strcmp(mod, "Amide (C-term)") && strchr("c", label) != NULL) {
  if(monoisotopic)
     return (double) 16.0187;
  else
     return (double) 16.0225;
}
if(! strcmp(mod, "aminotyrosine (Y)") && strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 178.074;
  else
     return (double) 178.188;
}
if(! strcmp(mod, "Arg_heavy (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 166.10111;
  else
     return (double) 166.1875;
}
if(! strcmp(mod, "Argbiotinhydrazide (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 355.168;
  else
     return (double) 355.456;
}
if(! strcmp(mod, "Argglutamicsealde (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 113.048;
  else
     return (double) 113.115;
}
if(! strcmp(mod, "b-methylthiol (D)") && strchr("D", label) != NULL) {
  if(monoisotopic)
     return (double) 161.015;
  else
     return (double) 161.179;
}
if(! strcmp(mod, "Biotin (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 354.173;
  else
     return (double) 354.468;
}
if(! strcmp(mod, "Biotin (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 227.085;
  else
     return (double) 227.303;
}
if(! strcmp(mod, "Biot_LC (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 467.257;
  else
     return (double) 467.625;
}
if(! strcmp(mod, "Biot_LC (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 340.169;
  else
     return (double) 340.461;
}
if(! strcmp(mod, "Carbamidomethyl (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 160.031;
  else
     return (double) 160.194;
}
if(! strcmp(mod, "Carbamidomethyl (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 185.116;
  else
     return (double) 185.224;
}
if(! strcmp(mod, "Carbamidomethyl (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 58.0293;
  else
     return (double) 58.0592;
}
if(! strcmp(mod, "Carbamidomethyl (H)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 194.08;
  else
     return (double) 194.191;
}
if(! strcmp(mod, "Carbamidomethyl (D)") && strchr("D", label) != NULL) {
  if(monoisotopic)
     return (double) 172.048;
  else
     return (double) 172.139;
}
if(! strcmp(mod, "Carbamidomethyl (E)") && strchr("E", label) != NULL) {
  if(monoisotopic)
     return (double) 186.064;
  else
     return (double) 186.165;
}
if(! strcmp(mod, "Carbamyl (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 171.101;
  else
     return (double) 171.197;
}
if(! strcmp(mod, "Carbamyl (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 44.0136;
  else
     return (double) 44.0326;
}
if(! strcmp(mod, "Carbamyl (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 199.107;
  else
     return (double) 199.21;
}
if(! strcmp(mod, "Carbamyl (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 146.015;
  else
     return (double) 146.168;
}
if(! strcmp(mod, "Carboxymethyl (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 161.015;
  else
     return (double) 161.179;
}
if(! strcmp(mod, "Citrullination (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 157.085;
  else
     return (double) 157.171;
}
if(! strcmp(mod, "Cysteic_acid (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 150.994;
  else
     return (double) 151.141;
}
if(! strcmp(mod, "Deamidation (NQ)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 115.027;
  else
     return (double) 115.087;
}
if(! strcmp(mod, "Deamidation (NQ)") && strchr("Q", label) != NULL) {
  if(monoisotopic)
     return (double) 129.043;
  else
     return (double) 129.114;
}
if(! strcmp(mod, "Deamidation_O18 (NQ)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 117.031;
  else
     return (double) 117.087;
}
if(! strcmp(mod, "Deamidation_O18 (NQ)") && strchr("Q", label) != NULL) {
  if(monoisotopic)
     return (double) 131.047;
  else
     return (double) 131.114;
}
if(! strcmp(mod, "Deamidation_O18 (N)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 117.031;
  else
     return (double) 117.087;
}
if(! strcmp(mod, "Deamidation_O18 (Q)") && strchr("Q", label) != NULL) {
  if(monoisotopic)
     return (double) 131.047;
  else
     return (double) 131.114;
}
if(! strcmp(mod, "di-Methylation (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 156.126;
  else
     return (double) 156.225;
}
if(! strcmp(mod, "di-Methylation (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 184.132;
  else
     return (double) 184.239;
}
if(! strcmp(mod, "dihydroxy-phe (F)") && strchr("F", label) != NULL) {
  if(monoisotopic)
     return (double) 179.058;
  else
     return (double) 179.173;
}
if(! strcmp(mod, "DSS (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 138.08;
  else
     return (double) 138.24;
}
if(! strcmp(mod, "DSS_OH (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 156.1;
  else
     return (double) 156.26;
}
if(! strcmp(mod, "EDT-i-biotin (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 577.206;
  else
     return (double) 577.781;
}
if(! strcmp(mod, "EDT-i-biotin (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 591.222;
  else
     return (double) 591.807;
}
if(! strcmp(mod, "EDT-m-biotin (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 688.238;
  else
     return (double) 688.879;
}
if(! strcmp(mod, "EDT-m-biotin (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 702.254;
  else
     return (double) 702.906;
}
if(! strcmp(mod, "ESP-Tag_heavy (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 476.335;
  else
     return (double) 476.702;
}
if(! strcmp(mod, "ESP-Tag_heavy (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 349.248;
  else
     return (double) 349.538;
}
if(! strcmp(mod, "ESP-Tag_light (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 466.273;
  else
     return (double) 466.641;
}
if(! strcmp(mod, "ESP-Tag_light (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 339.185;
  else
     return (double) 339.476;
}
if(! strcmp(mod, "FAD (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 886.151;
  else
     return (double) 886.677;
}
if(! strcmp(mod, "FAD (H)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 920.2;
  else
     return (double) 920.673;
}
if(! strcmp(mod, "Farnesylation (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 307.197;
  else
     return (double) 307.494;
}
if(! strcmp(mod, "formylkynurenin (W)") && strchr("W", label) != NULL) {
  if(monoisotopic)
     return (double) 218.069;
  else
     return (double) 218.209;
}
if(! strcmp(mod, "Gamma-carboxyl (D)") && strchr("D", label) != NULL) {
  if(monoisotopic)
     return (double) 159.017;
  else
     return (double) 159.097;
}
if(! strcmp(mod, "Gamma-carboxyl (E)") && strchr("E", label) != NULL) {
  if(monoisotopic)
     return (double) 173.032;
  else
     return (double) 173.124;
}
if(! strcmp(mod, "Geranyl-geranyl (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 375.26;
  else
     return (double) 375.611;
}
if(! strcmp(mod, "glucuronyl (N-term G)") && strchr("nG", label) != NULL) {
  if(monoisotopic)
     return (double) 233.053548;
  else
     return (double) 233.176;
}
if(! strcmp(mod, "Glutathione (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 408.077;
  else
     return (double) 408.451;
}
if(! strcmp(mod, "glyc_asn_asp (N)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 116.035;
  else
     return (double) 116.095;
}
if(! strcmp(mod, "Guanidination (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 170.117;
  else
     return (double) 170.212;
}
if(! strcmp(mod, "Gygi_ICATd0 (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 589.26;
  else
     return (double) 589.768;
}
if(! strcmp(mod, "Gygi_ICATd8 (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 597.311;
  else
     return (double) 597.818;
}
if(! strcmp(mod, "Hex (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 290.148;
  else
     return (double) 290.313;
}
if(! strcmp(mod, "Hex (N)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 276.096;
  else
     return (double) 276.243;
}
if(! strcmp(mod, "Hex (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 163.061;
  else
     return (double) 163.148;
}
if(! strcmp(mod, "Hex (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 263.1;
  else
     return (double) 263.245;
}
if(! strcmp(mod, "Hex (W)") && strchr("W", label) != NULL) {
  if(monoisotopic)
     return (double) 348.132;
  else
     return (double) 348.35;
}
if(! strcmp(mod, "HexNAc (N)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 317.122;
  else
     return (double) 317.295;
}
if(! strcmp(mod, "HexNAc (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 290.111;
  else
     return (double) 290.27;
}
if(! strcmp(mod, "HexNAc (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 304.127;
  else
     return (double) 304.296;
}
if(! strcmp(mod, "His2Asn (H)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 114.043;
  else
     return (double) 114.103;
}
if(! strcmp(mod, "His2Asp (H)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 115.027;
  else
     return (double) 115.087;
}
if(! strcmp(mod, "HNE (CHK)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 259.124;
  else
     return (double) 259.365;
}
if(! strcmp(mod, "HNE (CHK)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 293.174;
  else
     return (double) 293.361;
}
if(! strcmp(mod, "HNE (CHK)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 284.21;
  else
     return (double) 284.394;
}
if(! strcmp(mod, "HSe (C-term M)") && strchr("cM", label) != NULL) {
  if(monoisotopic)
     return (double) 101.047684;
  else
     return (double) 101.1004;
}
if(! strcmp(mod, "Hse_lact (C-term M)") && strchr("cM", label) != NULL) {
  if(monoisotopic)
     return (double) 83.037119;
  else
     return (double) 83.189229;
}
if(! strcmp(mod, "hydroxykynurenin (W)") && strchr("W", label) != NULL) {
  if(monoisotopic)
     return (double) 206.069;
  else
     return (double) 206.198;
}
if(! strcmp(mod, "Hydroxylation (D)") && strchr("D", label) != NULL) {
  if(monoisotopic)
     return (double) 131.022;
  else
     return (double) 131.087;
}
if(! strcmp(mod, "Hydroxylation (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 144.09;
  else
     return (double) 144.172;
}
if(! strcmp(mod, "Hydroxylation (N)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 130.038;
  else
     return (double) 130.102;
}
if(! strcmp(mod, "Hydroxylation (P)") && strchr("P", label) != NULL) {
  if(monoisotopic)
     return (double) 113.048;
  else
     return (double) 113.115;
}
if(! strcmp(mod, "Hydroxylation (F)") && strchr("F", label) != NULL) {
  if(monoisotopic)
     return (double) 163.063;
  else
     return (double) 163.173;
}
if(! strcmp(mod, "Hydroxylation (Y)") && strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 179.058;
  else
     return (double) 179.173;
}
if(! strcmp(mod, "IBTP (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 419.147;
  else
     return (double) 419.519;
}
if(! strcmp(mod, "ICAT_heavy") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 339.166;
  else
     return (double) 339.337;
}
if(! strcmp(mod, "ICAT_light") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 330.136;
  else
     return (double) 330.403;
}
if(! strcmp(mod, "IMID_heavy (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 200.158;
  else
     return (double) 200.274;
}
if(! strcmp(mod, "IMID_light (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 196.132;
  else
     return (double) 196.25;
}
if(! strcmp(mod, "Im_biotin (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 353.189;
  else
     return (double) 353.483;
}
if(! strcmp(mod, "Im_biotin (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 226.101;
  else
     return (double) 226.318;
}
if(! strcmp(mod, "I-traq (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 143.901;
  else
     return (double) 144.318;
}
if(! strcmp(mod, "MMTS (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 148.997;
  else
     return (double) 149.233;
}
if(! strcmp(mod, "I-Traq (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 272.197;
  else
     return (double) 272.359;
}
if(! strcmp(mod, "kynurenin (W)") && strchr("W", label) != NULL) {
  if(monoisotopic)
     return (double) 190.074;
  else
     return (double) 190.199;
}
if(! strcmp(mod, "Lipoyl (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 316.128;
  else
     return (double) 316.483;
}
if(! strcmp(mod, "Lysaminoadipicsealde (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 127.063;
  else
     return (double) 127.141;
}
if(! strcmp(mod, "Lysbiotinhydrazide (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 369.183;
  else
     return (double) 369.482;
}
if(! strcmp(mod, "Me-ester (C-term)") && strchr("c", label) != NULL) {
  if(monoisotopic)
     return (double) 31.0184;
  else
     return (double) 31.0339;
}
if(! strcmp(mod, "Me-ester (DE)") && strchr("D", label) != NULL) {
  if(monoisotopic)
     return (double) 129.043;
  else
     return (double) 129.114;
}
if(! strcmp(mod, "Me-ester (DE)") && strchr("E", label) != NULL) {
  if(monoisotopic)
     return (double) 143.058;
  else
     return (double) 143.141;
}
if(! strcmp(mod, "Me-ester (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 101.048;
  else
     return (double) 101.104;
}
if(! strcmp(mod, "Me-ester (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 115.063;
  else
     return (double) 115.131;
}
if(! strcmp(mod, "Methyl (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 117.025;
  else
     return (double) 117.17;
}
if(! strcmp(mod, "Methyl (H)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 151.075;
  else
     return (double) 151.166;
}
if(! strcmp(mod, "Methyl (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 142.111;
  else
     return (double) 142.199;
}
if(! strcmp(mod, "Methyl (N)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 128.059;
  else
     return (double) 128.129;
}
if(! strcmp(mod, "Methyl (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 15.0235;
  else
     return (double) 15.0345;
}
if(! strcmp(mod, "Methyl (Q)") && strchr("Q", label) != NULL) {
  if(monoisotopic)
     return (double) 142.074;
  else
     return (double) 142.156;
}
if(! strcmp(mod, "Methyl (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 170.117;
  else
     return (double) 170.212;
}
if(! strcmp(mod, "Myristoylation (N-term G)") && strchr("nG", label) != NULL) {
  if(monoisotopic)
     return (double) 267.219826;
  else
     return (double) 267.4075;
}
if(! strcmp(mod, "Myristoylation (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 338.293;
  else
     return (double) 338.528;
}
if(! strcmp(mod, "N-Acetyl (Protein)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 43.0184;
  else
     return (double) 43.0446;
}
if(! strcmp(mod, "N-Formyl (Protein)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 29.0027;
  else
     return (double) 29.018;
}
if(! strcmp(mod, "NEM (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 228.057;
  else
     return (double) 228.268;
}
if(! strcmp(mod, "NIPCAM (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 202.078;
  else
     return (double) 202.274;
}
if(! strcmp(mod, "Nitro (W)") && strchr("W", label) != NULL) {
  if(monoisotopic)
     return (double) 231.064;
  else
     return (double) 231.208;
}
if(! strcmp(mod, "Nitro (Y)") && strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 208.048;
  else
     return (double) 208.171;
}
/*else*/ if(! strcmp(mod, "N->D (N->D)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 115.02694;
  else
     return (double) 115.0886;
}
if(! strcmp(mod, "double_O18 (C-term)") && strchr("c", label) != NULL) {
  if(monoisotopic)
     return (double) 21.011;
  else
     return (double) 21.0078;
}
if(! strcmp(mod, "O18 (C-term)") && strchr("c", label) != NULL) {
  if(monoisotopic)
     return (double) 19.007;
  else
     return (double) 19.0071;
}
if(! strcmp(mod, "OxArgBiotin (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 466.236;
  else
     return (double) 466.598;
}
if(! strcmp(mod, "OxArgBiotinRed (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 468.252;
  else
     return (double) 468.613;
}
if(! strcmp(mod, "OxGln1 (Q)") && strchr("Q", label) != NULL) {
  if(monoisotopic)
     return (double) 70.0055;
  else
     return (double) 70.0468;
}
if(! strcmp(mod, "OxGln2 (Q)") && strchr("Q", label) != NULL) {
  if(monoisotopic)
     return (double) 71.9847;
  else
     return (double) 72.0196;
}
if(! strcmp(mod, "Oxidation (M)") && strchr("M", label) != NULL) {
  if(monoisotopic)
     return (double) 147.035;
  else
     return (double) 147.195;
}
if(! strcmp(mod, "Oxidation (HW)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 153.054;
  else
     return (double) 153.139;
}
if(! strcmp(mod, "Oxidation (HW)") && strchr("W", label) != NULL) {
  if(monoisotopic)
     return (double) 202.074;
  else
     return (double) 202.209;
}
if(! strcmp(mod, "OxLysBiotin (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 480.252;
  else
     return (double) 480.624;
}
if(! strcmp(mod, "OxLysBiotinRed (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 482.268;
  else
     return (double) 482.64;
}
if(! strcmp(mod, "OxProBiotin (P)") && strchr("P", label) != NULL) {
  if(monoisotopic)
     return (double) 466.236;
  else
     return (double) 466.598;
}
if(! strcmp(mod, "OxProBiotinRed (P)") && strchr("P", label) != NULL) {
  if(monoisotopic)
     return (double) 468.252;
  else
     return (double) 468.613;
}
if(! strcmp(mod, "p-pantetheine (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 426.11;
  else
     return (double) 426.402;
}
if(! strcmp(mod, "Palmitoylation (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 341.239;
  else
     return (double) 341.552;
}
if(! strcmp(mod, "Palmitoylation (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 366.325;
  else
     return (double) 366.581;
}
if(! strcmp(mod, "Palmitoylation (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 325.262;
  else
     return (double) 325.486;
}
if(! strcmp(mod, "Palmitoylation (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 339.277;
  else
     return (double) 339.513;
}
if(! strcmp(mod, "PEO-Biotin (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 517.203;
  else
     return (double) 517.662;
}
if(! strcmp(mod, "Phospho (ST)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 166.998;
  else
     return (double) 167.057;
}
if(! strcmp(mod, "Phospho (ST)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 181.014;
  else
     return (double) 181.084;
}
if(! strcmp(mod, "Phospho (Y)") && strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 243.03;
  else
     return (double) 243.153;
}
if(! strcmp(mod, "Phospho (STY)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 166.998;
  else
     return (double) 167.057;
}
if(! strcmp(mod, "Phospho (STY)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 181.014;
  else
     return (double) 181.084;
}
if(! strcmp(mod, "Phospho (STY)") && strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 243.03;
  else
     return (double) 243.153;
}
if(! strcmp(mod, "Phospho+PL (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 69.0215;
  else
     return (double) 69.062;
}
if(! strcmp(mod, "Phospho+PL (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 83.0371;
  else
     return (double) 83.0886;
}
if(! strcmp(mod, "Phospho+PL (Y)") && strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 145.053;
  else
     return (double) 145.158;
}
if(! strcmp(mod, "Phospho-NL (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 182.976;
  else
     return (double) 183.123;
}
if(! strcmp(mod, "Phospho-NL (D)") && strchr("D", label) != NULL) {
  if(monoisotopic)
     return (double) 194.993;
  else
     return (double) 195.067;
}
if(! strcmp(mod, "Phospho-NL (H)") && strchr("H", label) != NULL) {
  if(monoisotopic)
     return (double) 217.025;
  else
     return (double) 217.119;
}
if(! strcmp(mod, "Phospho-NL (S)") && strchr("S", label) != NULL) {
  if(monoisotopic)
     return (double) 166.998;
  else
     return (double) 167.057;
}
if(! strcmp(mod, "Phospho-NL (T)") && strchr("T", label) != NULL) {
  if(monoisotopic)
     return (double) 181.014;
  else
     return (double) 181.084;
}
if(! strcmp(mod, "Probiotinhydrazide (P)") && strchr("P", label) != NULL) {
  if(monoisotopic)
     return (double) 355.168;
  else
     return (double) 355.456;
}
if(! strcmp(mod, "Proglutamicsealde (P)") && strchr("P", label) != NULL) {
  if(monoisotopic)
     return (double) 127.063;
  else
     return (double) 127.141;
}
if(! strcmp(mod, "Propionamide (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 174.046;
  else
     return (double) 174.221;
}
if(! strcmp(mod, "Propionyl_heavy (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 187.131;
  else
     return (double) 187.214;
}
if(! strcmp(mod, "Propionyl_heavy (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 60.0441;
  else
     return (double) 60.0491;
}
if(! strcmp(mod, "Propionyl_light (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 184.121;
  else
     return (double) 184.236;
}
if(! strcmp(mod, "Propionyl_light (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 57.034;
  else
     return (double) 57.0712;
}
if(! strcmp(mod, "Pyridoxal-phos (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 357.109;
  else
     return (double) 357.299;
}
if(! strcmp(mod, "Pyridyl (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 247.132;
  else
     return (double) 247.293;
}
if(! strcmp(mod, "Pyridyl (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 120.045;
  else
     return (double) 120.129;
}
if(! strcmp(mod, "Pyro-cmC (N-term camC)") && strchr("nC", label) != NULL) {
  if(monoisotopic)
     return (double) 85.982641;
  else
     return (double) 86.1083;
}
if(! strcmp(mod, "Pyro-glu (N-term E)") && strchr("nE", label) != NULL) {
  if(monoisotopic)
     return (double) 111.032025;
  else
     return (double) 111.1002;
}
if(! strcmp(mod, "Pyro-glu (N-term Q)") && strchr("nQ", label) != NULL) { //Modification on Q that has to be terminal
  if(monoisotopic)
    return (double) 111.032031; 
  else
    return (double) 111.1002;
}
if(! strcmp(mod, "Pyroglutamic (P)") && strchr("P", label) != NULL) {
  if(monoisotopic)
     return (double) 111.032;
  else
     return (double) 111.099;
}
if(! strcmp(mod, "Pyrrolidinone (P)") && strchr("P", label) != NULL) {
  if(monoisotopic)
     return (double) 67.0422;
  else
     return (double) 67.0892;
}
if(! strcmp(mod, "Quat_0 (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 255.195;
  else
     return (double) 255.357;
}
if(! strcmp(mod, "Quat_0 (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 128.108;
  else
     return (double) 128.192;
}
if(! strcmp(mod, "Quat_3 (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 258.214;
  else
     return (double) 258.375;
}
if(! strcmp(mod, "Quat_3 (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 131.126;
  else
     return (double) 131.211;
}
if(! strcmp(mod, "Quat_6 (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 261.232;
  else
     return (double) 261.393;
}
if(! strcmp(mod, "Quat_6 (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 134.145;
  else
     return (double) 134.229;
}
if(! strcmp(mod, "Quat_9 (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 264.251;
  else
     return (double) 264.412;
}
if(! strcmp(mod, "Quat_9 (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 137.164;
  else
     return (double) 137.248;
}
if(! strcmp(mod, "S-pyridylethyl (C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 208.067;
  else
     return (double) 208.28;
}
if(! strcmp(mod, "SMA (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 255.158;
  else
     return (double) 255.314;
}
if(! strcmp(mod, "SMA (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 128.071;
  else
     return (double) 128.149;
}
if(! strcmp(mod, "Sodiated (C-term)") && strchr("c", label) != NULL) {
  if(monoisotopic)
     return (double) 38.9847;
  else
     return (double) 38.9891;
}
if(! strcmp(mod, "Sodiated (DE)") && strchr("D", label) != NULL) {
  if(monoisotopic)
     return (double) 137.009;
  else
     return (double) 137.069;
}
if(! strcmp(mod, "Sodiated (DE)") && strchr("E", label) != NULL) {
  if(monoisotopic)
     return (double) 151.025;
  else
     return (double) 151.096;
}
if(! strcmp(mod, "Suc_anh+4C13 (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 232.124;
  else
     return (double) 232.216;
}
if(! strcmp(mod, "Suc_anh+4C13 (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 105.037;
  else
     return (double) 105.051;
}
if(! strcmp(mod, "Suc_anh+4H2 (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 232.136;
  else
     return (double) 232.27;
}
if(! strcmp(mod, "Suc_anh+4H2 (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 105.049;
  else
     return (double) 105.105;
}
if(! strcmp(mod, "Suc_anh_light (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 228.111;
  else
     return (double) 228.245;
}
if(! strcmp(mod, "Suc_anh_light (N-term)") && strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 101.024;
  else
     return (double) 101.081;
}
if(! strcmp(mod, "Sulfation (Y)") && strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 243.02;
  else
     return (double) 243.236;
}
if(! strcmp(mod, "Sulphone (M)") && strchr("M", label) != NULL) {
  if(monoisotopic)
     return (double) 163.03;
  else
     return (double) 163.195;
}
if(! strcmp(mod, "tri-Methylation (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 170.142;
  else
     return (double) 170.252;
}
if(! strcmp(mod, "Tripalmitate (N-term C)") && strchr("nC", label) != NULL) {
  if(monoisotopic)
     return (double) 891.734967;
  else
     return (double) 892.4437;
}
if(! strcmp(mod, "Glycosylate (N)") && strchr("N", label) != NULL) {
  if(monoisotopic)
     return (double) 115.02693;
  else
     return (double) 115.08;
}
if(! strcmp(mod, "SILAC (K)") && strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 134.09426;
  else
     return (double) 134.1741;
}
if(! strcmp(mod, "SILAC (R)") && strchr("R", label) != NULL) {
  if(monoisotopic)
     return (double) 166.10111;
  else
     return (double) 166.1875;
}
if((! strcmp(mod, "(N-term)_iTRAQ") || ! strcmp(mod, "iTRAQ (N-term)")) 
	&& strchr("n", label) != NULL) {
  if(monoisotopic)
     return (double) 145.109888;
  else
     return (double) 145.162225;
}
if((! strcmp(mod, "Lysine(K)_iTRAQ") || ! strcmp(mod, "iTRAQ (K)")) 
	&& strchr("K", label) != NULL) {
  if(monoisotopic)
     return (double) 272.197026;
  else
     return (double) 272.3267;
}
if((! strcmp(mod, "Tyrosine(Y)_iTRAQ") || ! strcmp(mod, "iTRAQ (Y)")) 
	&& strchr("Y", label) != NULL) {
  if(monoisotopic)
     return (double) 307.165392;
  else
     return (double) 307.3277;
}
if(! strcmp(mod, "MMTS(C)") && strchr("C", label) != NULL) {
  if(monoisotopic)
     return (double) 148.9969;
  else
     return (double) 149.2345;
}
return (double) 0.0;
}

const char* ResidueMass::getStdModResidues(const char* mod) {
  Boolean n_term, c_term;
  return getStdModResidues(mod, n_term, c_term);
}

const char* ResidueMass::getStdModResidues(const char* mod, Boolean& n_term_aa_mod, Boolean& c_term_aa_mod) {
   load_table(); // in case we didn't already
   std::pair<rmap::const_iterator, rmap::const_iterator> p =
      residues.equal_range(std::string(mod)); //  set up iterator
   for (rmap::const_iterator i = p.first; i != p.second; ++i) {
      const cResidue &r = (*i).second;
      // note that this typically leaves uninit ?_term_aa_mod values, but
      // this is per the default implementation in tpp
      if (r.m_masses[0].m_nterm) {
         n_term_aa_mod = true;
      } else if (r.m_masses[0].m_cterm) {
         c_term_aa_mod = true;
      }
      return r.m_residue.c_str();
   }
  if(! strcmp(mod, "+N-formyl-met (Protein)")) {
    return "n";
  } if (! strcmp(mod, "13C6-15N2 (K)")) {
    return "K";
  } if (! strcmp(mod, "13C6-15N4 (R)")) {
    return "R";
  } if(! strcmp(mod, "2-amino-3-oxo-butanoic_acid (T)")) {
    return "T";
  } if(! strcmp(mod, "AB_old_ICATd0 (C)")) {
    return "C";
  } if(! strcmp(mod, "AB_old_ICATd8 (C)")) {
    return "C";
  } if(! strcmp(mod, "Acetyl (K)")) {
    return "K";
  } if(! strcmp(mod, "Acetyl (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Acetyl_heavy (K)")) {
    return "K";
  } if(! strcmp(mod, "Acetyl_heavy (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Acetyl_light (K)")) {
    return "K";
  } if(! strcmp(mod, "Acetyl_light (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Acrylamide_heavy (C)")) {
    return "C";
  } if(! strcmp(mod, "Acrylamide D0 (C)")) {
    return "C";
  } if(! strcmp(mod, "Acrylamide D3 (C)")) {
    return "C";
  } if(! strcmp(mod, "Amide (C-term)")) {
    return "c";
  } if(! strcmp(mod, "aminotyrosine (Y)")) {
    return "Y";
  } if(! strcmp(mod, "Argbiotinhydrazide (R)")) {
    return "R";
  } if(! strcmp(mod, "Argglutamicsealde (R)")) {
    return "R";
  } if(! strcmp(mod, "Arg_heavy (R)")) {
    return "R";
  } if(! strcmp(mod, "b-methylthiol (D)")) {
    return "D";
  } if(! strcmp(mod, "Biotin (K)")) {
    return "K";
  } if(! strcmp(mod, "Biotin (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Biot_LC (K)")) {
    return "K";
  } if(! strcmp(mod, "Biot_LC (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Carbamidomethyl (C)")) {
    return "C";
  } if(! strcmp(mod, "Carbamidomethyl (K)")) {
    return "K";
  } if(! strcmp(mod, "Carbamidomethyl (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Carbamidomethyl (H)")) {
    return "H";
  } if(! strcmp(mod, "Carbamidomethyl (D)")) {
    return "D";
  } if(! strcmp(mod, "Carbamidomethyl (E)")) {
    return "E";
  } if(! strcmp(mod, "Carbamyl (K)")) {
    return "K";
  } if(! strcmp(mod, "Carbamyl (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Carbamyl (R)")) {
    return "R";
  } if(! strcmp(mod, "Carbamyl (C)")) {
    return "C";
  } if(! strcmp(mod, "Carboxymethyl (C)")) {
    return "C";
  } if(! strcmp(mod, "Citrullination (R)")) {
    return "R";
  } if(! strcmp(mod, "Cysteic_acid (C)")) {
    return "C";
  } if(! strcmp(mod, "Deamidation (NQ)")) {
    return "NQ";
  } if(! strcmp(mod, "Deamidation_O18 (NQ)")) {
    return "NQ";
  } if(! strcmp(mod, "Deamidation_O18 (Q)")) {
    return "Q";
  } if(! strcmp(mod, "di-Methylation (K)")) {
    return "K";
  } if(! strcmp(mod, "di-Methylation (R)")) {
    return "R";
  } if(! strcmp(mod, "dihydroxy-phe (F)")) {
    return "F";
  } if(! strcmp(mod, "DSS (K)")) {
    return "K";
  } if(! strcmp(mod, "DSS_OH (K)")) {
    return "K";
  } if(! strcmp(mod, "EDT-i-biotin (S)")) {
    return "S";
  } if(! strcmp(mod, "EDT-i-biotin (T)")) {
    return "T";
  } if(! strcmp(mod, "EDT-m-biotin (S)")) {
    return "S";
  } if(! strcmp(mod, "EDT-m-biotin (T)")) {
    return "T";
  } if(! strcmp(mod, "ESP-Tag_heavy (K)")) {
    return "K";
  } if(! strcmp(mod, "ESP-Tag_heavy (N-term)")) {
    return "n";
  } if(! strcmp(mod, "ESP-Tag_light (K)")) {
    return "K";
  } if(! strcmp(mod, "ESP-Tag_light (N-term)")) {
    return "n";
  } if(! strcmp(mod, "FAD (C)")) {
    return "C";
  } if(! strcmp(mod, "FAD (H)")) {
    return "H";
  } if(! strcmp(mod, "Farnesylation (C)")) {
    return "C";
  } if(! strcmp(mod, "formylkynurenin (W)")) {
    return "W";
  } if(! strcmp(mod, "Gamma-carboxyl (D)")) {
    return "D";
  } if(! strcmp(mod, "Gamma-carboxyl (E)")) {
    return "E";
  } if(! strcmp(mod, "Geranyl-geranyl (C)")) {
    return "C";
  } if(! strcmp(mod, "glucuronyl (N-term G)")) {
    n_term_aa_mod = true;
    return "G";
  } if(! strcmp(mod, "Glutathione (C)")) {
    return "C";
  } if(! strcmp(mod, "glyc_asn_asp (N)")) {
    return "N";
  } if(! strcmp(mod, "Guanidination (K)")) {
    return "K";
  } if(! strcmp(mod, "Gygi_ICATd0 (C)")) {
    return "C";
  } if(! strcmp(mod, "Gygi_ICATd8 (C)")) {
    return "C";
  } if(! strcmp(mod, "Hex (K)")) {
    return "K";
  } if(! strcmp(mod, "Hex (N)")) {
    return "N";
  } if(! strcmp(mod, "Hex (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Hex (T)")) {
    return "T";
  } if(! strcmp(mod, "Hex (W)")) {
    return "W";
  } if(! strcmp(mod, "HexNAc (N)")) {
    return "N";
  } if(! strcmp(mod, "HexNAc (S)")) {
    return "S";
  } if(! strcmp(mod, "HexNAc (T)")) {
    return "T";
  } if(! strcmp(mod, "His2Asn (H)")) {
    return "H";
  } if(! strcmp(mod, "His2Asp (H)")) {
    return "H";
  } if(! strcmp(mod, "HNE (CHK)")) {
    return "CHK";
  } if(! strcmp(mod, "HSe (C-term M)")) {
    c_term_aa_mod = true;
    return "M";
  } if(! strcmp(mod, "Hse_lact (C-term M)")) {
    c_term_aa_mod = true;
    return "M";
  } if(! strcmp(mod, "hydroxykynurenin (W)")) {
    return "W";
  } if(! strcmp(mod, "Hydroxylation (D)")) {
    return "D";
  } if(! strcmp(mod, "Hydroxylation (K)")) {
    return "K";
  } if(! strcmp(mod, "Hydroxylation (N)")) {
    return "N";
  } if(! strcmp(mod, "Hydroxylation (P)")) {
    return "P";
  } if(! strcmp(mod, "Hydroxylation (F)")) {
    return "F";
  } if(! strcmp(mod, "Hydroxylation (Y)")) {
    return "Y";
  } if(! strcmp(mod, "IBTP (C)")) {
    return "C";
  } if(! strcmp(mod, "ICAT_heavy")) {
    return "C";
  } if(! strcmp(mod, "ICAT_light")) {
    return "C";
  } if(! strcmp(mod, "IMID_heavy (K)")) {
    return "K";
  } if(! strcmp(mod, "IMID_light (K)")) {
    return "K";
  } if(! strcmp(mod, "Im_biotin (K)")) {
    return "K";
  } if(! strcmp(mod, "Im_biotin (N-term)")) {
    return "n";
  } if(! strcmp(mod, "I-traq (N-term)")) {
    return "n";
  } if(! strcmp(mod, "MMTS (C)")) {
    return "C";
  } if(! strcmp(mod, "I-Traq (K)")) {
    return "K";
  } if(! strcmp(mod, "kynurenin (W)")) {
    return "W";
  } if(! strcmp(mod, "Lipoyl (K)")) {
    return "K";
  } if(! strcmp(mod, "Lysaminoadipicsealde (K)")) {
    return "K";
  } if(! strcmp(mod, "Lysbiotinhydrazide (K)")) {
    return "K";
  } if(! strcmp(mod, "Me-ester (C-term)")) {
    return "c";
  } if(! strcmp(mod, "Me-ester (DE)")) {
    return "DE";
  } if(! strcmp(mod, "Me-ester (S)")) {
    return "S";
  } if(! strcmp(mod, "Me-ester (T)")) {
    return "T";
  } if(! strcmp(mod, "Methyl (C)")) {
    return "C";
  } if(! strcmp(mod, "Methyl (H)")) {
    return "H";
  } if(! strcmp(mod, "Methyl (K)")) {
    return "K";
  } if(! strcmp(mod, "Methyl (N)")) {
    return "N";
  } if(! strcmp(mod, "Methyl (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Methyl (Q)")) {
    return "Q";
  } if(! strcmp(mod, "Methyl (R)")) {
    return "R";
  } if(! strcmp(mod, "Myristoylation (N-term G)")) {
    n_term_aa_mod = true;
    return "G";
  } if(! strcmp(mod, "Myristoylation (K)")) {
    return "K";
  } if(! strcmp(mod, "N-Acetyl (Protein)")) {
    return "n";
  } if(! strcmp(mod, "N-Formyl (Protein)")) {
    return "n";
  } if(! strcmp(mod, "NEM (C)")) {
    return "C";
  } if(! strcmp(mod, "NIPCAM (C)")) {
    return "C";
  } if(! strcmp(mod, "Nitro (W)")) {
    return "W";
  } if(! strcmp(mod, "Nitro (Y)")) {
    return "Y";
  } if(! strcmp(mod, "N->D (N->D)")) {
    return "N";
  } if(! strcmp(mod, "O18 (C-term)")) {
    return "c";
  } if(! strcmp(mod, "double_O18 (C-term)")) {
    return "c";
  } if(! strcmp(mod, "OxArgBiotin (R)")) {
    return "R";
  } if(! strcmp(mod, "OxArgBiotinRed (R)")) {
    return "R";
  } if(! strcmp(mod, "OxGln1 (Q)")) {
    return "Q";
  } if(! strcmp(mod, "OxGln2 (Q)")) {
    return "Q";
  } if(! strcmp(mod, "Oxidation (M)")) {
    return "M";
  } if(! strcmp(mod, "Oxidation (HW)")) {
    return "HW";
//  } if(! strcmp(mod, "Oxidation (HW)")) {
//return "W";
  } if(! strcmp(mod, "OxLysBiotin (K)")) {
    return "K";
  } if(! strcmp(mod, "OxLysBiotinRed (K)")) {
    return "K";
  } if(! strcmp(mod, "OxProBiotin (P)")) {
    return "P";
  } if(! strcmp(mod, "OxProBiotinRed (P)")) {
    return "P";
  } if(! strcmp(mod, "p-pantetheine (S)")) {
    return "S";
  } if(! strcmp(mod, "Palmitoylation (C)")) {
    return "C";
  } if(! strcmp(mod, "Palmitoylation (K)")) {
    return "K";
  } if(! strcmp(mod, "Palmitoylation (S)")) {
    return "S";
  } if(! strcmp(mod, "Palmitoylation (T)")) {
    return "T";
  } if(! strcmp(mod, "PEO-Biotin (C)")) {
    return "C";
  } if(! strcmp(mod, "Phospho (ST)")) {
    return "ST";
  } if(! strcmp(mod, "Phospho (Y)")) {
    return "Y";
  } if(! strcmp(mod, "Phospho (STY)")) {
    return "STY";
  } if(! strcmp(mod, "Phospho+PL (S)")) {
    return "S";
  } if(! strcmp(mod, "Phospho+PL (T)")) {
    return "T";
  } if(! strcmp(mod, "Phospho+PL (Y)")) {
    return "Y";
  } if(! strcmp(mod, "Phospho-NL (C)")) {
    return "C";
  } if(! strcmp(mod, "Phospho-NL (D)")) {
    return "D";
  } if(! strcmp(mod, "Phospho-NL (H)")) {
    return "H";
  } if(! strcmp(mod, "Phospho-NL (S)")) {
    return "S";
  } if(! strcmp(mod, "Phospho-NL (T)")) {
    return "T";
  } if(! strcmp(mod, "Probiotinhydrazide (P)")) {
    return "P";
  } if(! strcmp(mod, "Proglutamicsealde (P)")) {
    return "P";
  } if(! strcmp(mod, "Propionamide (C)")) {
    return "C";
  } if(! strcmp(mod, "Propionyl_heavy (K)")) {
    return "K";
  } if(! strcmp(mod, "Propionyl_heavy (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Propionyl_light (K)")) {
    return "K";
  } if(! strcmp(mod, "Propionyl_light (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Pyridoxal-phos (K)")) {
    return "K";
  } if(! strcmp(mod, "Pyridyl (K)")) {
    return "K";
  } if(! strcmp(mod, "Pyridyl (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Pyro-cmC (N-term camC)")) {
    n_term_aa_mod = true;
    return "C";
  } if(! strcmp(mod, "Pyro-glu (N-term E)")) {
    n_term_aa_mod = true;
    return "E";
  } if(! strcmp(mod, "Pyro-glu (N-term Q)")) {
    n_term_aa_mod = true;
    return "Q";
  } if(! strcmp(mod, "Pyroglutamic (P)")) {
    return "P";
  } if(! strcmp(mod, "Pyrrolidinone (P)")) {
    return "P";
  } if(! strcmp(mod, "Quat_0 (K)")) {
    return "K";
  } if(! strcmp(mod, "Quat_0 (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Quat_3 (K)")) {
    return "K";
  } if(! strcmp(mod, "Quat_3 (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Quat_6 (K)")) {
    return "K";
  } if(! strcmp(mod, "Quat_6 (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Quat_9 (K)")) {
    return "K";
  } if(! strcmp(mod, "Quat_9 (N-term)")) {
    return "n";
  } if(! strcmp(mod, "S-pyridylethyl (C)")) {
    return "C";
  } if(! strcmp(mod, "SMA (K)")) {
    return "K";
  } if(! strcmp(mod, "SMA (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Sodiated (C-term)")) {
    return "c";
  } if(! strcmp(mod, "Sodiated (DE)")) {
    return "DE";
//} if(! strcmp(mod, "Sodiated (DE)")) {
//  return "E";
  } if(! strcmp(mod, "Suc_anh+4C13 (K)")) {
    return "K";
  } if(! strcmp(mod, "Suc_anh+4C13 (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Suc_anh+4H2 (K)")) {
    return "K";
  } if(! strcmp(mod, "Suc_anh+4H2 (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Suc_anh_light (K)")) {
    return "K";
  } if(! strcmp(mod, "Suc_anh_light (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Sulfation (Y)")) {
    return "Y";
  } if(! strcmp(mod, "Sulphone (M)")) {
    return "M";
  } if(! strcmp(mod, "tri-Methylation (K)")) {
    return "K";
  } if(! strcmp(mod, "Tripalmitate (N-term C)")) {
    n_term_aa_mod = true;
    return "C";
  } if(! strcmp(mod, "Glycosylate (N)")) {
    return "N";
  } if(! strcmp(mod, "SILAC (K)")) {
    return "K";
  } if(! strcmp(mod, "SILAC (R)")) {
    return "R";
  } if(! strcmp(mod, "(N-term)_iTRAQ") || ! strcmp(mod, "iTRAQ (N-term)")) {
    return "n";
  } if(! strcmp(mod, "Lysine(K)_iTRAQ")|| ! strcmp(mod, "iTRAQ (K)")) {
    return "K";
  } if(! strcmp(mod, "Tyrosine(Y)_iTRAQ")|| ! strcmp(mod, "iTRAQ (Y)")) {
    return "Y";
  } if(! strcmp(mod, "MMTS(C)")) {
    return "C";
  }
  return "";
}

//
// stuff to load user-specified residue masses
// Copyright (C) Insilicos LLC 2006-2010 All Rights Reserved
//

static bool table_load_attempted = false;

// strip leading spaces, trailing newline from fgets
static void tidyrecord(char *str) {
   while (' '==*str) { // eat leading spaces
      memmove(str,str+1,strlen(str));
   }
   for (char *cp=str;*cp;) { // eat any newlines
      if (strchr("\r\n",*cp)) {
         memmove(cp,cp+1,strlen(cp));
      } else {
         cp++;
      }
   }
}

// parse masses from input line of form 
// "<keyword>:<res> <mono> <avg>" or "<keyword>:<mono> <avg>"
static bool readmasses(char *buf,resinfo_t &ma) {
   ma.m_cterm = false;
   ma.m_nterm = false;
   char *space1 = strchr(buf,' ');
   char *space2 = space1?strchr(space1+1,' '):NULL;
   if (space1 && !space2) { // probably ProteinNterm:<val> <val>
      space1 = strchr(buf,':');
      space2 = space1?strchr(buf,' '):NULL;
   }
   if (space1 && space2) {
      ma.m_monomass = atof(++space1);
      ma.m_mass = atof(++space2);
   }
   return (space1 && space2);
}

static void load_table() {
   if (!table_load_attempted) {
      std::string fname = getConfPath();
	  fname += "residue_modifications.cfg";
      FILE *fp=fopen(fname.c_str(),"r");
      bool datahit = false;
      if (fp) {
         char buf[1024];
         char *name = NULL;
         cResidue *residue=NULL;
         while(fgets(buf,sizeof(buf),fp)) {

            tidyrecord(buf);
            if (buf[0] != '#') { // not a comment
               struct stat statbuf;
               if (!datahit && !stat(buf,&statbuf)) { // external file decl?
                  FILE *ext=fopen(buf,"r");
                  if (ext) {
                     fclose(fp);
                     fp = ext;
                     fname = buf;
                     if (fgets(buf,sizeof(buf),fp)) { // read from ext file
                         tidyrecord(buf);
                     }
                     if (buf[0] == '#') { // a comment
                         continue;
                     }
                  }
               }
               if (!strncmp(buf,"Title:",6)) { // start of record
                  if (!datahit) {
                     datahit = true;
                     std::cerr << "loading residue modifications from " << fname.c_str() << std::endl;
                  }
                  name = strdup(buf+6); // grab everything after title
                  residue = new cResidue;
               } else if ('*'==*buf) {   // end of record
                  residues.insert(rmap::value_type(std::string(name),*residue));
                  // watch for case like this:
                  // Title:Label:13C(6)15N(2) (K)
                  // Residues:K 136.109162 136.1150
                  // *
                  // but mascot .dat file has
                  // 13C6-15N2 (K)
                  bool weirdness = false;
                  if (!strncasecmp(name,"Label:",6)) {
                     memmove(name,name+6,strlen(name+5));
                     weirdness = true;
                  }
                  // convert 13C(6)15N(2) (K) to 13C6-15N2 (K)
                  char *cp;
                  while ((cp=strchr(name,'(')) != strrchr(name,'(')) {
                     weirdness = true;
                     memmove(cp,cp+1,strlen(cp));
                  }
                  while ((cp=strchr(name,')')) != strrchr(name,')')) {
                     weirdness = true;
                     if (' '==*(cp+1)) {
                        memmove(cp,cp+1,strlen(cp));
                     } else {
                        *cp = '-';
                     }
                  }
                  if (weirdness) {
                     residues.insert(rmap::value_type(std::string(name),*residue));
                  }
                  delete residue;
                  residue = NULL;
                  free(name);
               } else if (strstr(buf,"Residues:")) {
                  resinfo_t ma;
                  if (readmasses(buf,ma)) {
                     char res[2];
                     res[0] = *(strchr(buf,':')+1);
                     res[1] = 0;
                     residue->m_residue += res;
                     residue->m_masses.push_back(ma);
                  }
               } else if (strstr(buf,"ResiduesNterm:")) {
                  resinfo_t ma;
                  if (readmasses(buf,ma)) {
                     char res[2];
                     res[0] = *(strchr(buf,':')+1);
                     res[1] = 0;
                     residue->m_residue += res;
                     ma.m_nterm = true;
                     residue->m_masses.push_back(ma);
                  }
               } else if (strstr(buf,"ResiduesCterm:")) {
                  resinfo_t ma;
                  if (readmasses(buf,ma)) {
                     char res[2];
                     res[0] = *(strchr(buf,':')+1);
                     res[1] = 0;
                     residue->m_residue += res;
                     ma.m_cterm = true;
                     residue->m_masses.push_back(ma);
                  }
               } else if (strstr(buf,"ProteinNterm:")) {
                  resinfo_t ma;
                  if (readmasses(buf,ma)) {
                     residue->m_residue += "n";
                     ma.m_nterm = true;
                     residue->m_masses.push_back(ma);
                  }
               } else if (strstr(buf,"ProteinCterm:")) {
                  resinfo_t ma;
                  if (readmasses(buf,ma)) {
                     residue->m_residue += "c";
                     ma.m_cterm = true;
                     residue->m_masses.push_back(ma);
                  }
               } else if (strstr(buf,"Nterm:")) {
                  resinfo_t ma;
                  if (readmasses(buf,ma)) {
                     residue->m_residue += "n";
                     ma.m_nterm = true;
                     residue->m_masses.push_back(ma);
                  }
               } else if (strstr(buf,"Cterm:")) {
                  resinfo_t ma;
                  if (readmasses(buf,ma)) {
                     residue->m_residue += "c";
                     ma.m_cterm = true;
                     residue->m_masses.push_back(ma);
                  }
               }
            }
         }
         fclose(fp);
      }
      table_load_attempted = true;
   }
}
