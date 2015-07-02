#include "AccurateMonoMass.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;


static bool map_generated_ = false;
static map<char, MolecularFormula> amino_acid_to_formula_;

#define H 0
#define C 1
#define N 2
#define O 3
#define P 4
#define S 5

double mass_table[] = {
  1.00782503207, //1H
  12.0000000,    //12C
  14.0030740048, //14N
  15.99491461956, //16O
  30.97376163,    //31P
  31.97207100  //32S
  };


MolecularFormula::MolecularFormula() {
  H_ = 0;
  C_ = 0;
  N_ = 0;
  O_ = 0;
  P_ = 0;
  S_ = 0;
}

MolecularFormula::~MolecularFormula() {
  ;
}

MolecularFormula::MolecularFormula(int nH, int nC, int nN, int nO, int nP, int nS) {
  H_ = nH;
  C_ = nC;
  N_ = nN;
  O_ = nO;
  P_ = nP;
  S_ = nS;
}

void MolecularFormula::add(MolecularFormula& formula) {
  H_ += formula.H_;
  C_ += formula.C_;
  N_ += formula.N_;
  O_ += formula.O_;
  P_ += formula.P_;
  S_ += formula.S_;
}

double MolecularFormula::calculateMass() {
  double ans = 0.0;

  ans += mass_table[H] * H_;
  ans += mass_table[C] * C_;
  ans += mass_table[N] * N_;
  ans += mass_table[O] * O_;
  ans += mass_table[P] * P_;
  ans += mass_table[S] * S_;

  return ans;
}


static MolecularFormula H2O(2, 0, 0, 1, 0, 0);

void AccurateMonoMass::generateMap() {

  amino_acid_to_formula_['A'] = MolecularFormula(5,  3, 1, 1, 0, 0);
  amino_acid_to_formula_['C'] = MolecularFormula(5,  3, 1, 1, 0, 1);
  amino_acid_to_formula_['c'] = MolecularFormula(8,  5, 2, 2, 0, 1); //iodoacetmide adds C2H3NO
  amino_acid_to_formula_['D'] = MolecularFormula(5,  4, 1, 3, 0, 0);
  amino_acid_to_formula_['E'] = MolecularFormula(7,  5, 1, 3, 0, 0);
  amino_acid_to_formula_['F'] = MolecularFormula(9,  9, 1, 1, 0, 0);
  amino_acid_to_formula_['G'] = MolecularFormula(3,  2, 1, 1, 0, 0);
  amino_acid_to_formula_['H'] = MolecularFormula(7,  6, 3, 1, 0, 0);
  amino_acid_to_formula_['I'] = MolecularFormula(11, 6, 1, 1, 0, 0);
  amino_acid_to_formula_['K'] = MolecularFormula(12, 6, 2, 1, 0, 0);
  amino_acid_to_formula_['L'] = MolecularFormula(11, 6, 1, 1, 0, 0);
  amino_acid_to_formula_['M'] = MolecularFormula(9,  5, 1, 1, 0, 1);
  //amino_acid_to_formula['m'] = MolecularFormula(9, 5, 1, 2, 0, 1); # oxidized methionine
  amino_acid_to_formula_['N'] = MolecularFormula(6,  4, 2, 2, 0, 0);
  amino_acid_to_formula_['P'] = MolecularFormula(7,  5, 1, 1, 0, 0);
  amino_acid_to_formula_['Q'] = MolecularFormula(8,  5, 2, 2, 0, 0);
  amino_acid_to_formula_['R'] = MolecularFormula(12, 6, 4, 1, 0, 0);
  amino_acid_to_formula_['S'] = MolecularFormula(5,  3, 1, 2, 0, 0);
  //amino_acid_to_formula['s'] = MolecularFormula(4,  3, 1, 5, 1, 0); # phosphate
  amino_acid_to_formula_['T'] = MolecularFormula(7,  4, 1, 2, 0, 0);
  //amino_acid_to_formula['t'] = MolecularFormula(6,  4, 1, 5, 1, 0); # phosphate
  amino_acid_to_formula_['V'] = MolecularFormula(9,  5, 1, 1, 0, 0);
  amino_acid_to_formula_['W'] = MolecularFormula(10, 11, 2, 1, 0, 0);
  amino_acid_to_formula_['Y'] = MolecularFormula(9,  9,  1, 2, 0, 0);
  //amino_acid_to_formula['y'] = MolecularFormula(8,  9,  1, 5, 1, 0);  # phosphate
  map_generated_ = true;
}

double AccurateMonoMass::calculateMass(const string& sequence) {

  if (!map_generated_) generateMap();
  MolecularFormula formula;

  for (int i=0;i<sequence.length();i++) {

    formula.add(amino_acid_to_formula_[sequence[i]]);

  }
  formula.add(H2O);


  return formula.calculateMass();
}

int main(int argc, char**argv) {
  
  AccurateMonoMass::generateMap();

  map<char, MolecularFormula>::iterator iter;

  cout << fixed << setprecision(11);

  for (iter = amino_acid_to_formula_.begin();
    iter != amino_acid_to_formula_.end();
    ++iter) {
    cout << "Mass of :"<<iter -> first <<" = " << iter -> second.calculateMass() << endl;
  }

  double diff = amino_acid_to_formula_['c'].calculateMass() - amino_acid_to_formula_['C'].calculateMass();

  cout <<"Diff of acetomide:"<<diff<<endl;

  cout <<"Mass of AA:"<<AccurateMonoMass::calculateMass(string("AA"))<<endl;
  cout <<"Mass of A + A:"<<(AccurateMonoMass::calculateMass(string("A")) + AccurateMonoMass::calculateMass(string("A"))-H2O.calculateMass())<<endl;
  double massA = AccurateMonoMass::calculateMass("A");

  double massAA = AccurateMonoMass::calculateMass("AAAAAAA");

  double massApA = 7*massA - 6 * H2O.calculateMass();

  double ppm = fabs(massAA - massApA) / massAA * 1e6;

  cout <<"ppm:"<<ppm<<endl;

  cout <<"H2O:"<<H2O.calculateMass()<<endl;

}


