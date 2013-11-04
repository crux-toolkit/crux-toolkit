#include "CPeriodicTable.h"
#include <stdio.h>
#include <cstring>

using namespace std;



void CPeriodicTable::add(const char* symbol, double mass) {
  element elem;
  elem.atomicNum = (int)table.size();
  elem.mass = mass;
  table.push_back(elem);
  strcpy(table.at(table.size()-1).symbol, symbol);
}

void CPeriodicTable::init() {
  //cerr << "CPeriodicTable::init()"<<endl;
  table.clear();
  add("X", 0);
  add("H", 1.0078246);
  add("He", 3.01603);
  add("Li", 6.015121);
  add("Be", 9.012182);
  add("B", 10.012937);
  add("C", 12.0000000);
  add("N", 14.0030732);
  add("O", 15.9949141);
  add("F", 18.9984032);
  add("Ne", 19.992435);
  add("Na", 22.989767);
  add("Mg", 23.985042);
  add("Al", 26.981539);
  add("Si", 27.976927);
  add("P",  30.973762);
  add("S",  31.972070);
  add("Cl", 34.9688531);
  add("Ar", 35.967545);
  add("K",  38.963707);
  add("Ca", 39.962591);
  add("Sc", 44.955910);
  add("Ti", 45.952629);
  add("V",  49.947161);
  add("Cr", 49.946046);
  add("Mn", 54.938047);
  add("Fe", 53.939612);
  add("Co", 58.933198);
  add("Ni", 57.935346);
  add("Cu", 62.939598);
  add("Zn", 63.929145);
  add("Ga", 68.925580);
  add("Ge", 69.924250);
  add("As", 74.921594);
  add("Se", 73.922475);
  add("Br", 78.918336);
  add("Kr", 77.914);
  add("Rb", 84.911794);
  add("Sr", 83.913430);
  add("Y",  88.905849);
  add("Zr", 89.904703);
  add("Nb", 92.906377);
  add("Mo", 91.906808);
  add("Tc", 98.0);
  add("Ru", 95.907599);
  add("Rh", 102.905500);
  add("Pd", 101.905634);
  add("Ag", 106.905092);
  add("Cd", 105.906461);
  add("In", 112.904061);
  add("Sn", 111.904826);
  add("Sb", 120.903821);
  add("Te", 119.904048);
  add("I",  126.904473);
  add("Xe", 123.905894);
  add("Cs", 132.905429);
  add("Ba", 129.906282);
  add("La", 137.90711);
  add("Ce", 135.907140);
  add("Pr", 140.907647);
  add("Nd", 141.907719);
  add("Pm", 145.0);
  add("Sm", 143.911998);
  add("Eu", 150.919847);
  add("Gd", 151.919786);
  add("Tb", 158.925342);
  add("Dy", 155.925277);
  add("Ho", 164.930319);
  add("Er", 161.928775);
  add("Tm", 168.934212);
  add("Yb", 167.933894);
  add("Lu", 174.940770);
  add("Hf", 173.940044);
  add("Ta", 179.947462);
  add("W", 179.946701);
  add("Re", 184.952951);
  add("Os", 183.952488);
  add("Ir", 190.960584);
  add("Pt", 189.959917);
  add("Au", 196.966543);
  add("Hg", 195.965807);
  add("Tl", 202.972320);
  add("Pb", 203.973020);
  add("Bi", 208.980374);
  add("Po", 209.0);
  add("At", 210.0);
  add("Rn", 222.0);
  add("Fr", 223.0);
  add("Ra", 226.025);
  add("Ac", 227.028);
  add("Th", 232.038054);
  add("Pa", 231.0359);
  add("U",  234.040946);
  add("Np", 237.048);
  add("Pu", 244.0);
  add("Am", 243.0);
  add("Cm", 247.0);
  add("Bk", 247.0);
  add("Cf", 251.0);
  add("Es", 252.0);
  add("Fm", 257.0);
  add("Md", 258.0);
  add("No", 259.0);
  add("Lr", 260.0);
  add("Hx", 1.0078246);
  add("Cx", 12.0000000);
  add("Nx", 14.0030732);
  add("Ox", 15.9949141);
  add("Sx", 31.972070);
}

CPeriodicTable::CPeriodicTable(){

  init();
}



CPeriodicTable::CPeriodicTable(char* c){
  if (c == NULL || c[0] == '\0') {
    init();
  } else {
    loadTable(c);
  }
}

CPeriodicTable::~CPeriodicTable(){
}

element& CPeriodicTable::at(int i){
  return table.at(i);
}

 void CPeriodicTable::loadTable(char* c){
  table.clear();
  FILE* fptr;
  element e;
  
  fptr = fopen(c,"rt");
  if(fptr==NULL) {
    cerr << "Cannot open periodic table!" << endl;
    cerr << "Path:"<<c<<endl;
    cerr << "Initializing with defaults"<< endl;
    init();
    return;
  }
  
  /* 
     This loop reads in table entries, line by line.
     It has basic error detection (missing data), but additional
     checks should be made to confirm the content of data
  */
  while(!feof(fptr)){

		fscanf(fptr,"%d\t%s\t%lf\n",&e.atomicNum,e.symbol,&e.mass);   
    table.push_back(e);
    
  }

  fclose(fptr);
  
}

int CPeriodicTable::size(){
  return table.size();
}
