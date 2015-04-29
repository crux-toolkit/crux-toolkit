#ifndef _CPERIODICTABLE_H
#define _CPERIODICTABLE_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

typedef struct {
  int atomicNum;
  char symbol[3];
  double mass;
} element;

class CPeriodicTable {
 public:
   //Constructors & Destructors
   //CPeriodicTable();
   CPeriodicTable(char* c="Hardklor.dat");
   ~CPeriodicTable();

   //Methods:
   element& at(int);
   int size();

 protected:
 private:
   //Methods:
   void defaultValues();
   void loadTable(char*);

   //Data Members:
   vector<element> table;

};

#endif
