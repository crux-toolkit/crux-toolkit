// Tandem2xml.h
//     Utility functions for converting X!Tandem output to
//     pepXML.

#include <string>
#include "Common/util.h"
using namespace std;

#ifndef _countof
#define _countof(array) (sizeof(array) / sizeof(array[0]))
#endif

extern const double dProtonMass;

// Output helpers
string& nl();
string& nlIn();
string& nlOut();

