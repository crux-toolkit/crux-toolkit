/*
Program       : StringConvertor
Author        : Nichole King      
Date          : 10.16.05 

Handy functions for string type conversion

*/

#include "StringConvertor.h"
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace std;

/**
* convert integer to string 
*/
string itos(int is)
{
    stringstream ss;

    string str;

    ss << is;

    str=ss.str();

    return str;
}

/**
*  convert double to string - choose precision and 
* scientifc or fixed notation .  
* (1) scientific notation with precision of 1
* (2) fixed point notation with precision of 1
* (3) fixed point notation with precision of 2
* (4) fixed point notation with precision of 3
* @param double number
* @param pattern for format
*/
string ftos(double fs, int sciorfix)
{
    stringstream ss;

    string str;

    if (sciorfix == 0) 
    {
       ss << setprecision(2) << setiosflags(ios::scientific) << fs;
    } else if (sciorfix == 1) 
    {
       ss << setprecision(1) << setiosflags(ios::fixed) << fs;
    } else if (sciorfix == 2) 
    {
       ss << setprecision(2) << setiosflags(ios::fixed) << fs;
    } else 
    {
       ss << setprecision(3) << setiosflags(ios::fixed) << fs;
    }

    str=ss.str();

    return str;
}
