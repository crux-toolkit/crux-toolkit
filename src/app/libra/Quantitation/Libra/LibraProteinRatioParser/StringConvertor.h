#ifndef SCONVERTOR_H
#define SCONVERTOR_H

#include <string>

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
std::string ftos(double, int); 

/**
* convert integer to string:
*/
std::string itos(int);

#endif 
