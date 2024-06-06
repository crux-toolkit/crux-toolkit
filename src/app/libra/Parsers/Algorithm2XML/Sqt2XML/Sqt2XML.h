#ifndef SQT2XML_H_
#define SQT2XML_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>

using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::vector;
using std::getline;
using boost::lexical_cast;

typedef string L_entry;

struct M_entry
{
	int rankByXcorr;
	int rankBySp;
	float calculatedMassPlus1;
	float deltacn;
	float xcorr;
	float sp;
	int matchedIons;
	int predictedIons;
	string sequence;
	char validation;

	vector<L_entry> loci;
};

struct S_entry
{
	int firstScan;
	int lastScan;
	int chargeState;
	float processingTime;
	string processingHostname;
	float observedMassPlus1;
	float totalIntensity;
	float lowestSp;
	int numSequenceComparisons;

	vector<M_entry> matches;
};

// optimized specializations of lexical_cast for string->number conversions
namespace boost
{
	template<>
	inline float lexical_cast( const std::string& str )
	{
		return (float) strtod( str.c_str(), NULL );
	}

	template<>
	inline double lexical_cast( const std::string& str )
	{
		return strtod( str.c_str(), NULL );
	}

	template<>
	inline int lexical_cast( const std::string& str )
	{
		return atoi( str.c_str() );
	}

	template<>
	inline long lexical_cast( const std::string& str )
	{
		return atol( str.c_str() );
	}

	template<>
	inline unsigned int lexical_cast( const std::string& str )
	{
		return (unsigned int) atoi( str.c_str() );
	}

	template<>
	inline unsigned long lexical_cast( const std::string& str )
	{
		return (unsigned long) atol( str.c_str() );
	}
}

#endif
