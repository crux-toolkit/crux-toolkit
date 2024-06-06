#include <ctype.h>
#include "masscalc.h"
#include <stdlib.h>

/*
 * Mass calculator class
 *     Calculates molecular masses based on atomic masses.
 *     Atomic masses come from http://www.unimod.org/unimod_help.html.
 */
masscalc::masscalc(bool n15)
{
#ifdef COMET_EXACT
	addMass("H", 1.0078250, 1.00794);
	addMass("O", 15.9949146, 15.9994);
	if (n15)
	  addMass("N", 15.000074, 15.00374);
	else
	  addMass("N", 14.003074, 14.00674);
	addMass("C", 12.0, 12.0107);
	addMass("S", 31.9720718, 32.066);
	addMass("P", 30.9737633, 30.973761);
#else
	addMass("H", 1.007825035, 1.00794);
	addMass("O", 15.99491463, 15.9994);
	if (n15)
	  addMass("N", 15.000074, 15.00374);
	else
	  addMass("N", 14.003074, 14.0067);
	addMass("C", 12.0, 12.0107);
	addMass("S", 31.9720707, 32.065);
	addMass("P", 30.973762, 30.973761);
#endif
}

double masscalc::calcMass(const char* _m, massType _t)
{
	double totalMass = 0.0;

	const char* pchAtom = _m;
	const char* pchCount = pchAtom;

	string atom;
	int count;

	while (*pchCount != '\0')
	{
		// Advance past atom name
		pchCount++;
		while (isalpha(*pchCount) && !isupper(*pchCount))
			pchCount++;

		// Get count, 1 if not present
		count = 1;
		if (isdigit(*pchCount))
			count = atoi(pchCount);

		// Add atomic mass * count
		atom.assign(pchAtom, pchCount - pchAtom);
		totalMass += getMass(atom.data(), _t) * count;

		// Advance past count, if there is one
		while (*pchCount != '\0' && !isalpha(*pchCount))
			pchCount++;

		pchAtom = pchCount;
	}

	return totalMass;
}

double masscalc::getMass(const char* _m, massType _t)
{
		map<string, massPair>::iterator it = m_masses.find(_m);
		if (it == m_masses.end())
			return 0.0;

		if (_t == monoisotopic)
			return it->second.monoisotopic;

		return it->second.average;
}

void masscalc::addMass(const char* _m, double _mono, double _ave)
{
	string atom = _m;
	massPair masses;
	masses.monoisotopic = _mono;
	masses.average = _ave;

	m_masses.insert(make_pair(atom, masses));
}
