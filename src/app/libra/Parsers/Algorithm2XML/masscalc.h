#include <string>
#include <map>

using namespace std;

/*
 * Mass calculator class
 *     Calculates molecular masses based on atomic masses.
 *     Atomic masses come from http://www.unimod.org/unimod_help.html.
 */
class masscalc
{
public:
	enum massType
	{	monoisotopic, average };

public:
  	masscalc(bool n15 = false);

	double calcMass(const char* _m, massType _t = monoisotopic);

private:
	void addMass(const char* _m, double _mono, double _ave);
	double getMass(const char* _m, massType _t);

private:
	struct massPair
	{
		double monoisotopic;
		double average;
	};

	map<string, massPair> m_masses;
};
