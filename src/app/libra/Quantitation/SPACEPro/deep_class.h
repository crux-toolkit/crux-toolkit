#ifndef deep_class_h
#define deep_class_h

#include "deep_structs.h"
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector> 
#include <algorithm>

class peptide_lists {
public:
	int total_psms_in_xml;
	std::vector<dsPSM> all_psm; 
	std::vector<dsPeptide> all_peptides; 
	std::vector<dsProtein> all_proteins;

	//Functions
	void calc(metrics& my_metrics);
	bool delete_dup();
	bool enzymatic_calc(my_parameters& my_params);
	void json(std::string fn);
	bool miss_cleave(my_parameters& my_params);
	void print(metrics& my_metrics, my_parameters& params, std::string marquee);
	bool prot_stats();
	bool reader();
	void xml_parse(my_parameters& my_params);

private:
	void addString(std::vector<std::string>& v, const char* str){
		std::string s=str;
		addString(v,s);
	}
  void addString(std::vector<std::string> & v, std::string& str) {v.push_back(str);}
	bool cleanNoise(std::vector<dsXIC>& v);
	float calcPeakArea(std::vector<dsXIC>& v);
	dsPeptide convert_PSM_to_peptide(dsPSM& psm);

	static bool compareInten(const dsPair& a, const dsPair& b) { return a.ft_areaXIC > b.ft_areaXIC; };
	static bool compareSeqZ(const dsPSM& a, const dsPSM& b) {
		int i = a.pep_seq.compare(b.pep_seq);
		if (i == 0) return (a.charge < b.charge);
		else return (i < 0);
	}
	static bool compareTrypIndex(const dsPair& a, const dsPair& b) { return a.trypIndex < b.trypIndex; }
	static bool compareMissIndex(const dsPair& a, const dsPair& b) { return a.missIndex < b.missIndex; }
	static bool compareTotal(const dsProtein& a, const dsProtein& b) { return a.total > b.total; }
	static bool compareProt(const dsPeptide& a, const dsPeptide& b) { return a.prot_seq < b.prot_seq; }

};


#endif