#include "Out2XML.h"
#include "Common/TPPVersion.h"

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc
  if(argc < 3) {
    cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << endl;
    cerr << " usage: Out2XML <path to directory with out files> <# of top hits to report [1, 10]>" 
	 << " (OPTIONS)" << endl 
	 << endl
	 << " OPTIONS:" << endl
         << "   -m: use monoisotopic precursor weight (default: setting specified in sequest.params)" << endl
         << "   -a: use average precursor weight (default: setting specified in sequest.params)" << endl
         << "   -M: maldi mode" << endl
	 << "   -all: output all peptides, don't filter out X containing peptides " << endl
	 << "   -pI: compute and report peptide pI's " << endl
	 << "   -P<path to -including- sequest.params file>: (default) <path to directory with out files>/sequest.params " << endl
	 << "   -E<enzyme> " << endl
	 << "     Where <enzyme> is: " << endl
	 << "       trypsin - Cut: KR, No Cut: P, Sense: C-term (default)" << endl
	 << "       ralphtrypsin - Cut: STKR, No Cut: P, Sense: C-term " << endl
	 << "       stricttrypsin - Cut: KR, No Cut: none, Sense: C-term " << endl
	 << "       argc - Cut: R, No Cut: P, Sense: C-term " << endl
	 << "       aspn - Cut: D, No Cut: none, Sense: N-term " << endl
	 << "       chymotrypsin - Cut: YWFL, No Cut: P, Sense: C-term " << endl
	 << "       cnbr - Cut: M, No Cut: P, Sense: C-term " << endl
	 << "       elastase - Cut: GVLIA, No Cut: P, Sense: C-term " << endl
	 << "       formicacid - Cut: D, No Cut: P, Sense: C-term " << endl
	 << "       gluc - Cut: DE, No Cut: P, Sense: C-term " << endl
	 << "       gluc_bicarb - Cut: E, No Cut: P, Sense: C-term " << endl
	 << "       iodosobenzoate - Cut: W, No Cut: terminal, Sense: C-term " << endl
	 << "       lysc - Cut: K, No Cut: P, Sense: C-term " << endl
	 << "       lysc-p - Cut: K, No Cut: none, Sense: C-term " << endl
	 << "       lysn - Cut: K, No Cut: none, Sense: N-term " << endl
	 << "       lysn_promisc - Cut: KR, No Cut: none, Sense: N-term " << endl
	 << "       nonspecific - Cut: all, No Cut: none, Sense: N/A " << endl
	 << "       pepsina - Cut: FL, No Cut: terminal, Sense: C-term " << endl
	 << "       protein_endopeptidase - Cut: P, No Cut: terminal, Sense: C-term " << endl
	 << "       staph_protease - Cut: E, No Cut: terminal, Sense: C-term " << endl
	 << "       tca - Cut: KR, No Cut: P, Sense: C-term " << endl
	 << "           - Cut: YWFM, No Cut: P, Sense: C-term " << endl
	 << "           - Cut: D, No Cut: none, Sense: N-term " << endl
	 << "       trypsin/cnbr - Cut: KR, No Cut: P, Sense: C-term " << endl
	 << "       trypsin_gluc - Cut: DEKR, No Cut: P, Sense: C-term " << endl
	 << "       trypsin_k - Cut: K, No Cut: P, Sense: C-term " << endl
	 << "       trypsin_r - Cut: R, No Cut: P, Sense: C-term " << endl;

    exit(1);
  }
  Boolean write_all = False;

  Out2XML* conv = new Out2XML(argv[1], atoi(argv[2]), argv, argc);
  conv->processData();
  delete conv;
}
