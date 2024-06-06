#include "RTCalculator.h"
#include "Common/TPPVersion.h"
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc
  
  //TODO: Add input error handling
  if(argc < 2) {
    cerr <<  argv[0] << " (" << szTPPVersionInfo << ")" << endl;
    cerr << "USAGE: " << argv[0] <<" CV ANN=[neural_net_file] COEFF=[coeff_file] TRAIN=[training_file] WITHCAT=[RTCatalog] CATONLY=[RTCatalog] TRAINCAT=[training_RTCatalog] PEPS=[peptide_file] PEPXML=[pepxml_file] OUTFILE=[output_results_here]" << endl;
    exit(1);
  }
  
  string coeffile = "";

  RTCalculator* rtcalc =  new RTCalculator(0);
  int testarg = 0;
  bool train = false;
  bool traincat = false;
  bool ann = false;
  bool catonly = false;
  bool withcat = false;
  bool cv = false;

  string pepxmlfile = "";
  string trainfile = "";
  string traincatfile = "";

  string catonlyfile = "";

  string withcatfile = "";

  string annfile = "";
  string testfile = "";

  string outfile = "";
  
  ostream* outstream = &cout;
  istream* pepstream = &cin;
  ofstream pepoutfile;
  ifstream pepinfile;
  //TODO Param passing needs work!!!
  for (int argidx = testarg+1; argidx < argc; argidx++) {
    string arg = argv[argidx];

    if (!strncmp(argv[argidx], "COEFF=", 6)) {
      coeffile = string(argv[argidx] + 6);
    }

    if (!strncmp(argv[argidx], "ANN=", 4)) {
      annfile = string(argv[argidx] + 4);
      ann=true;
    }


    if (!strncmp(argv[argidx], "CV", 2)) {
      cv=true;
    }

    if (!strncmp(argv[argidx], "TRAIN=", 6)) {
      trainfile = string(argv[argidx] + 6);
      train = true;
    }
    if (!strncmp(argv[argidx], "TRAINCAT=", 9)) {
      traincatfile = string(argv[argidx] + 9);
      traincat = true;
    }
    if (!strncmp(argv[argidx], "CATONLY=", 8)) {
      catonlyfile = string(argv[argidx] + 8);
      catonly = true;
    }
    if (!strncmp(argv[argidx], "WITHCAT=", 8)) {
      withcatfile = string(argv[argidx] + 8);
      withcat = true;
    }
    if (!strncmp(argv[argidx], "PEPS=", 5)) {
      testfile = string(argv[argidx] + 5);
    }
    if (!strncmp(argv[argidx], "PEPXML=", 7)) {
      pepxmlfile = string(argv[argidx] + 7);
    }
    
    if (!strncmp(argv[argidx], "OUTFILE=", 8)) {
      outfile = string(argv[argidx] + 8);
    }

  }
  

  if (train && !ann) {
    ifstream fin(trainfile.c_str());
    if(! fin) {
      cerr << "cannot read training file " << trainfile << endl;
      exit(1);
    }
    
    rtcalc->train_RTcoeff(fin); 
    // rtcalc->learnNeuralNet();
    
    rtcalc->linearRegressRT();
    rtcalc->recalc_RTstats();
    rtcalc->write_RTstats(cout);
    rtcalc->write_RTcoeff(cout);

    if (!coeffile.empty()) {
      ofstream fout(coeffile.c_str());
      if(! fout) {
	cerr << "cannot write coefficients file " << coeffile << endl;
	exit(1);
      }
      
      rtcalc->write_RTcoeff(fout);
      fout.close();
    }
    return 0;

  } 
  else if (traincat && !ann) {

    
    rtcalc->train_RTCatalog(traincatfile.c_str()); 
    // rtcalc->learnNeuralNet();
    
    rtcalc->linearRegressRT();
    rtcalc->recalc_RTstats();
    rtcalc->write_RTstats(cout);
    rtcalc->write_RTcoeff(cout);

    if (!coeffile.empty()) {
      ofstream fout(coeffile.c_str());
      if(! fout) {
	cerr << "cannot write coefficients file " << coeffile << endl;
	exit(1);
      }
      
      rtcalc->write_RTcoeff(fout);
      fout.close();
    }
    return 0;

  }
  else if (train && ann) {
    ifstream fin(trainfile.c_str());
    if(! fin) {
      cerr << "cannot read training file " << trainfile << endl;
      exit(1);
    }
    
    rtcalc->train_RTcoeff(fin); 
    if (!cv) {
      rtcalc->learnNeuralNet2();
    }
    else {
      rtcalc->learnNeuralNetCV();
    }

    //rtcalc->learnNeuralNet();
    
    //rtcalc->linearRegressRT();
    rtcalc->recalc_RTstats();
    rtcalc->write_RTstats(cout);
    
    rtcalc->write_RTann(annfile);
    rtcalc->write_paramIndex(cout);

    if (!annfile.empty()) {
      ofstream fout(annfile.c_str());
      if(! fout) {
	cerr << "cannot write neural net file " << annfile << endl;
	exit(1);
      }
      
      rtcalc->write_paramIndex(fout);
      fout.close();
    }
    return 0;

  }
  else if (traincat && ann) {

    rtcalc->train_RTCatalog(traincatfile.c_str()); 
    if (!cv) {
      rtcalc->learnNeuralNet2();
    }
    else {
      rtcalc->learnNeuralNetCV();
    }
    //rtcalc->learnNeuralNet();
    
    //rtcalc->linearRegressRT();
    rtcalc->recalc_RTstats();
    rtcalc->write_RTstats(cout);
    
    rtcalc->write_RTann(annfile);
    rtcalc->write_paramIndex(cout);

    if (!annfile.empty()) {
      ofstream fout(annfile.c_str());
      if(! fout) {
	cerr << "cannot write neural net file " << annfile << endl;
	exit(1);
      }
      
      rtcalc->write_paramIndex(fout);
      fout.close();
    }
    return 0;

  }
  else if (catonly || withcat) {
    if (catonly) {
      rtcalc->read_RTCatalog(catonlyfile.c_str());
    }
    else {
      rtcalc->read_RTCatalog(withcatfile.c_str());
    }
  }
  
  if (!train && !ann) {
    if (coeffile.empty()) {
      coeffile = getConfPath() + (string)"RTCalc.coeff";
    }
    ifstream fin(coeffile.c_str());
    if(! fin) {
      cerr << "cannot read coefficients file " << coeffile << endl;
      exit(1);
    }
    
    rtcalc->read_RTcoeff(fin); 

  }

  else if (!train) {
   if (annfile.empty()) {
      annfile = getConfPath() + (string)"RTCalc.ann";
    }
    ifstream fin(annfile.c_str());
    if(! fin) {
      cerr << "cannot read neural net param index file " << annfile << endl;
      exit(1);
    }
    
    rtcalc->read_paramIndex(fin); 

  }


  if (!pepxmlfile.empty()) {
    //TODO: DDS Implement rtcalc->parse_pepXML(pepxmlfile);
    return 0;
  }

  string peptide;
  
  if (!testfile.empty()) {
    pepinfile.open(testfile.c_str());
    if(! pepinfile) {
      cerr << "cannot read peptide file " << testfile << endl;
      exit(1);
    }
    pepstream = &pepinfile;
    
  }

  
  if (!outfile.empty()) {
    pepoutfile.open(outfile.c_str());
    if(! pepoutfile) {
      cerr << "cannot write output file " << testfile << endl;
      exit(1);
    }
    outstream = &pepoutfile;
  }
 

  while (1) {
    (*pepstream) >> peptide;

    if (pepstream->fail()) {
      break;
    }
    double RTout = -1;
    if (catonly) {
      (*outstream) << peptide << "\t" << rtcalc->get_PepRTCatalog(peptide) << "\n";
    }
    else if (withcat) {
      RTout = rtcalc->get_PepRTCatalog(peptide);
      if (RTout < -1000000000 || RTout > 1000000000) {
	if (!ann) {
	  RTout = rtcalc->calc_PepRT(peptide);	  
	}
	else {
	  RTout = rtcalc->calcANN_PepRT2(peptide);
	}
      }
       (*outstream) << peptide << "\t" << RTout << "\n";
    }
    else if (!ann) {
      (*outstream) << peptide << "\t" << rtcalc->calc_PepRT(peptide) << "\n";
    }
    else {
      //      cout << peptide << "\t" << rtcalc->calcANN_PepRT(peptide) << "\n";
      (*outstream) << peptide << "\t" << rtcalc->calcANN_PepRT2(peptide) << "\n";
    }
    //compute peptide's RT

  }


  return 0;
}
