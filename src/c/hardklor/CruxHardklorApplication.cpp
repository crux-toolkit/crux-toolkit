/**
 * \file CruxHardklorApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "CruxHardklorApplication.h"
#include "DelimitedFileWriter.h"

using namespace std;

/**
 * \returns a blank CruxHardklorApplication object
 */
CruxHardklorApplication::CruxHardklorApplication() {

}

/**
 * Destructor
 */
CruxHardklorApplication::~CruxHardklorApplication() {
}

/**
 * main method for CruxHardklorApplication
 */
int CruxHardklorApplication::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "hardklor-algorithm",
    "cdm",
    "min-charge",
    "max-charge",
    "corr",
    "depth",
    "distribution-area",
    "averagine-mod",
    "mzxml-filter",
    "no-base",
    "max-p",
    "resolution",
    "instrument",
    "centroided",
    "scan-number",
    "sensitivity",
    "signal-to-noise",
    "sn-window",
    "mz-window",
    "max-width",
    "parameter-file",
    "verbosity"
  };

  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"spectra"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  string input_spectra(get_string_parameter("spectra"));

  return main(input_spectra);

}

/**
 * \brief runs hardklor on the input spectra
 * \returns whether hardklor was successful or not
 */
int CruxHardklorApplication::main(
  const string& input_spectra ///< file path of spectra to process.
  ) {

  /* Write the dat files */
  string output_filename = make_file_path("hardklor.mono.txt");

  // TODO remove this dependency.
  string isotope_dat_filename = make_file_path("ISOTOPE.DAT");
  writeIsotopeDat(isotope_dat_filename);

  /* build argument list */
  vector<string> hk_args_vec;
  hk_args_vec.push_back("hardklor");
  hk_args_vec.push_back(input_spectra);
  hk_args_vec.push_back(output_filename);

  //Add options

  //TODO remove this dependency.
  hk_args_vec.push_back("-mdat");
  hk_args_vec.push_back(isotope_dat_filename);

  hk_args_vec.push_back("-a");

  HARDKLOR_ALGORITHM_T hk_algorithm = get_hardklor_algorithm("hardklor-algorithm");
  char* hk_algorithm_str = 
    hardklor_hardklor_algorithm_type_to_string(hk_algorithm);

  hk_args_vec.push_back(hk_algorithm_str);
  free(hk_algorithm_str);

  hk_args_vec.push_back("-cdm");
  hk_args_vec.push_back(get_string_parameter_pointer("cdm"));

  hk_args_vec.push_back("-chMin");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_int_parameter("min-charge")));

  hk_args_vec.push_back("-chMax");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_int_parameter("max-charge")));

  hk_args_vec.push_back("-corr");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_double_parameter("corr")));

  hk_args_vec.push_back("-d");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_int_parameter("depth")));

  hk_args_vec.push_back("-da");

  if (get_boolean_parameter("distribution-area")) {
    hk_args_vec.push_back("true");
  } else {
    hk_args_vec.push_back("false");
  }

  if (string(get_string_parameter_pointer("averagine-mod")) !=  "__NULL_STR") {
    string par_value = get_string_parameter_pointer("averagine-mod");
    carp(CARP_DEBUG,"averagine-mod=%s",par_value.c_str());
    vector<string> tokens1;
    tokenize(par_value, tokens1, ';');

    for (unsigned int idx1 = 0; idx1 < tokens1.size();idx1++) {
      hk_args_vec.push_back("-m");
      vector<string> tokens2;
      tokenize(tokens1[idx1], tokens2, ':');
      ostringstream oss;
      oss << tokens2[0];
      for (unsigned int idx2 = 1;idx2 < tokens2.size(); idx2++) {
        oss << " " << tokens2[idx2];
      }
      hk_args_vec.push_back(oss.str());
    }
  }

  if (string(get_string_parameter_pointer("mzxml-filter")) != "none") { 
    hk_args_vec.push_back("-mF");
    hk_args_vec.push_back(get_string_parameter_pointer("mzxml-filter"));
  }

  if (get_boolean_parameter("no-base")) {
    hk_args_vec.push_back("-nb");
  }  
  
  hk_args_vec.push_back("-p");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_int_parameter("max-p")));

  hk_args_vec.push_back("-res");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_double_parameter("resolution")));
  hk_args_vec.push_back(get_string_parameter("instrument"));
  
  if (get_boolean_parameter("centroided")) {
    hk_args_vec.push_back("-c true");
  } else {
    hk_args_vec.push_back("-c false");
  }

  if (string(get_string_parameter_pointer("scan-number")) != "__NULL_STR") {
    const char* scan_numbers=get_string_parameter_pointer("scan-number");
    int first_scan;
    int last_scan;
    get_range_from_string(scan_numbers, first_scan, last_scan);
    
    hk_args_vec.push_back("-sc");
    hk_args_vec.push_back(DelimitedFileWriter::to_string(first_scan));
    hk_args_vec.push_back(DelimitedFileWriter::to_string(last_scan));
  }

  hk_args_vec.push_back("-sl");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_int_parameter("sensitivity")));

  hk_args_vec.push_back("-sn");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_double_parameter("signal-to-noise")));

  hk_args_vec.push_back("-snWin");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_double_parameter("sn-window")));
  
  if (get_boolean_parameter("static-sn")) {
    hk_args_vec.push_back("true");
  } else {
    hk_args_vec.push_back("false");
  }

  if (string(get_string_parameter_pointer("mz-window")) != "__NULL_STR") {
    double first;
    double last;
    if (get_range_from_string(get_string_parameter_pointer("mz-window"), first, last)) {
      hk_args_vec.push_back("-w");
      hk_args_vec.push_back(DelimitedFileWriter::to_string(first));
      hk_args_vec.push_back(DelimitedFileWriter::to_string(last));
    }
  }    

  hk_args_vec.push_back("-win");
  hk_args_vec.push_back(DelimitedFileWriter::to_string(get_double_parameter("max-width")));
  
  if (string(get_string_parameter_pointer("hardklor-options")) != "__NULL_STR") {
    string hardklor_options = get_string_parameter_pointer("hardklor-options");
    vector<string> tokens;
    tokenize(hardklor_options, tokens, ',');
    for (size_t idx=0;idx < tokens.size();idx++) {
      hk_args_vec.push_back(tokens.at(idx));
    }
  }

  /* build argv line */
  int hk_argc = hk_args_vec.size();

  char** hk_argv = new char*[hk_argc];

  hk_argv[0] = (char*)hk_args_vec[0].c_str();
  for (int idx = 1;idx < hk_argc ; idx++) {
    hk_argv[idx] = (char*)hk_args_vec[idx].c_str();
    carp(CARP_DEBUG, "hk_argv[%d]=%s", idx, hk_argv[idx]);
  }

  /* Call hardklorMain */
  int ret = hardklorMain(hk_argc, hk_argv);

  delete []hk_argv;

  return ret;
}

/**
 * \returns the command name for CruxHardklorApplication
 */
string CruxHardklorApplication::getName() {
  return "hardklor";
}

/**
 * \returns the description for CruxHardklorApplication
 */
string CruxHardklorApplication::getDescription() {

  return "Identify isotopic distributions from high-resolution mass spectra.";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CruxHardklorApplication::needsOutputDirectory() {
  return true;
}

/**
 * writes the ISOTOPE.DAT file for hardklor
 */
void CruxHardklorApplication::writeIsotopeDat(
  string& filename ///<path for dat file
) {

  ofstream fout(filename.c_str());
  writeIsotopeDat(fout);
  fout.close();
}

/**
 * write the ISOTOPE.DAT to an output stream
 */
void CruxHardklorApplication::writeIsotopeDat(
  ostream& os ///< the output stream to use
  ) {

  os <<"X  2" << endl;
  os <<"1  0.9" << endl;
  os <<"2  0.1" << endl;
  os << endl;
  os <<"H  2" << endl;
  os <<"1.0078246  0.999855"<<endl;
  os <<"2.0141021  0.000145"<<endl;
  os << endl;
  os <<"He  2"<<endl;
  os <<"3.01603    0.00000138"<<endl;
  os <<"4.00260    0.99999862"<<endl;
  os <<endl;
  os <<"Li  2" << endl;
  os <<"6.015121   0.075"<<endl;
  os <<"7.016003   0.925"<<endl;
  os << endl;
  os <<"Be  1" << endl;
  os <<"9.012182   1.0"<<endl;
  os << endl;
  os <<"B  2"<<endl;
  os <<"10.012937  0.199"<<endl;
  os <<"11.009305  0.801"<<endl;
  os << endl;
  os <<"C  2"<<endl;
  os <<"12.0000000 0.98916"<<endl;
  os <<"13.0033554 0.01084"<<endl;
  os << endl;
  os <<"N  2"<<endl;
  os <<"14.0030732 0.99633"<<endl;
  os <<"15.0001088 0.00366"<<endl;
  os << endl;
  os <<"O  3"<<endl;
  os <<"15.9949141 0.997576009706"<<endl;
  os <<"16.9991322 0.000378998479"<<endl;
  os <<"17.9991616 0.002044991815"<<endl;
  os << endl;
  os <<"F  1"<<endl;
  os <<"18.9984032 1.0"<<endl;
  os << endl;
  os <<"Ne  3"<<endl;
  os <<"19.992435  0.9048"<<endl;
  os <<"20.993843  0.0027"<<endl;
  os <<"21.991383  0.0925"<<endl;
  os << endl;
  os << "Na  1"<<endl;
  os <<"22.989767  1.0"<<endl;
  os << endl;
  os <<"Mg  3"<<endl;
  os <<"23.985042  0.7899"<<endl;
  os <<"24.985837  0.1000"<<endl;
  os <<"25.982593  0.1101"<<endl;
  os << endl;
  os <<"Al  1"<<endl;
  os <<"26.981539  1.0"<<endl;
  os << endl;
  os <<"Si  3"<<endl;
  os <<"27.976927  0.9223"<<endl;
  os <<"28.976495  0.0467"<<endl;
  os <<"29.973770  0.0310"<<endl;
  os << endl;
  os <<"P  1"<<endl;
  os <<"30.973762  1.0"<<endl;
  os << endl;
  os <<"S  4"<<endl;
  os <<"31.972070  0.95021"<<endl;
  os <<"32.971456  0.00745"<<endl;
  os <<"33.967866  0.04221"<<endl;
  os <<"35.967080  0.00013"<<endl;
  os << endl;
  os <<"Cl  2"<<endl;
  os <<"34.9688531 0.755290"<<endl;
  os <<"36.9659034 0.244710"<<endl;
  os << endl;
  os <<"Ar  3"<<endl;
  os <<"35.967545  0.00337"<<endl;
  os <<"37.962732  0.00063"<<endl;
  os <<"39.962384  0.99600"<<endl;
  os << endl;
  os <<"K  3"<<endl;
  os <<"38.963707  0.932581"<<endl;
  os <<"39.963999  0.000117"<<endl;
  os <<"40.961825  0.067302"<<endl;
  os << endl;
  os <<"Ca  6"<<endl;
  os <<"39.962591  0.96941"<<endl;
  os <<"41.958618  0.00647"<<endl;
  os <<"42.958766  0.00135"<<endl;
  os <<"43.955480  0.02086"<<endl;
  os <<"45.953689  0.00004"<<endl;
  os <<"47.952533  0.00187"<<endl;
  os << endl;
  os <<"Sc  1"<<endl;
  os <<"44.955910  1.0"<<endl;
  os << endl;
  os <<"Ti  5" << endl;
  os <<"45.952629  0.080" <<endl;
  os <<"46.951764  0.073" <<endl;
  os <<"47.947947  0.738" <<endl;
  os <<"48.947871  0.055" <<endl;
  os <<"49.944792  0.054" <<endl;
  os <<endl;
  os <<"V  2" << endl;
  os <<"49.947161  0.00250"<<endl;
  os <<"50.943962  0.99750"<<endl;
  os <<endl;
  os <<"Cr  4" << endl;
  os <<"49.946046  0.04345"<<endl;
  os <<"51.940509  0.83790"<<endl;
  os <<"52.940651  0.09500"<<endl;
  os <<"53.938882  0.02365"<<endl;
  os <<endl;
  os <<"Mn  1"<<endl;
  os <<"54.938047  1.0"<<endl;
  os << endl;
  os <<"Fe  4"<<endl;
  os <<"53.939612  0.0590"<<endl;
  os <<"55.934939  0.9172"<<endl;
  os <<"56.935396  0.0210"<<endl;
  os <<"57.933277  0.0028"<<endl;
  os << endl;
  os <<"Co  1"<<endl;
  os <<"58.933198  1.0"<<endl;
  os << endl;
  os <<"Ni  5"<<endl;
  os <<"57.935346  0.6827"<<endl;
  os <<"59.930788  0.2610"<<endl;
  os << "60.931058  0.0113"<<endl;
  os << "61.928346  0.0359"<<endl;
  os << "63.927968  0.0091"<<endl;
  os << endl;
  os << "Cu  2" << endl;
  os << "62.939598  0.6917"<<endl;
  os << "64.927793  0.3083"<<endl;
  os << endl;
  os <<"Zn  5"<<endl;
  os <<"63.929145  0.486"<<endl;
  os <<"65.926034  0.279"<<endl;
  os <<"66.927129  0.041"<<endl;
  os <<"67.924846  0.188"<<endl;
  os <<"69.925325  0.006"<<endl;
  os << endl;
  os <<"Ga  2"<<endl;
  os <<"68.925580  0.60108"<<endl;
  os <<"70.924700  0.39892"<<endl;
  os << endl;
  os <<"Ge  5"<<endl;
  os <<"69.924250  0.205"<<endl;
  os <<"71.922079  0.274"<<endl;
  os <<"72.923463  0.078"<<endl;
  os <<"73.921177  0.365"<<endl;
  os <<"75.921401  0.078"<<endl;
  os << endl;
  os << "As  1"<<endl;
  os << "74.921594  1.0"<<endl;
  os << endl;
  os << "Se  6"<<endl;
  os <<"73.922475  0.009"<<endl;
  os <<"75.919212  0.091"<<endl;
  os <<"76.919912  0.076"<<endl;
  os <<"77.9190    0.236"<<endl;
  os <<"79.916520  0.499"<<endl;
  os <<"81.916698  0.089"<<endl;
  os << endl;
  os <<"Br  2"<<endl;
  os <<"78.918336  0.5069"<<endl;
  os <<"80.916289  0.4931"<<endl;
  os << endl;
  os <<"Kr  6"<<endl;
  os <<"77.914     0.0035"<<endl;
  os <<"79.916380  0.0225"<<endl;
  os <<"81.913482  0.116"<<endl;
  os <<"82.914135  0.115"<<endl;
  os <<"83.911507  0.570"<<endl;
  os <<"85.910616  0.173"<<endl;
  os << endl;
  os <<"Rb  2"<<endl;
  os <<"84.911794  0.7217"<<endl;
  os <<"86.909187  0.2783"<<endl;
  os << endl;
  os << "Sr  4"<<endl;
  os <<"83.913430  0.0056"<<endl;
  os <<"85.909267  0.0986"<<endl;
  os <<"86.908884  0.0700"<<endl;
  os <<"87.905619  0.8258"<<endl;
  os << endl;
  os <<"Y  1"<<endl;
  os <<"88.905849  1.0"<<endl;
  os << endl;
  os <<"Zr  5"<<endl;
  os <<"89.904703  0.5145"<<endl;
  os <<"90.905644  0.1122"<<endl;
  os <<"91.905039  0.1715"<<endl;
  os <<"93.906314  0.1738"<<endl;
  os <<"95.908275  0.0280"<<endl;
  os << endl;
  os <<"Nb  1"<<endl;
  os <<"92.906377  1.0"<<endl;
  os << endl;
  os <<"Mo  7"<<endl;
  os <<"91.906808  0.1484"<<endl;
  os <<"93.905085  0.0925"<<endl;
  os <<"94.905840  0.1592"<<endl;
  os <<"95.904678  0.1668"<<endl;
  os <<"96.906020  0.0955"<<endl;
  os <<"97.905406  0.2413"<<endl;
  os <<"99.907477  0.0963"<<endl;
  os << endl;
  os <<"Tc  1"<<endl;
  os <<"98.0   1.0"<<endl;
  os << endl;
  os << "Ru  7"<<endl;
  os <<"95.907599  0.0554"<<endl;
  os <<"97.905287  0.0186"<<endl;
  os <<"98.905939  0.127"<<endl;
  os <<"99.904219  0.126"<<endl;
  os <<"100.905582  0.171"<<endl;
  os <<"101.904348  0.316"<<endl;
  os <<"103.905424  0.186"<<endl;
  os << endl;
  os <<"Rh  1"<<endl;
  os <<"102.905500  1.0"<<endl;
  os << endl;
  os <<"Pd  6"<<endl;
  os <<"101.905634  0.0102"<<endl;
  os <<"103.904029  0.1114"<<endl;
  os <<"104.905079  0.2233"<<endl;
  os <<"105.903478  0.2733"<<endl;
  os <<"107.903895  0.2646"<<endl;
  os <<"109.905167  0.1172"<<endl;
  os << endl;
  os <<"Ag  2"<<endl;
  os <<"106.905092  0.51839"<<endl;
  os <<"108.904757  0.48161"<<endl;
  os << endl;
  os <<"Cd  8" << endl;
  os <<"105.906461  0.0125"<<endl;
  os <<"107.904176  0.0089"<<endl;
  os <<"109.903005  0.1249"<<endl;
  os <<"110.904182  0.1280"<<endl;
  os <<"111.902758  0.2413"<<endl;
  os <<"112.904400  0.1222"<<endl;
  os <<"113.903357  0.2873"<<endl;
  os <<"115.904754  0.0749"<<endl;
  os << endl;
  os <<"In  2" << endl;
  os <<"112.904061  0.043"<<endl;
  os <<"114.903880  0.957"<<endl;
  os << endl;
  os <<"Sn  10"<<endl;
  os <<"111.904826  0.0097"<<endl;
  os <<"113.902784  0.0065"<<endl;
  os <<"114.903348  0.0036"<<endl;
  os <<"115.901747  0.1453"<<endl;
  os <<"116.902956  0.0768"<<endl;
  os <<"117.901609  0.2422"<<endl;
  os <<"118.903310  0.0858"<<endl;
  os <<"119.902200  0.3259"<<endl;
  os <<"121.903440  0.0463"<<endl;
  os <<"123.905274  0.0579"<<endl;
  os << endl;
  os << "Sb  2" << endl;
  os <<"120.903821  0.574"<<endl;
  os <<"122.904216  0.426"<<endl;
  os << endl;
  os << "Te  8" << endl;
  os << "119.904048  0.00095" << endl;
  os << "121.903054  0.0259" << endl;
  os << "122.904271  0.00905" << endl;
  os << "123.902823  0.0479" << endl;
  os << "124.904433  0.0712" << endl;
  os << "125.903314  0.1893" << endl;
  os << "127.904463  0.3170" << endl;
  os << "129.906229  0.3387" << endl;
  os << endl;
  os << "I  1" << endl;
  os << "126.904473  1.0"<<endl;
  os << endl;
  os << "Xe  9" << endl;
  os << "123.905894  0.0010"<<endl;
  os << "125.904281  0.0009"<<endl;
  os << "127.903531  0.0191"<<endl;
  os << "128.904780  0.264" << endl;
  os << "129.903509  0.041" << endl;
  os << "130.905072  0.212" << endl;
  os << "131.904144  0.269" << endl;
  os << "133.905395  0.104" << endl;
  os << "135.907214  0.089" << endl;
  os << endl;
  os << "Cs  1" << endl;
  os << "132.905429  1.0" << endl;
  os << endl;
  os << "Ba  7" << endl;
  os << "129.906282  0.00106" << endl;
  os << "131.905042  0.00101" << endl;
  os << "133.904486  0.0242" << endl;
  os << "134.905665  0.06593" << endl;
  os << "135.904553  0.0785" << endl;
  os << "136.905812  0.1123" << endl;
  os << "137.905232  0.7170" << endl;
  os << endl;
  os << "La  2" << endl;
  os << "137.90711   0.00090" << endl;
  os << "138.906347  0.99910" << endl;
  os << endl;
  os << "Ce  4" << endl;
  os << "135.907140  0.0019" << endl;
  os << "137.905985  0.0025" << endl;
  os << "139.905433  0.8843" << endl;
  os << "141.909241  0.1113" << endl;
  os << endl;
  os << "Pr  1" << endl;
  os << "140.907647  1.0" << endl;
  os << endl;
  os << "Nd  7" << endl;
  os << "141.907719  0.2713"<<endl;
  os << "142.909810  0.1218"<<endl;
  os << "143.910083  0.2380"<<endl;
  os << "144.912570  0.0830"<<endl;
  os << "145.913113  0.1719"<<endl;
  os << "147.916889  0.0576"<<endl;
  os << "149.920887  0.0564"<<endl;
  os << endl;
  os << "Pm  1" << endl;
  os << "145.0  1.0" << endl;
  os << endl;
  os << "Sm  7" << endl;
  os << "143.911998  0.031"<<endl;
  os << "146.914895  0.150"<<endl;
  os << "147.914820  0.113"<<endl;
  os << "148.917181  0.138"<<endl;
  os << "149.917273  0.074"<<endl;
  os << "151.919729  0.267"<<endl;
  os << "153.922206  0.227"<<endl;
  os << endl;
  os << "Eu  2" << endl;
  os <<"150.919847  0.478"<<endl;
  os <<"152.921225  0.522"<<endl;
  os << endl;
  os << "Gd  7"<<endl;
  os <<"151.919786  0.0020"<<endl;
  os <<"153.920861  0.0218"<<endl;
  os <<"154.922618  0.1480"<<endl;
  os <<"155.922118  0.2047"<<endl;
  os <<"156.923956  0.1565"<<endl;
  os <<"157.924099  0.2484"<<endl;
  os <<"159.927049  0.2186"<<endl;
  os << endl;
  os <<"Tb  1"<<endl;
  os <<"158.925342  1.0"<<endl;
  os << endl;
  os <<"Dy  7"<<endl;
  os <<"155.925277  0.0006"<<endl;
  os <<"157.924403  0.0010"<<endl;
  os <<"159.925193  0.0234"<<endl;
  os <<"160.926930  0.189"<<endl;
  os <<"161.926795  0.255"<<endl;
  os <<"162.928728  0.249"<<endl;
  os <<"163.929171  0.282"<<endl;
  os << endl;
  os << "Ho  1"<<endl;
  os << "164.930319  1.0"<<endl;
  os << endl;
  os << "Er  6" << endl;
  os <<"161.928775  0.0014" << endl;
  os <<"163.929198  0.0161" << endl;
  os <<"165.930290  0.336" << endl;
  os <<"166.932046  0.2295" << endl;
  os <<"167.932368  0.268" << endl;
  os <<"169.935461  0.149" << endl;
  os << endl;
  os <<"Tm  1"<<endl;
  os <<"168.934212  1.0"<<endl;
  os << endl;
  os <<"Yb  7"<<endl;
  os <<"167.933894  0.0013"<<endl;
  os <<"169.934759  0.0305"<<endl;
  os <<"170.936323  0.143"<<endl;
  os <<"171.936378  0.219"<<endl;
  os <<"172.938208  0.1612"<<endl;
  os <<"173.938859  0.318"<<endl;
  os <<"175.942564  0.127"<<endl;
  os << endl;
  os <<"Lu  2"<<endl;
  os <<"174.940770  0.9741"<<endl;
  os <<"175.942679  0.0259"<<endl;
  os << endl;
  os << "Hf  6"<<endl;
  os <<"173.940044  0.00162"<<endl;
  os <<"175.941406  0.05206"<<endl;
  os <<"176.943217  0.18606"<<endl;
  os <<"177.943696  0.27297"<<endl;
  os <<"178.945812  0.13629"<<endl;
  os <<"179.946545  0.35100"<<endl;
  os << endl;
  os << "Ta  2"<<endl;
  os <<"179.947462  0.00012"<<endl;
  os <<"180.947992  0.99988"<<endl;
  os << endl;
  os << "W  5" << endl;
  os <<"179.946701  0.0012" << endl;
  os <<"181.948202  0.263" << endl;
  os <<"182.950220  0.1428" << endl;
  os <<"183.950928  0.307" << endl;
  os <<"185.954357  0.286" << endl;
  os << endl;
  os << "Re  2" << endl;
  os << "184.952951  0.3740" << endl;
  os << "186.955744  0.6260" << endl;
  os << endl;
  os << "Os  7" << endl;
  os << "183.952488  0.0002" << endl;
  os << "185.953830  0.0158" << endl;
  os << "186.955741  0.016" << endl;
  os << "187.955860  0.133" << endl;
  os << "188.958137  0.161" << endl;
  os << "189.958436  0.264" << endl;
  os << "191.961467  0.410" << endl;
  os << endl;
  os << "Ir  2" << endl;
  os << "190.960584  0.373" << endl;
  os << "192.962917  0.627" << endl;
  os << endl;
  os << "Pt  6" << endl;
  os << "189.959917  0.0001" << endl;
  os << "191.961019  0.0079" << endl;
  os << "193.962655  0.329" << endl;
  os << "194.964766  0.338" << endl;
  os << "195.964926  0.253" << endl;
  os << "197.967869  0.072" << endl;
  os << endl;
  os << "Au  1" << endl;
  os <<"196.966543  1.0" << endl;
  os << endl;
  os << "Hg  7" << endl;
  os << "195.965807  0.0015" << endl;
  os << "197.966743  0.100" << endl;
  os << "198.968254  0.169" << endl;
  os << "199.968300  0.231" << endl;
  os << "200.970277  0.132" << endl;
  os << "201.970617  0.298" << endl;
  os << "203.973467  0.0685" << endl;
  os << endl;
  os << "Tl  2" << endl;
  os << "202.972320  0.29524" << endl;
  os << "204.974401  0.70476" << endl;
  os << endl;
  os << "Pb  4" << endl;
  os <<"203.973020  0.014" << endl;
  os <<"205.974440  0.241" << endl;
  os <<"206.975872  0.221" << endl;
  os <<"207.976627  0.524" << endl;
  os << endl;
  os <<"Bi  1" << endl;
  os <<"208.980374  1.0"<<endl;
  os << endl;
  os <<"Po  1"<<endl;
  os <<"209.0  1.0"<<endl;
  os << endl;
  os << "At  1" << endl;
  os <<"210.0  1.0" << endl;
  os << endl;
  os << "Rn  1" << endl;
  os << "222.0  1.0" << endl;
  os << endl;
  os << "Fr  1" << endl;
  os << "223.0  1.0" << endl;
  os << endl;
  os <<"Ra  1" << endl;
  os <<"226.025  1.0"<<endl;
  os << endl;
  os << "Ac  1" << endl;
  os << "227.028  1.0"<<endl;
  os << endl;
  os <<"Th  1"<<endl;
  os <<"232.038054  1.0"<<endl;
  os << endl;
  os << "Pa  1" << endl;
  os << "231.0359  1.0" << endl;
  os << endl; 
  os << "U  3" << endl;
  os << "234.040946  0.000055" << endl;
  os << "235.043924  0.00720" << endl;
  os << "238.050784  0.992745" << endl;
  os << endl;
  os << "Np  1" << endl;
  os <<"237.048  1.0" << endl;
  os << endl;
  os << "Pu  1" << endl;
  os <<"244.0  1.0" << endl;
  os << endl;
  os << "Am  1" << endl;
  os <<"243.0  1.0" << endl;
  os << endl;
  os << "Cm  1" << endl;
  os << "247.0  1.0" << endl;
  os << endl;
  os << "Bk  1" << endl;
  os << "247.0  1.0" << endl;
  os << endl;
  os << "Cf  1" << endl;
  os << "251.0  1.0" << endl;
  os << endl;
  os << "Es  1" << endl;
  os << "252.0  1.0" << endl;
  os << endl;
  os << "Fm  1" << endl;
  os << "257.0  1.0" << endl;
  os << endl;
  os << "Md  1" << endl;
  os << "258.0  1.0" << endl;
  os << endl;
  os << "No  1" << endl;
  os << "259.0  1.0" << endl;
  os << endl;
  os << "Lr  1" << endl;
  os << "260.0  1.0" << endl;
  os << endl;
  os << "Hx  2" << endl;
  os << "1.0078246  0.999855" << endl;
  os << "2.0141021  0.000145" << endl;
  os << endl;
  os << "Cx  2" << endl;
  os << "12.0000000 0.98916" << endl;
  os << "13.0033554 0.01084" << endl;
  os << endl;
  os << "Nx  2" << endl;
  os <<"14.0030732 0.99633" << endl;
  os <<"15.0001088 0.00366" << endl;
  os << endl;
  os <<"Ox  3" << endl;
  os <<"15.9949141 0.997576009706" << endl;
  os <<"16.9991322 0.000378998479" << endl;
  os <<"17.9991616 0.002044991815" << endl;
  os << endl;
  os << "Sx  4" << endl;
  os <<"31.972070  0.95021"<<endl;
  os <<"32.971456  0.00745"<<endl;
  os <<"33.967866  0.04221"<<endl;
  os <<"35.967080  0.00013"<<endl;
}


/**
 * writes the Hardklor.dat to a path
 */
void CruxHardklorApplication::writeHardklorDat(
  string& filename ///<path to write the Hardklor.dat to
  ) {

  ofstream fout(filename.c_str());
  writeHardklorDat(fout);
  fout.close();
}

/**
 * writes the Hardklor.dat to a stream
 */
void CruxHardklorApplication::writeHardklorDat(
  ostream& os ///< stream to write to.
  ) {

  os << "0" << "\t" << "X"  << "\t" << "0"         << endl;
  os << "1" << "\t" << "H"  << "\t" << "1.0078246" << endl;
  os << "2" << "\t" << "He" << "\t" << "3.01603"   << endl; 
  os << "3" << "\t" << "Li" << "\t" << "6.015121"  << endl;
  os << "4" << "\t" << "Be" << "\t" << "9.012182"  << endl;
  os << "5" << "\t" << "B"  << "\t" << "10.012937"  << endl; 
  os << "6" << "\t" << "C"  << "\t" << "12.0000000" << endl;
  os << "7" << "\t" << "N"  << "\t" << "14.0030732" << endl;
  os << "8" << "\t" << "O"  << "\t" << "15.9949141" << endl; 
  os << "9" << "\t" << "F"  << "\t" << "18.9984032" << endl; 
  os <<"10" << "\t" << "Ne" << "\t" << "19.992435"  << endl;
  os <<"11" << "\t" << "Na" << "\t" << "22.989767"  << endl; 
  os <<"12" << "\t" << "Mg" << "\t" << "23.985042"  << endl; 
  os <<"13" << "\t" << "Al" << "\t" << "26.981539"  << endl; 
  os <<"14" << "\t" << "Si" << "\t" << "27.976927"  << endl; 
  os <<"15" << "\t" << "P"  << "\t" << "30.973762"  << endl; 
  os <<"16" << "\t" << "S"  << "\t" << "31.972070"  << endl;  
  os <<"17" << "\t" << "Cl" << "\t" << "34.9688531" << endl; 
  os <<"18" << "\t" << "Ar" << "\t" << "35.967545"  << endl; 
  os <<"19" << "\t" << "K"  << "\t" << "38.963707"  << endl;
  os <<"20" << "\t" << "Ca" << "\t" << "39.962591"  << endl;
  os <<"21" << "\t" << "Sc" << "\t" << "44.955910"  << endl; 
  os <<"22" << "\t" << "Ti" << "\t" << "45.952629"  << endl; 
  os <<"23" << "\t" << "V"  << "\t" << "49.947161"  << endl; 
  os <<"24" << "\t" << "Cr" << "\t" << "49.946046"  << endl; 
  os <<"25" << "\t" << "Mn" << "\t" << "54.938047"  << endl;
  os <<"26" << "\t" << "Fe" << "\t" << "53.939612"  << endl; 
  os <<"27" << "\t" << "Co" << "\t" << "58.933198"  << endl;
  os <<"28" << "\t" << "Ni" << "\t" << "57.935346"  << endl; 
  os <<"29" << "\t" << "Cu" << "\t" << "62.939598"  << endl; 
  os <<"30" << "\t" << "Zn" << "\t" << "63.929145"  << endl; 
  os <<"31" << "\t" << "Ga" << "\t" << "68.925580"  << endl; 
  os <<"32" << "\t" << "Ge" << "\t" << "69.924250"  << endl;
  os <<"33" << "\t" << "As" << "\t" << "74.921594"  << endl;
  os <<"34" << "\t" << "Se" << "\t" << "73.922475"  << endl; 
  os <<"35" << "\t" << "Br" << "\t" << "78.918336"  << endl; 
  os <<"36" << "\t" << "Kr" << "\t" << "77.914"     << endl;    
  os <<"37" << "\t" << "Rb" << "\t" << "84.911794"  << endl;
  os <<"38" << "\t" << "Sr" << "\t" << "83.913430"  << endl; 
  os <<"39" << "\t" << "Y"  << "\t" << "88.905849"  << endl;
  os <<"40" << "\t" << "Zr" << "\t" << "89.904703"  << endl;
  os <<"41" << "\t" << "Nb" << "\t" << "92.906377"  << endl;
  os <<"42" << "\t" << "Mo" << "\t" << "91.906808"  << endl;
  os <<"43" << "\t" << "Tc" << "\t" << "98.0"       << endl; 
  os <<"44" << "\t" << "Ru" << "\t" << "95.907599"  << endl;
  os <<"45" << "\t" << "Rh" << "\t" << "102.905500" << endl;
  os <<"46" << "\t" << "Pd" << "\t" << "101.905634" << endl;
  os <<"47" << "\t" << "Ag" << "\t" << "106.905092" << endl;
  os <<"48" << "\t" << "Cd" << "\t" << "105.906461" << endl;
  os <<"49" << "\t" << "In" << "\t" << "112.904061" << endl;
  os <<"50" << "\t" << "Sn" << "\t" << "111.904826" << endl;
  os <<"51" << "\t" << "Sb" << "\t" << "120.903821" << endl;
  os <<"52" << "\t" << "Te" << "\t" << "119.904048" << endl;
  os <<"53" << "\t" << "I"  << "\t" << "126.904473" << endl;
  os <<"54" << "\t" << "Xe" << "\t" << "123.905894" << endl;
  os <<"55" << "\t" << "Cs" << "\t" << "132.905429" << endl;
  os <<"56" << "\t" << "Ba" << "\t" << "129.906282" << endl;
  os <<"57" << "\t" << "La" << "\t" << "137.90711"  << endl;
  os <<"58" << "\t" << "Ce" << "\t" << "135.907140" << endl;
  os <<"59" << "\t" << "Pr" << "\t" << "140.907647" << endl;
  os <<"60" << "\t" << "Nd" << "\t" << "141.907719" << endl;
  os <<"61" << "\t" << "Pm" << "\t" << "145.0"      << endl;
  os <<"62" << "\t" << "Sm" << "\t" << "143.911998" << endl; 
  os <<"63" << "\t" << "Eu" << "\t" << "150.919847" << endl;
  os <<"64" << "\t" << "Gd" << "\t" << "151.919786" << endl;
  os <<"65" << "\t" << "Tb" << "\t" << "158.925342" << endl;
  os <<"66" << "\t" << "Dy" << "\t" << "155.925277" << endl;
  os <<"67" << "\t" << "Ho" << "\t" << "164.930319" << endl;
  os <<"68" << "\t" << "Er" << "\t" << "161.928775" << endl;
  os << "69" << "\t" << "Tm" << "\t" << "168.934212" << endl;
  os << "70" << "\t" << "Yb" << "\t" << "167.933894" << endl;
  os << "71" << "\t" << "Lu" << "\t" << "174.940770" << endl;
  os << "72" << "\t" << "Hf" << "\t" << "173.940044" << endl;
  os << "73" << "\t" << "Ta" << "\t" << "179.947462" << endl;
  os << "74" << "\t" << "W"  << "\t" << "179.946701" << endl;
  os << "75" << "\t" << "Re" << "\t" << "184.952951" << endl;
  os << "76" << "\t" << "Os" << "\t" << "183.952488" << endl;
  os << "77" << "\t" << "Ir" << "\t" << "190.960584" << endl;
  os << "78" << "\t" << "Pt" << "\t" << "189.959917" << endl;
  os << "79" << "\t" << "Au" << "\t" << "196.966543" << endl;
  os << "80" << "\t" << "Hg" << "\t" << "195.965807" << endl;
  os << "81" << "\t" << "Tl" << "\t" << "202.972320" << endl;
  os << "82" << "\t" << "Pb" << "\t" << "203.973020" << endl;
  os << "83" << "\t" << "Bi" << "\t" << "208.980374" << endl;
  os << "84" << "\t" << "Po" << "\t" << "209.0" << endl;
  os << "85" << "\t" << "At" << "\t" << "210.0" << endl;
  os << "86" << "\t" << "Rn" << "\t" << "222.0" << endl;
  os << "87" << "\t" << "Fr" << "\t" << "223.0" << endl;
  os << "88" << "\t" << "Ra" << "\t" << "226.025" << endl;
  os << "89" << "\t" << "Ac" << "\t" << "227.028" << endl;
  os << "90" << "\t" << "Th" << "\t" << "232.038054" << endl;
  os << "91" << "\t" << "Pa" << "\t" << "231.0359" << endl;
  os << "92" << "\t" << "U"  << "\t" << "234.040946" << endl;
  os << "93" << "\t" << "Np" << "\t" << "237.048" << endl;
  os << "94" << "\t" << "Pu" << "\t" << "244.0" << endl;
  os << "95" << "\t" << "Am" << "\t" << "243.0" << endl;
  os << "96" << "\t" << "Cm" << "\t" << "247.0" << endl;
  os << "97" << "\t" << "Bk" << "\t" << "247.0" << endl;
  os << "98" << "\t" << "Cf" << "\t" << "251.0" << endl;
  os << "99" << "\t" << "Es" << "\t" << "252.0" << endl;
  os <<"100" << "\t" << "Fm" << "\t" << "257.0" << endl;
  os <<"101" << "\t" << "Md" << "\t" << "258.0" << endl;
  os <<"102" << "\t" << "No" << "\t" << "259.0" << endl;
  os <<"103" << "\t" << "Lr" << "\t" << "260.0" << endl;
  os <<"104" << "\t" << "Hx" << "\t" << "1.0078246" << endl;
  os <<"105" << "\t" << "Cx" << "\t" << "12.0000000" << endl;
  os <<"106" << "\t" << "Nx" << "\t" << "14.0030732" << endl;
  os <<"107" << "\t" << "Ox" << "\t" << "15.9949141" << endl;
  os <<"108" << "\t" << "Sx" << "\t" << "31.972070" << endl;  

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
