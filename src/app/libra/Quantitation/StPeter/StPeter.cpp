#include "StPeter.h"
#include "PepXMLParser.h"
#include "ProtXMLParser.h"
#include "Parsers/mzParser/mzParser.h"

#define VERSION "1.6.0"
#define BDATE "March 8 2024"

using namespace mzParser;

//Global Functions
bool    cmdProcess      (char* argv[], int start, int stop, double& dFDR, double& dSL, double& dTol, double& dMinProb, bool& bInt, bool& bDeg, bool& bRun, bool& bExp);
bool    exportToProtXML (char* fn, stpGroup& group);
bool    exportToProtXML (char* fn, vector<stpGroup>& group, vector<string>& files, bool bExp = false);
int     findPeak        (BasicSpectrum& s, double mass, double tol);
int     findPSM         (string& s, vector<stpPSMIndex>& v);
double  matchIons       (BasicSpectrum& s, stpPeptide& p, double tol);
bool    readDB          (char* fn, vector<stpDB>& v);
int     readPepXML      (string fn, double dFDR, double minProb, stpGroup& group, vector<string>& ms2Files, vector<string>& xLabels, map<string,bool>& xFiles );
bool    readProtXML     (char* fn, double dFDR, double minProb, bool bDeg, stpGroup& group, string& pepXML);
bool    readSpectra     (double tol, stpGroup& group, vector<string>& ms2Files, string& outFile, double micrograms, bool bInt, int fID=-1, string* xLbl=NULL, map<string,bool>* xFiles=NULL);
bool    comparePSMIndex (const stpPSMIndex& a, const stpPSMIndex& b);

//Global Variables
double aaMass[128];
double protProb; 
char cDeg[64];
char cFDR[64];
char cSL[64];
char cRun[64];
char cExp[64];
char cTol[64];
char cMinProb[64];

int main(int argc, char* argv[]){

  time_t timeNow;
  int    ret;

  //Export identification
  cout << "\nStPeter: MS2 intensity-based label-free quantification" << endl;
  cout << "Copyright 2015-2024, Jason Winget, Michael Hoopmann, David Shteynberg Institute for Systems Biology" << endl;
  cout << "Version " << VERSION << " " << BDATE << endl;

  if(argc==1){
    cout << "\nUsage: StPeter [options] <protXML>" << endl;
    cout << "\tprotXML = a protXML file containing identifications, scored with iProphet or PeptideProphet" << endl;
    cout << "\tOptions:" << endl;
    cout << "\t\t-d, --degen = allow degenerate peptides in protein quantitation. Default is off." << endl;
    cout << "\t\t-f, --fdr <value> = an FDR cutoff value, e.g. 0.01. Default is 0.01" << endl;
    cout << "\t\t-p, --probabilityMin <value> = minimum probability to use, regardless of FDR cutoff. Default is 0." << endl;
    cout << "\t\t-r, --run = separate protein quantities by run. Cannot be used together with -x. Default is off." << endl;
    cout << "\t\t-x, --exp = separate protein quantities by experiment_label. Cannot be used together with -r. Default is off." << endl;
    cout << "\t\t-s, --sampleLoad <value> = the amount of protein sample measured, in micrograms, e.g. 0.5. Default is 0." << endl;
    cout << "\t\t-t, --tolerance <value> = mass tolerance for matching MS2 peaks (Daltons). Default is 0.4" << endl;
    
    return 1;
  }

  //Set global amino acid masses
  for(int i=0;i<128;i++) aaMass[i]=0;
  aaMass['A']=71.0371103;
  aaMass['C']=103.0091803;
  aaMass['D']=115.0269385;
  aaMass['E']=129.0425877;
  aaMass['F']=147.0684087;
  aaMass['G']=57.0214611;
  aaMass['H']=137.0589059;
  aaMass['I']=113.0840579;
  aaMass['K']=128.0949557;
  aaMass['L']=113.0840579;
  aaMass['M']=131.0404787;
  aaMass['N']=114.0429222;
  aaMass['P']=97.0527595;
  aaMass['Q']=128.0585714;
  aaMass['R']=156.1011021;
  aaMass['S']=87.0320244;
  aaMass['T']=101.0476736;
  aaMass['V']=99.0684087;
  aaMass['W']=186.0793065;
  aaMass['Y']=163.0633228;

  //Data structures needed
  vector<string> ms2Files;
  vector<string> xLabels;
  map<string,bool> xFiles;
 //vector<stpProtein> proteins;
  stpGroup proteins;
  string outFile;
  string pepXMLFile;
  double dFDR;
  double dSL;
  double dTol;
  double dMinProb;
  bool bInt;
  bool bDeg;
  bool bRun;
  bool bExp;
  
  cout << "****** BEGIN StPeter ANALYSIS ******" << endl;
  
  //Timestamp
  time(&timeNow);
  cout << " Time at start of analysis: " << ctime(&timeNow) << endl;

  //Process command-line
  if(!cmdProcess(argv,1,argc-1,dFDR,dSL,dTol,dMinProb,bInt,bDeg,bRun,bExp)) {
    cout << "ERROR: Could not process command-line arguments." << endl;
    return -1;
  }

  //Create output file name - note this may be temporary
  outFile="null";
  //outFile=argv[argc-1];
  //if(bInt) outFile.replace(outFile.size()-9,9,"_intensities.csv");
  //else outFile.replace(outFile.size()-9,9,"_nsi.csv");

  //Read proteins and peptides from protXML at user-defined FDR.
  if(!readProtXML(argv[argc-1],dFDR,dMinProb,bDeg,proteins,pepXMLFile)) {
    cout << "ERROR: Could not read protXML file: " << argv[argc-1] << endl;
    return -2;
  }
  if(proteins.proteins->size()==0) {
    cout << "No proteins in protXML fit criteria for StPeter quantitation. Please address any warning messages for possible solution." << endl;
    return -2;
  }
  
  //Extract scan numbers for the observed peptides from the pepXML files.
  ret=readPepXML(pepXMLFile,dFDR,dMinProb,proteins,ms2Files,xLabels,xFiles);
  switch(ret){
  case 1:
    cout << "ERROR: Could not read pepXML file: " << &pepXMLFile[0] << endl;
    return -3;
  case 2:
    cout << "ERROR: Could not read data files referenced in: " << &pepXMLFile[0] << endl;
    cout << "Please update paths (preferred) or move missing data files to current working directory." << endl;
    return -3;
  default:
    break;
  }

  if(bExp){
    vector<stpGroup> prots;
    for(size_t x=0;x<xLabels.size();x++){
      prots.push_back(proteins);
      vector<string> file;
      if (!readSpectra(dTol, prots[x], ms2Files, outFile, dSL, bInt, -1, &xLabels[x], &xFiles) ) {
        cout << "ERROR: Could not read spectra." << endl;
        return -4;
      }
    }

    //Export back to protXML
    if (!exportToProtXML(argv[argc - 1], prots,xLabels,true)) {
      cout << "ERROR: Could not export SIns to protXML file: " << argv[argc-1] << endl;
      return -5;
    }

  } else if(bRun){
    vector<stpGroup> prots;
    for(size_t a=0;a<ms2Files.size();a++){
      prots.push_back(proteins);
      vector<string> file;
      file.push_back(ms2Files[a]);
      if (!readSpectra(dTol, prots[a], file, outFile, dSL, bInt, (int)a)) {
        cout << "ERROR: Could not read spectra." << endl;
        return -4;
      }
    }

    //Export back to protXML
    if (!exportToProtXML(argv[argc - 1], prots,ms2Files)) {
      cout << "ERROR: Could not export SIns to protXML file: " << argv[argc-1] << endl;
      return -5;
    }

  } else {

    //Extract spectra, match fragment ions, and calculate SIn.
      if(!readSpectra(dTol,proteins,ms2Files,outFile,dSL,bInt)) {
      cout << "ERROR: Could not read spectra." << endl;
      return -4;
    }

    //Export back to protXML
    if(!exportToProtXML(argv[argc-1],proteins)){
      cout << "ERROR: Could not export SIns to protXML file: " << argv[argc - 1] << endl;
      return -5;
    }
  }

  //Timestamp
  time(&timeNow);
  cout << "\n Time at end of analysis: " << ctime(&timeNow) << endl;
 
  cout << "****** FINISHED StPeter ANALYSIS ******" << endl;

  return 0;

}

bool cmdProcess(char* argv[], int start, int stop, double& dFDR, double& dSL, double& dTol, double& dMinProb, bool& bInt, bool& bDeg, bool& bRun, bool& bExp){
  //set defaults
  dFDR=0.01;
  dSL=0;
  dTol=0.4;
  dMinProb=0;
  bInt=false;
  bDeg=false;
  bRun=false;
  bExp=false;
  strcpy(cFDR,"0.01");
  strcpy(cSL,"0");
  strcpy(cTol,"0.4");
  strcpy(cMinProb, "0");
  strcpy(cDeg, "no");
  strcpy(cRun, "no");
  strcpy(cExp, "no");

  //Parse options
  for(int i=start;i<stop;i+=2){
    if(strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--fdr")==0) {
      strcpy(cFDR,argv[i+1]);
      dFDR=atof(argv[i+1]);
	  } else if(strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--sampleLoad")==0) {
	    strcpy(cSL,argv[i+1]);
	    dSL=atof(argv[i+1]);
    } else if(strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--tolerance")==0) {
      strcpy(cTol,argv[i+1]);
      dTol=atof(argv[i+1]);
    } else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--probabilityMin") == 0) {
      strcpy(cMinProb, argv[i + 1]);
      dMinProb = atof(argv[i + 1]);
    } else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--degen") == 0) {
      bDeg = true;
      strcpy(cDeg, "yes");
      i--;
    } else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--run") == 0) {
      bRun= true;
      strcpy(cRun, "yes");
      i--;
    } else if (strcmp(argv[i], "-x") == 0 || strcmp(argv[i], "--exp") == 0) {
      bExp= true;
      strcpy(cExp, "yes");
      i--;
    } else {
      cout << "ERROR: unknown option: " << argv[i] << endl;
      return false;
    }
  }

  if (bExp && bRun) {
    cout << "[ERROR:] separate by run (-r,--run) and separate by experiment options (-x,--exp) should not BOTH be enabled " << endl;
    return false;
  }
  //Report parameters for user summary
  cout << " Parameters:" << endl;
  cout << "  degenerate peptides = " << cDeg << endl;
  cout << "  fdr = " << dFDR << endl;
  cout << "  minimum probability = " << dMinProb << endl;
  cout << "  separate by run = " << cRun << endl;
  cout << "  separate by experiment = " << cExp << endl;
  cout << "  sample load = " << dSL << endl;
  cout << "  tolerance = " << dTol << endl;
  //cout << "  intensities = ";
  //if (bInt) cout << "yes" << endl;
  //else cout << "no" << endl;
  
  return true;
}

bool exportToProtXML(char* fn, stpGroup& group){
  char str[120000];
  char* tok;
  string s;
  size_t i, j, k,n;
  vector<string> buf;
  FILE* f;
  time_t timeNow;

  cout << "\n Exporting SIn to: " << fn << endl;

  //Read file to buffer
  f = fopen(fn, "rt");
  if (f == NULL) {
    cout << "ERROR: could not open: " << fn << endl;
    return false;
  }
  while (!feof(f)){
    if (fgets(str, 120000, f) == NULL) continue;
    tok = strtok(str, "\n\r"); //strip line endings
    if (tok == NULL) continue; //and skip empty lines
    s = tok;
    buf.push_back(s);
  }
  fclose(f);

  //Write file, while integrating StPeter (and removing any prior StPeter results)
  f = fopen(fn, "wt");
  if (f == NULL) {
    cout << "ERROR: could not open: " << fn << endl;
    return false;
  }
  for (i = 0; i<buf.size(); i++){
    //skip previous annotations
    if (buf[i].find("<analysis_result analysis=\"stpeter\">") == 0) {
      i++;
      while (true){
        if (buf[i].find("</analysis_result>") == 0) break;
        i++;
      }
      continue;
    }
    if (buf[i].find("<analysis_summary analysis=\"stpeter\"") != string::npos) {
      i++;
      while (true){
        if (buf[i].find("</analysis_summary>") == 0) break;
        i++;
      }
      continue;
    }

    //add new annotations
    if (buf[i].find("</protein_summary_header>") == 0){
      char timebuf[80];
      time(&timeNow);
      strftime(timebuf, 80, "%Y-%m-%dT%H:%M:%S", localtime(&timeNow));
      fprintf(f, "%s\n", &buf[i][0]);
      fprintf(f, "<analysis_summary analysis=\"stpeter\" time=\"%s\" id=\"1\">\n", timebuf);
      fprintf(f, "<StPeter_analysis_summary version=\"%s\" probability=\"%.4lf\" minimum_probability=\"%s\" degenerate_peptides=\"%s\" FDR=\"%s\" separate_by_run=\"%s\" sampleLoad=\"%s\" tolerance=\"%s\"/>\n", VERSION, protProb, cMinProb, cDeg, cFDR, cRun,cSL, cTol);
      fprintf(f, "</analysis_summary>\n");
      continue;
    }
    if (buf[i].find("<protein protein_name=") != string::npos) {
      fprintf(f, "%s\n", &buf[i][0]);
      //extract the protein name
      j = buf[i].find('"');
      k = buf[i].find('"', j + 1);
      s = buf[i].substr(j + 1, k - j - 1);
      //see if the protein is in our list
      for(k=0;k<group.proteins->size();k++){
        if (group.proteins->at(k).name.compare(s) == 0) break;
      }
      if(k==group.proteins->size()) continue;
      if(group.proteins->at(k).SIN == 0) continue;

      //export SIn
      fprintf(f, "<analysis_result analysis=\"stpeter\">\n");
      if (cDeg[0] == 'n') fprintf(f, "<StPeterQuant SI=\"%.10lf\"", group.proteins->at(k).SI);
      else fprintf(f, "<StPeterQuant dSI=\"%.10lf\"", group.proteins->at(k).SI);
      if(cDeg[0]=='n') fprintf(f, " SIn=\"%.10lf\"", group.proteins->at(k).logSIN);
      else fprintf(f, " dSIn=\"%.10lf\"", group.proteins->at(k).logSIN);
      if (cDeg[0] == 'n') fprintf(f, " counts=\"%d\"", (int)group.proteins->at(k).SC);
      else fprintf(f, " dCounts=\"%.4lf\"", group.proteins->at(k).SC);
      if (cDeg[0] == 'n') fprintf(f, " NSAF=\"%.10lf\"", group.proteins->at(k).logNSAF);
      else fprintf(f, " dNSAF=\"%.10lf\"", group.proteins->at(k).logNSAF);
      fprintf(f, " ng=\"%.10lf\"", group.proteins->at(k).ng);
      fprintf(f, " ngC=\"%.10lf\">\n", group.proteins->at(k).ngC);
      for (n = 0; n < group.proteins->at(k).peptides.size(); n++){
        if (group.proteins->at(k).peptides[n].dSI == 0) continue;
        fprintf(f, "<StPeterQuant_peptide sequence=\"%s\" charge=\"%d\" ", &group.peptides->at(group.proteins->at(k).peptides[n].index).sequence[0], group.peptides->at(group.proteins->at(k).peptides[n].index).charge);
        if (cDeg[0] == 'n') fprintf(f, "SI=\"%.1lf\" pSIn=\"%.6lf\" SC=\"%d\"/>\n", group.proteins->at(k).peptides[n].dSI, group.proteins->at(k).peptides[n].dSIn, (int)group.proteins->at(k).peptides[n].dSC);
        else fprintf(f, "dSI=\"%.1lf\" dpSIn=\"%.6lf\" dSC=\"%.4lf\"/>\n", group.proteins->at(k).peptides[n].dSI, group.proteins->at(k).peptides[n].dSIn, group.proteins->at(k).peptides[n].dSC);
      }
      fprintf(f, "</StPeterQuant>\n");
      fprintf(f, "</analysis_result>\n");
      continue;
    }

    //simply export the line
    fprintf(f, "%s\n", &buf[i][0]);
  }
  fclose(f);
  cout << "  Export successful." << endl;

  /* diagnostics
  f=fopen("dummy.txt","wt");
  fprintf(f,"Protein\tSpecCounts\tSpecIndex\tlog2NSAF\tlog2SIn\tngSC\tngSI\n");
  for(k=0;k<group.proteins->size();k++){
    fprintf(f,"%s",&group.proteins->at(k).name[0]);
    fprintf(f, "\t%.4lf", group.proteins->at(k).SC);
    fprintf(f, "\t%.4lf", group.proteins->at(k).SI);
    fprintf(f, "\t%.8lf", group.proteins->at(k).logNSAF);
    fprintf(f, "\t%.8lf", group.proteins->at(k).logSIN);
    fprintf(f, "\t%.8lf", group.proteins->at(k).ngC);
    fprintf(f, "\t%.8lf", group.proteins->at(k).ng);
    fprintf(f,"\n");
  }
  fclose(f);
  */

  return true;
}

bool exportToProtXML(char* fn, vector<stpGroup>& group, vector<string>& files, bool bExp) {
  char str[120000];
  char* tok;
  string s;
  size_t i, j, k, n;
  vector<string> buf;
  FILE* f;
  time_t timeNow;

  cout << "\n Exporting SIn to: " << fn << endl;

  //Read file to buffer
  f = fopen(fn, "rt");
  if (f == NULL) {
    cout << "ERROR: could not open: " << fn << endl;
    return false;
  }
  while (!feof(f)) {
    if (fgets(str, 120000, f) == NULL) continue;
    tok = strtok(str, "\n\r"); //strip line endings
    if (tok == NULL) continue; //and skip empty lines
    s = tok;
    buf.push_back(s);
  }
  fclose(f);

  //Write file, while integrating StPeter (and removing any prior StPeter results)
  f = fopen(fn, "wt");
  if (f == NULL) {
    cout << "ERROR: could not open: " << fn << endl;
    return false;
  }
  for (i = 0; i < buf.size(); i++) {
    //skip previous annotations
    if (buf[i].find("<analysis_result analysis=\"stpeter\">") == 0) {
      i++;
      while (true) {
        if (buf[i].find("</analysis_result>") == 0) break;
        i++;
      }
      continue;
    }
    if (buf[i].find("<analysis_summary analysis=\"stpeter\"") != string::npos) {
      i++;
      while (true) {
        if (buf[i].find("</analysis_summary>") == 0) break;
        i++;
      }
      continue;
    }

    //add new annotations
    if (buf[i].find("</protein_summary_header>") == 0) {
      char timebuf[80];
      time(&timeNow);
      strftime(timebuf, 80, "%Y-%m-%dT%H:%M:%S", localtime(&timeNow));
      fprintf(f, "%s\n", &buf[i][0]);
      fprintf(f, "<analysis_summary analysis=\"stpeter\" time=\"%s\" id=\"1\">\n", timebuf);
      fprintf(f, "<StPeter_analysis_summary version=\"%s\" probability=\"%.4lf\" minimum_probability=\"%s\" degenerate_peptides=\"%s\" FDR=\"%s\" separate_by_run=\"%s\" sampleLoad=\"%s\" tolerance=\"%s\"/>\n", VERSION, protProb, cMinProb, cDeg, cFDR, cRun, cSL, cTol);
      if (!bExp) {
	for(size_t a=0;a<files.size();a++){
	  fprintf(f,"<StPeter_msrun file=\"%s\" id=\"%d\" />\n", files[a].c_str(),(int)a);
	}
      }
      else {
	for(size_t a=0;a<files.size();a++){
	  fprintf(f,"<StPeter_experiment label=\"%s\" id=\"%d\" />\n", files[a].c_str(),(int)a);
	}
      }
      
      fprintf(f, "</analysis_summary>\n");
      continue;
    }
    if (buf[i].find("<protein protein_name=") != string::npos) {
      fprintf(f, "%s\n", &buf[i][0]);
      //extract the protein name
      j = buf[i].find('"');
      k = buf[i].find('"', j + 1);
      s = buf[i].substr(j + 1, k - j - 1);

      //see if the protein is in our list
      for(size_t a=0;a<group.size();a++){
        for (k = 0; k < group[a].proteins->size(); k++) {
          if (group[a].proteins->at(k).name.compare(s) == 0) break;
        }
        if (k == group[a].proteins->size()) continue;
        if (group[a].proteins->at(k).SIN == 0) continue;
        break;
      }

      //export SIn
      fprintf(f, "<analysis_result analysis=\"stpeter\">\n");

      for (size_t a = 0; a < group.size(); a++) {
        for (k = 0; k < group[a].proteins->size(); k++) {
          if (group[a].proteins->at(k).name.compare(s) == 0) break;
        }
        if (k == group[a].proteins->size()) continue;
        if (group[a].proteins->at(k).SIN == 0) continue;
    
      
        if (cDeg[0] == 'n')
	   if (!bExp) {
	     fprintf(f, "<StPeterQuant run_id=\"%d\" SI=\"%.10lf\"", (int)a,group[a].proteins->at(k).SI);
	   }
	   else {
	     fprintf(f, "<StPeterQuant exp_id=\"%d\" SI=\"%.10lf\"", (int)a,group[a].proteins->at(k).SI);
	   }
	else fprintf(f, "<StPeterQuant dSI=\"%.10lf\"", group[a].proteins->at(k).SI);
        if (cDeg[0] == 'n') fprintf(f, " SIn=\"%.10lf\"", group[a].proteins->at(k).logSIN);
        else fprintf(f, " dSIn=\"%.10lf\"", group[a].proteins->at(k).logSIN);
        if (cDeg[0] == 'n') fprintf(f, " counts=\"%d\"", (int)group[a].proteins->at(k).SC);
        else fprintf(f, " dCounts=\"%.4lf\"", group[a].proteins->at(k).SC);
        if (cDeg[0] == 'n') fprintf(f, " NSAF=\"%.10lf\"", group[a].proteins->at(k).logNSAF);
        else fprintf(f, " dNSAF=\"%.10lf\"", group[a].proteins->at(k).logNSAF);
        fprintf(f, " ng=\"%.10lf\"", group[a].proteins->at(k).ng);
        fprintf(f, " ngC=\"%.10lf\">\n", group[a].proteins->at(k).ngC);
        for (n = 0; n < group[a].proteins->at(k).peptides.size(); n++) {
          if (group[a].proteins->at(k).peptides[n].dSI == 0) continue;
          fprintf(f, "<StPeterQuant_peptide sequence=\"%s\" charge=\"%d\" ", &group[a].peptides->at(group[a].proteins->at(k).peptides[n].index).sequence[0], group[a].peptides->at(group[a].proteins->at(k).peptides[n].index).charge);
          if (cDeg[0] == 'n') fprintf(f, "SI=\"%.1lf\" pSIn=\"%.6lf\" SC=\"%d\"/>\n", group[a].proteins->at(k).peptides[n].dSI, group[a].proteins->at(k).peptides[n].dSIn, (int)group[a].proteins->at(k).peptides[n].dSC);
          else fprintf(f, "dSI=\"%.1lf\" dpSIn=\"%.6lf\" dSC=\"%.4lf\"/>\n", group[a].proteins->at(k).peptides[n].dSI, group[a].proteins->at(k).peptides[n].dSIn, group[a].proteins->at(k).peptides[n].dSC);
        }
        fprintf(f, "</StPeterQuant>\n");

      }
      fprintf(f, "</analysis_result>\n");
      continue;
    }

    //simply export the line
    fprintf(f, "%s\n", &buf[i][0]);
  }
  fclose(f);
  cout << "  Export successful." << endl;

  return true;
}

//Binary search to match closest peak to requested mass
int findPeak(BasicSpectrum& s, double mass, double tol){
  int sz=(int)s.size();
  int lower=0;
  int mid=sz/2;
  int upper=sz;

  double min=mass-tol;
  double max=mass+tol;

  //binary search to closest mass of lowest tolerance
  while(s[mid].mz!=min){
		if(lower>=upper) break;
    if(min<s[mid].mz){
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
		if(mid==sz) {
			mid--;
			break;
		}
	}

  //return tallest peak within mass tolerance
  double maxIntensity=0;
  int j=-1;
  for(int i=mid;i<sz;i++){
    if(s[i].mz<min) continue;
    if(s[i].mz>max) break;
    if(s[i].intensity>maxIntensity){
      maxIntensity=s[i].intensity;
      j=i;
    }
  }
  return j;
}

//Binary search to obtain index of requested PSM
int findPSM(string& s, vector<stpPSMIndex>& v){
  size_t sz=v.size();
  size_t lower=0;
  size_t mid=sz/2;
  size_t upper=sz;
  int i;

  //binary search to string match
  i=v[mid].id.compare(s);
  while(i!=0){
		if(lower>=upper) return -1;
    if(i>0){
      if(mid==0) return -1;
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
		if(mid==sz) return -1;
    i=v[mid].id.compare(s);
	}
  return (int)mid;
}

double matchIons(BasicSpectrum& s, stpPeptide& p, double tol){
  size_t i;
  int j;
  int z;
  double b;
  double y;
  double sum;
  string mod;
  vector<double> ions;
  //bool bKill = false;

  if (p.charge>2) z = 2;
  else z = 1;

  //calculate ion series (b and y)
  b = 0;
  for (i = 0; i<p.preciseSequence.size() - 1; i++){

    //special case for n-terminal mods...add them to first amino acid
    if(i==0 && p.preciseSequence[i]=='n'){
      mod = "";
      i += 2;
      while (p.preciseSequence[i] != ']') mod += p.preciseSequence[i++];
      b += atof(&mod[0]);
      continue;
    }
      
    b += aaMass[p.preciseSequence[i]];
    if (i<p.preciseSequence.size() - 1 && p.preciseSequence[i + 1] == '['){
      //bKill = true;
      mod = "";
      i += 2;
      while (p.preciseSequence[i] != ']') mod += p.preciseSequence[i++];
      b += atof(&mod[0]);
    }
    y = p.mass - b;
    for (j = 1; j <= z; j++){
      ions.push_back((b + j*1.007276466) / j);
      ions.push_back((y + j*1.007276466) / j);
    }
  }

  sum = 0;
  for (i = 0; i<ions.size(); i++){
    j = findPeak(s, ions[i], tol);
    if (j<0) continue;
    if (s[j].intensity<0) continue;
    sum += log(s[j].intensity);  //testing log!!
    s[j].intensity = -s[j].intensity; //don't double count peaks
  }

  return sum;

}

//Reads in a FASTA database, storing each entry as a paired set of name and sequence strings.
bool readDB(char* fn, vector<stpDB>& v){
  stpDB d;
  char  str[10240];
  char  str2[10240];
  char* tok;

  v.clear();
  FILE* f=fopen(fn,"rt");
  if(f==NULL) return false;

  d.name="NIL";
  while(!feof(f)){
    if(fgets(str2,10240,f)==NULL) continue;
    if(strlen(str2)>0){
      tok=strtok(str2,"\r\n");
      if(tok==NULL) continue;
      strcpy(str,tok);
    } else {
      continue;
    }
    if(str[0]=='>') {
      if(d.name.compare("NIL")!=0) v.push_back(d);
      d.name=str;
      d.sequence="";
    } else {
      d.sequence+=str;
    }
  }
  fclose(f);
  v.push_back(d);
  return true;

}

int readPepXML(string fn, double dFDR, double minProb, stpGroup& group, vector<string>& ms2Files, vector<string>& xLabels, map<string,bool>& xFiles){
  size_t j;
  int k, n, z;
  bool bIpro, bLoad;
  bool bWarn = false;
  double prob;
  stpMS2 m;
  string str;
  string str2;
  char cStr[256];
  FILE* f;

  stpPSMIndex psm;
  vector<stpPSMIndex> vPSM;

  map<string,bool> mLabels;
  m.sumIntensity = 0;

  //Read pepXML
  cout << "\n Extracting peptides from pepXML: " << fn << endl;
  PepXMLParser p;
  bLoad = false;
  if (p.parse(&fn[0]) == 0) bLoad = true;
  if (!bLoad) {
    str = fn.substr(fn.rfind('/') + 1, fn.size());
    if (p.parse(&str[0]) == 0) bLoad = true;
  }
  if (!bLoad) {
    str = fn.substr(fn.rfind('\\') + 1, fn.size());
    if (p.parse(&str[0]) == 0) bLoad = true;
  }
  if (!bLoad) return 1;
  prob = p.getProbability(dFDR);
  bIpro = p.getIprophet();
  if (bIpro) cout << "  Minimum PSM probability at " << dFDR << " FDR (iProphet): " << prob << endl;
  else cout << "  Minimum PSM probability at " << dFDR << " FDR (PeptideProphet): " << prob << endl;
  if(prob<minProb){
    cout << "  Probability at " << dFDR << " below user-defined threshold. Raising minimum PSM probability to " << minProb << endl;
    prob=minProb;
  }

  //Get list of contributing spectral files
  ms2Files.clear();
  for (k = 0; k<p.getFileCount(); k++) {
    p.getFileFromList(k, cStr);

    //MH: check if file can be opened. This is redundant with code above, but I'm pressed for time
    //and will put an elegant solution together later. Don't judge me.
    bLoad = false;
    f = fopen(cStr, "rt");
    if (f != NULL) {
      str = cStr;
      bLoad = true;
      fclose(f);
    }
    if (!bLoad && !bWarn){
      cout << "  WARNING: Data file(s) in pepXML not found. Please update paths in pepXML. Attempting to open from current working directory." << endl;
      bWarn = true;
    }
    if (!bLoad) { //strip linux path
      str2 = cStr;
      str = str2.substr(str2.rfind('/') + 1, str2.size());
      f = fopen(&str[0], "rt");
      if (f != NULL) {
        bLoad = true;
        fclose(f);
      }
    }
    if (!bLoad) { //strip windows path
      str2 = cStr;
      str = str2.substr(str2.rfind('\\') + 1, str2.size());
      f = fopen(&str[0], "rt");
      if (f != NULL) {
        bLoad = true;
        fclose(f);
      }
    }
    if (!bLoad) return 2;
    ms2Files.push_back(str);
  }

  //Create a sorted list of all PSMs for faster indexing
  for (j = 0; j<group.peptides->size(); j++){
    psm.id = group.peptides->at(j).sequence;
    sprintf(cStr, "%d", group.peptides->at(j).charge);
    psm.id += cStr;
    psm.pepID = j;
    vPSM.push_back(psm);
  }
  sort(vPSM.begin(), vPSM.end(), comparePSMIndex);

  //Read through all PSMs; keep those above probability cutoff; map them to all instances in list of proteins
  z = 0;
  for (k = 0; k<p.size(); k++){
    if (bIpro){
      if (p[k].iProphetProbability<prob) continue;
    } else {
      if (p[k].probability<prob) continue;
    }

    m.fileID = p[k].fileID;
    m.xLabel = p[k].xLabel;
    m.scanNum = p[k].scanNum;

    str = m.xLabel;
    sprintf(cStr, "_%d", p[k].fileID);
    str += cStr;

    xFiles[str] = true;

    if (mLabels.find(m.xLabel)==mLabels.end()) {
      xLabels.push_back(m.xLabel);
      mLabels[m.xLabel] = true;
    }
    
    //Match PSM to all peptide occurrences in protein list. This is not done optimally here.
    str = p[k].modifiedPeptidePlus;
    sprintf(cStr, "%d", p[k].charge);
    str += cStr;
    n = findPSM(str, vPSM);
    if (n>-1) {
      z++;
      group.peptides->at(vPSM[n].pepID).scans->push_back(m);
      if(group.peptides->at(vPSM[n].pepID).preciseSequence.size()<2){
        group.peptides->at(vPSM[n].pepID).preciseSequence = p[k].modifiedPeptide; //first instance gets the mod masses
      }
      //cout << str << "\t" << group.peptides->at(vPSM[n].pepID).sequence << "\t" << group.peptides->at(vPSM[n].pepID).preciseSequence << endl;
    } 
  }
  cout << "  PSMs above FDR mapped back to proteins: " << z << endl;

  //Sanity check: see which peptides in the protein list have no PSM.
  //k = 0;
  //for (i = 0; i<prots.size(); i++){
  //  for (j = 0; j<prots[i].peptides->size(); j++){
  //    if (prots[i].peptides->at(j).scans->size() == 0) {
  //      k++;
        //cout << &prots[i].name[0] << "\t" << &prots[i].peptides->at(j).sequence[0] << "\t+" << prots[i].peptides->at(j).charge << " has no PSMs" << endl;
  //    }
  //  }
  //}
  //cout << "  " << k << " peptides from protXML were excluded due to PSM FDR below threshold." << endl;
  return 0;
}

bool readProtXML(char* fn, double dFDR, double minProb, bool bDeg, stpGroup& group, string& pepXML) {
  bool bExtractLength = false;
  bool bDB = false;
  double prob;
  int i;
  size_t j, k, n, np, z;
  string db, str;
  stpProtein prot;
  stpPeptide pep;
  vector<stpDB> vDB;

  //Clear any existing data
  prot.SIN = 0;

  cout << "\n Reading protXML: " << fn << endl;

  //Read protXML
  ProtXMLParser p;
  if (!p.parse(fn)) return false;
  prob = p.getProbability(dFDR);
  protProb = prob;
  pepXML = p.getPepXML();

  cout << "  Probability at " << dFDR << " FDR: " << prob << endl;
  if(prob<minProb){
    cout << "  Probability at " << dFDR << " below user-defined threshold. Raising minimum protein probability to " << minProb << endl;
    prob=minProb;
  }

  //Iterate over all proteins in protXML
  n = 0;
  np=0;
  z = 0;
  for (i = 0; i<p.size(); i++){
    if (p[i].probability<prob) continue;
    z++;

    //check if we are in the same group, if not, then add new group
    //if(groups.size()==0 || groups.back().groupNumber!=p[i].groupID){
    //  group.groupNumber=p[i].groupID;
    //  group.peptides->clear();
    //  group.proteins->clear();
    //  groups.push_back(group);
    //}

    //Get the protein name and add it to the group
    prot.name = p[i].proteinName;
    prot.description = p[i].proteinDescription;
    if (p[i].length == 0) bExtractLength = true;
    prot.length = p[i].length;
    group.proteins->push_back(prot);

    //Iterate over peptides in protXML, keeping non-degenerates
    for (j = 0; j<p[i].peptides->size(); j++){
      if (!bDeg && !p[i].peptides->at(j).nonDegenerate) continue;

      //Note that protXML sequences do not match pepXML sequences; pepXML modification masses have more precision.
      pep.sequence = p[i].peptides->at(j).modSequence;
      pep.preciseSequence.clear();
      pep.charge = p[i].peptides->at(j).charge;
      pep.mass = p[i].peptides->at(j).calcNeutralPepMass;
      for(k=0;k<group.peptides->size();k++){
        if(group.peptides->at(k)==pep) break;
      }
      if (k == group.peptides->size()){
        pep.protIndex->clear();
        pep.protIndex->push_back(group.proteins->size()-1);
        group.peptides->push_back(pep);
      } else {
        group.peptides->at(k).protIndex->push_back(group.proteins->size()-1);
      }
      n++;
    }

    //if there is at least one peptide, count the protein as quantifiable
    if (group.peptides->size()>0) np++;
  }

  cout << "  Number of proteins above FDR cutoff: " << z << endl;
  cout << "  Number of quantifiable proteins above FDR cutoff: " << np << endl;
  cout << "  Number of peptides for proteins: " << n << endl;

  //If any protein lengths are zero, extract the real protein length from the database
  if (bExtractLength){
    db = p.getDatabase();
    cout << "  Extracting true protein lengths for zero length proteins. Database: " << &db[0] << endl;
    bDB = readDB(&db[0], vDB);
    if (!bDB){
      str = db.substr(db.rfind('/') + 1, db.size());
      bDB = readDB(&str[0], vDB);
    }
    if (!bDB){
      cout << "  WARNING: Could not read database file: " << &db[0] << " or " << &str[0] << endl;
      cout << "  Zero length proteins will be excluded from analysis." << endl;
      for (z = 0; z<group.proteins->size(); z++){
        if (group.proteins->at(z).length>0) np--;
      }
    } else {
      for(z=0;z<group.proteins->size();z++){
        if (group.proteins->at(z).length>0) continue;
        for (k = 0; k<vDB.size(); k++){
          n = vDB[k].name.find(group.proteins->at(z).name);
          if (n == 1) {
            group.proteins->at(z).length = (int)vDB[k].sequence.size();
            break;
          }
        }
        if (k == vDB.size()){
          cout << "  WARNING: Cannot find protein length in database. Cannot compute SIn for: " << &group.proteins->at(z).name[0] << endl;
          np--;
        }
      }
    }
  }

  cout << "  Number of proteins fitting criteria for StPeter analysis: " << np << endl;
  return true;
}

bool readSpectra(double tol, stpGroup& group, vector<string>& ms2Files, string& outFile, double micrograms, bool bInt, int fID, string* xLbl, map<string,bool>* xFiles){
  size_t i, k, n,z;
  int iTmp, iPercent;
  double score;
  double globalCount=0;
  double globalScore = 0;
  double totalSIn = 0;
  string str;
  char cStr[256];
  bool bLoad;
  double msScore;
  double count;

  BasicSpectrum s;
  MzParser r(&s);

  /*Alternative outputting from a previous version. Leaving out for now.
  //Open file to export summary (independent of ProtXML)
  string sample = outFile.substr(0, outFile.find('.'));
  FILE* f = fopen(&outFile[0], "wt");
  if (bInt) fprintf(f, "ProteinName,PeptideSequence,PrecursorCharge,Run,Sample,Intensity\n");
  else fprintf(f, "Protein Name,Protein Description,SIN\n");
  */

  cout << "\n Extracting peak intensities from " << ms2Files.size() << " spectral files:" << endl;

  //Iterate over all spectra files
  for (i = 0; i<ms2Files.size(); i++){

    if (xLbl != NULL && xFiles != NULL) {
      str = *xLbl;
      sprintf(cStr, "_%ld", i);
      str += cStr;
      if (xFiles->find(str) == xFiles->end() || (*xFiles)[str] != true) {
	continue;
      }

    }


    //Export progress
    iPercent = 0;
    printf("  %zu of %zu: %s   %2d%%", i + 1, ms2Files.size(), &ms2Files[i][0], iPercent);
    fflush(stdout);

    //If file cannot be found, try in current working directory before failing.
    str = ms2Files[i].substr(ms2Files[i].rfind('/') + 1, ms2Files[i].size());
    bLoad = r.load(&ms2Files[i][0]);
    if (!bLoad) bLoad = r.load(&str[0]);
    if (!bLoad) {
      str = ms2Files[i].substr(ms2Files[i].rfind('\\') + 1, ms2Files[i].size());
      bLoad = r.load(&str[0]);
    }
    if (!bLoad) {
      cout << "\n  ERROR reading file: " << &ms2Files[i][0] << " or " << &str[0] << endl;
      return false;
    }


    //Iterate over all peptides
    for (k = 0; k<group.peptides->size(); k++){
        
      //Update progress meter
      iTmp = (int)(k * 100 / group.peptides->size());
      if (iTmp>iPercent){
        iPercent = iTmp;
        printf("\b\b\b%2d%%", iPercent);
        fflush(stdout);
      }

      //Export peptide output, if requested by user
      //if (bInt) fprintf(f, "%s,%s,%d,%s,%s,", &prots[j].name[0], &prots[j].peptides->at(k).sequence[0], prots[j].peptides->at(k).charge, &str[0], &sample[0]);

      //Calculate sum of all PSMs for a peptide
      msScore = 0;
      count = 0;
      for (n = 0; n<group.peptides->at(k).scans->size(); n++){
        if(fID>-1){
          if (group.peptides->at(k).scans->at(n).fileID != fID) continue;
        } else {
          if (group.peptides->at(k).scans->at(n).fileID != i) continue;
        }
	if (xLbl != NULL) {
	  if (group.peptides->at(k).scans->at(n).xLabel != *xLbl) continue;
	}
        if (!r.readSpectrum(group.peptides->at(k).scans->at(n).scanNum))
	  return false;
        score = matchIons(s, group.peptides->at(k), tol);
        //How to deal with bad search/prophet results? If there are no matched ions, they're probably wrong
        if (score < 0.01) {
          //cout << "WARNING: no score for: " << group.peptides->at(k).sequence << "\t" << group.peptides->at(k).preciseSequence << "\t" << group.peptides->at(k).scans->at(n).scanNum << endl;
          //exit(1);
        }
        count++;
        globalCount++;
        msScore += score;
        globalScore += score;
      }

      group.peptides->at(k).spectralCount += count;
      group.peptides->at(k).spectralIndex += msScore;

      //Export peptide output,if requested by user
      //if (bInt) fprintf(f, "%.4lf\n", msScore);
    }

    //Finalize progress meter
    printf("\b\b\b%2d%%\n", 100);
    fflush(stdout);

  }

  string repHdr = "Global";

  if (xLbl) {
    repHdr = "Experiment ";
    repHdr += *xLbl;
  }

  
  
  //Export global SI for user summary
  cout << "  " << repHdr << " Counts: " << globalCount << endl;
  cout << "  " << repHdr << " SI: " << globalScore << endl;

  //Compute weighted SI for each protein
  double* uniProt=NULL;
  double* uniProtC=NULL;
  double siSum,cSum;
  size_t sz;
  stpdPep pepInfo;
  //double reSum=0;

  //reset all sums and memory
  siSum=0;
  sz = group.proteins->size();
  if(uniProt!=NULL) delete [] uniProt;
  uniProt = new double[sz];
  if (uniProtC != NULL) delete[] uniProtC;
  uniProtC = new double[sz];
  for(k=0;k<sz;k++) uniProt[k]=uniProtC[k]=0;

  //compute contribution from unique peptides
  for(k=0;k<group.peptides->size();k++){
    if(group.peptides->at(k).protIndex->size()==1) {
      z = group.peptides->at(k).protIndex->at(0);
      uniProt[z] += group.peptides->at(k).spectralIndex;
      uniProtC[z] += group.peptides->at(k).spectralCount;
      //siSum += groups[j].peptides->at(k).spectralIndex;
    }
  }

  //compute SI and SC
  for (k = 0; k<group.peptides->size(); k++){

    if (group.peptides->at(k).protIndex->size() == 1) { //non-degenerate peptides
      z = group.peptides->at(k).protIndex->at(0);
      group.proteins->at(z).SI+=group.peptides->at(k).spectralIndex;
      group.proteins->at(z).SC+=group.peptides->at(k).spectralCount;
      pepInfo.dSI = group.peptides->at(k).spectralIndex;
      pepInfo.dSC = group.peptides->at(k).spectralCount;
      pepInfo.dSIn = log2(pepInfo.dSI/globalScore);
      pepInfo.index = k;
      group.proteins->at(z).peptides.push_back(pepInfo);

    } else { //degenerate peptides
      siSum=0;
      cSum=0;
      pepInfo.clear();
      for (n = 0; n<group.peptides->at(k).protIndex->size(); n++){
        z = group.peptides->at(k).protIndex->at(n);
        siSum += uniProt[z];
        cSum += uniProtC[z];
      }
      if (siSum == 0) { //if all proteins have no degenerate peptides, just divvy them up equally
        for (n = 0; n<group.peptides->at(k).protIndex->size(); n++){
          z = group.peptides->at(k).protIndex->at(n);
          group.proteins->at(z).SI += group.peptides->at(k).spectralIndex / group.peptides->at(k).protIndex->size();
          group.proteins->at(z).SC += group.peptides->at(k).spectralCount / group.peptides->at(k).protIndex->size();
        }
        continue;
      }
      for (n = 0; n<group.peptides->at(k).protIndex->size();n++){
        z = group.peptides->at(k).protIndex->at(n);
        /* not real code
        group.proteins->at(z).SI += group.peptides->at(k).spectralIndex;
        group.proteins->at(z).SC += group.peptides->at(k).spectralCount;
        pepInfo.dSI = group.peptides->at(k).spectralIndex;
        pepInfo.dSC = group.peptides->at(k).spectralCount;
        pepInfo.index = k;
        group.proteins->at(z).peptides.push_back(pepInfo);
        */
        group.proteins->at(z).SI += group.peptides->at(k).spectralIndex*uniProt[z]/siSum;
        group.proteins->at(z).SC += group.peptides->at(k).spectralCount*uniProtC[z]/cSum;
        pepInfo.dSI = group.peptides->at(k).spectralIndex*uniProt[z] / siSum;
        pepInfo.dSC = group.peptides->at(k).spectralCount*uniProtC[z] / cSum;
        pepInfo.dSIn = log2(pepInfo.dSI / globalScore);
        pepInfo.index = k;
        group.proteins->at(z).peptides.push_back(pepInfo);
      }
    }
  }

  //clean up memory
  if (uniProt != NULL) delete[] uniProt;
  if (uniProtC != NULL) delete[] uniProtC;

  //sanity check
  double reGlobal = 0;
  double reGlobalC = 0;
  for(k=0;k<group.proteins->size();k++){
    reGlobal+=group.proteins->at(k).SI;
    reGlobalC+=group.proteins->at(k).SC;
  }

  cout << "\t" << repHdr << " Weighted Recounts: " << endl;
  cout << "\t\tscores: " << reGlobal << " vs. original: " << globalScore << endl;
  cout << "\t\tcounts: " << reGlobalC << " vs. original: " << globalCount << endl;

  //Normalize SI for each protein
  n = 0;
  for(k=0;k<group.proteins->size();k++){
    if (group.proteins->at(k).SI == 0) continue;
    score = group.proteins->at(k).SI / globalScore / group.proteins->at(k).length;
    totalSIn += score;
    group.proteins->at(k).SIN = score;
    group.proteins->at(k).logSIN = log(score) / log(2.0);
    n++;
    //Export SIn output, if requested by user
    //if (!bInt) fprintf(f, "%s,\"%s\",%.4lf\n", &prots[j].name[0], &prots[j].description[0], prots[j].logSIN);
  }
  cout << "  " << n << " proteins have SIn." << endl;
  cout << "  total SIn = " << totalSIn << endl;
  
  //Convert to nanograms
  for (k = 0; k<group.proteins->size(); k++){
    group.proteins->at(k).ng = group.proteins->at(k).SIN / totalSIn*micrograms * 1000;
    //cout << groups[j].proteins->at(k).ng << endl;
  }

  //Compute NSAF
  n = 0;
  totalSIn=0;
  for (k = 0; k<group.proteins->size(); k++){
    if (group.proteins->at(k).SC == 0) continue;
    group.proteins->at(k).SAF = group.proteins->at(k).SC / group.proteins->at(k).length;
    totalSIn += group.proteins->at(k).SAF;
  }
  for (k = 0; k<group.proteins->size(); k++){
    if (group.proteins->at(k).SC == 0) continue;
    group.proteins->at(k).NSAF = group.proteins->at(k).SAF/totalSIn;
    group.proteins->at(k).logNSAF = log(group.proteins->at(k).NSAF) / log(2.0);
    n++;
  }
  cout << "  " << n << " proteins have NSAF." << endl;
  //cout << "  total NSAF = " << score << endl;

  //Convert to nanograms
  for (k = 0; k<group.proteins->size(); k++){
    group.proteins->at(k).ngC = group.proteins->at(k).NSAF*micrograms * 1000;
    //cout << groups[j].proteins->at(k).ngC << endl;
  }

  //fclose(f);
  return true;

}

bool comparePSMIndex(const stpPSMIndex& a, const stpPSMIndex& b){
  return (a.id.compare(b.id)<0);
}


