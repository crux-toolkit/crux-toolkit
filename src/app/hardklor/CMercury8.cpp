/*=====================================================================*/
/* Program MERCURY2.C                                                  */
/*                                                                     */
/* MERCURY5 is a version of MERCURY2 using (mostly) double percision.  */
/* instead of floating point arithmatic. It gives more accurate        */
/* intensity values than MERCURY2.                                     */
/* MERCURY2 is an integer based version of MERCURY, although most of   */
/* the arithmetic is floating point. Using integer (changed to float)  */
/* values for isotopic masses, the calculation can be performed with   */
/* a much smaller data set and is extremely fast.  The ASCII output    */
/* file is a stick representation of the mass spectrum. There is no    */
/* ultrahigh resolution mode in this program.                          */
/*                                                                     */
/* Algorithm by       Alan L. Rockwood                                 */
/* Programming by     Steve Van Orden                                  */
/*=====================================================================*/

/*=====================================================================*/
/* C++ implementation (CMercury5) by Michael Hoopmann, 2004            */
/*                                                                     */
/* To use:                                                             */
/*   1. Create CMercury5 object.                                       */
/*   2. Call GoMercury(formula, [optional] charge, [optional] filename)*/
/*   3. Optionally call Echo(true) to display output to screen.        */
/*                                                                     */
/* Example:                                                            */
/*   #include "CMercury5.h"                                            */
/*   using namespace std;                                              */
/*   int main(){                                                       */
/*     CMercury5 dist;                                                 */
/*     dist.Echo(true);                                                */
/*     dist.GoMercury("C6H12O6",1);                                    */
/*     return 0;                                                       */
/*   };                                                                */
/*                                                                     */
/*  Including the optional filename to GoMercury outputs the           */
/*  distribution to file. The data can be manipulated in code in the   */
/*  FixedData vector. See the header files CMercury5.h and mercury.h   */
/*  for struct type.                                                   */
/*=====================================================================*/
 
#include "CMercury8.h"
#include <cmath>
#include <iostream>
using namespace std;

CMercury8::CMercury8(){
  InitializeData();
  showOutput = false;
  bAccMass = false;
  bRelAbun = true;
  zeroMass=0;
}

CMercury8::CMercury8(char* fn){
  InitializeData(fn);
  showOutput = false;
  bAccMass = false;
  bRelAbun = true;
  zeroMass=0;
}

CMercury8::~CMercury8(){
  int Z;

  for (Z=0; Z<=MAXAtomNo; Z++) {

    delete [] Element[Z].IsoMass;
    delete [] Element[Z].IntMass;
    delete [] Element[Z].IsoProb;
    //if(Element[Z].WrapMass!=NULL) delete [] Element[Z].WrapMass;

    delete [] Orig[Z].IsoMass;
    delete [] Orig[Z].IntMass;
    delete [] Orig[Z].IsoProb;
    //if(Orig[Z].WrapMass!=NULL) delete [] Orig[Z].WrapMass;
  }
}

void CMercury8::Echo(bool b){
  showOutput = b;
}

//quick hack for N-enrichment
//This needs to be expanded to be universal for all elements
void CMercury8::Enrich(int c,int e,double d){
  int i=0;
  int j=0;
  double ab=0;
	float f=(float)d;

  //c++;
  if(showOutput) cout << "Enrich: " << Element[c].Symbol << " " << c << "\t" << d << endl;

  //Find highest probability
  for(i=0;i<Element[c].NumIsotopes;i++){
    if(Element[c].IsoProb[i]>ab){
      j=i;
      ab=Element[c].IsoProb[i];
    }
  }

  //Normalize all isotope abundances
  for(i=0;i<Element[c].NumIsotopes;i++){
    Element[c].IsoProb[i]/=Element[c].IsoProb[j];
  }

  //Calculate enrichment
  for(i=0;i<Element[c].NumIsotopes;i++){
    if(i==e) Element[c].IsoProb[i]=(1-f)*Element[c].IsoProb[i]+f;
    else Element[c].IsoProb[i]=(1-f)*Element[c].IsoProb[i];
  }

  EnrichAtoms.push_back(c);

}
  

/*************************************************/
/* FUNCTION Intro - called by main()             */
/*************************************************/
void CMercury8::Intro() {
   printf("*********************************************************************\n");
   printf("*                       M E R C U R Y  V I I I                      *\n");
   printf("*                                                                   *\n");
   printf("*  An Integer based Fourier transform isotopic distibution program  *\n");
   printf("*          Now capable of calculating accurate masses!              *\n");
   printf("*********************************************************************\n");
   printf("\n");
   printf("      Algorithm by : Alan L. Rockwood\n\n");
   printf("      Program by   : Steven L. Van Orden - 1\n");
   printf("                     Michael R. Hoopmann - 2\n\n");
   printf("      Developed at : 1 - Pacific Northwest Laboratories / \n");
   printf("                         Battelle Northwest\n");
   printf("                         in the laboratory of Richard D. Smith\n");
   printf("                     2 - University of Washington\n");
   printf("                         Department of Genome Science\n");
   printf("                         Michael J. MacCoss laboratory\n\n\n");
 
}  /* End of Intro() */
 
/***************************************************/
/* FUNCTION InitializeData - called by constructor */
/***************************************************/
//This function reads the ISOTOPE.DAT file that must be in the same
//folder as the application.
void CMercury8::InitializeData(const char* fn) {


  if (fn == NULL) {
    //use hardcoded defaults.
    InitializeDataHardcoded();
  } else {

    FILE *ElementFile;
    int  i, Z;
 
    if ((ElementFile = fopen(fn, "rt")) == NULL) {
      printf("\nError - Cannot open File: ISOTOPE.DAT\n");
      InitializeDataHardcoded();
      return;
    }
     
    for (Z=0; Z<=MAXAtomNo; Z++) {
      Element[Z].Symbol[0]=Element[Z].Symbol[1]=Element[Z].Symbol[2]=0;
  
      fscanf(ElementFile,"%2s %d\n", Element[Z].Symbol,&Element[Z].NumIsotopes);
      strcpy(Orig[Z].Symbol,Element[Z].Symbol);
      Orig[Z].NumIsotopes = Element[Z].NumIsotopes;
  
      Element[Z].IsoMass = new float[Element[Z].NumIsotopes+1];
      Element[Z].IntMass = new int[Element[Z].NumIsotopes+1];
      Element[Z].IsoProb = new float[Element[Z].NumIsotopes+1];
      //Element[Z].WrapMass = NULL;
  
      Orig[Z].IsoMass = new float[Orig[Z].NumIsotopes+1];
      Orig[Z].IntMass = new int[Orig[Z].NumIsotopes+1];
      Orig[Z].IsoProb = new float[Orig[Z].NumIsotopes+1];
      //Orig[Z].WrapMass = NULL;
      
      for (i=0; i<Element[Z].NumIsotopes; i++) {
        fscanf(ElementFile, "%f \n", &Element[Z].IsoMass[i]);
        fscanf(ElementFile, "%f \n", &Element[Z].IsoProb[i]);
        Element[Z].IntMass[i] = (int)(Element[Z].IsoMass[i]+0.5);
        Orig[Z].IsoMass[i]=Element[Z].IsoMass[i];
        Orig[Z].IsoProb[i]=Element[Z].IsoProb[i];
        Orig[Z].IntMass[i]=Element[Z].IntMass[i];
      }
        
      Element[Z].NumAtoms = 0;
      Element[Z].IsoMass[Element[Z].NumIsotopes] = 0;
      Element[Z].IsoProb[Element[Z].NumIsotopes] = 0;
      Orig[Z].NumAtoms = 0;
      Orig[Z].IsoMass[Orig[Z].NumIsotopes] = 0;
      Orig[Z].IsoProb[Orig[Z].NumIsotopes] = 0;
  
      fscanf(ElementFile, " \n");
    }
  
    fclose(ElementFile);
  
  }
}

void CMercury8::InitializeDataHardcoded() {

  cerr << "Loading hardcoded table" << endl;
  for (int Z=0;Z<-MAXAtomNo;Z++) {
    Element[Z].Symbol[0]=Element[Z].Symbol[1]=Element[Z].Symbol[2]=0;
  }

  Element[0].Symbol[0] = 'X';
  Element[0].NumIsotopes = 2;
  Element[0].IsoMass = new float[3];
  Element[0].IsoMass[0] = 1;
  Element[0].IsoMass[1] = 2;
  Element[0].IsoProb = new float[3];
  Element[0].IsoProb[0] = 0.9;
  Element[0].IsoProb[1] = 0.1;

  Element[1].Symbol[0] = 'H';
  Element[1].NumIsotopes = 2;
  Element[1].IsoMass = new float[3];
  Element[1].IsoMass[0] = 1.0078246;
  Element[1].IsoMass[1] = 2.0141021;
  Element[1].IsoProb = new float[3];
  Element[1].IsoProb[0] = 0.999855;
  Element[1].IsoProb[1] = 0.000145;

  Element[2].Symbol[0] = 'H';
  Element[2].Symbol[1] = 'e';
  Element[2].NumIsotopes = 2;
  Element[2].IsoMass = new float[3];
  Element[2].IsoMass[0] = 3.01603;
  Element[2].IsoMass[1] = 4.00260;
  Element[2].IsoProb = new float[3];
  Element[2].IsoProb[0] = 0.00000138;
  Element[2].IsoProb[1] = 0.99999862;

  Element[3].Symbol[0] = 'L';
  Element[3].Symbol[1] = 'i';
  Element[3].NumIsotopes = 2;
  Element[3].IsoMass = new float[3];
  Element[3].IsoMass[0] = 6.015121;
  Element[3].IsoMass[1] = 7.016003;
  Element[3].IsoProb = new float[3];
  Element[3].IsoProb[0] = 0.075;
  Element[3].IsoProb[1] = 0.925;

  Element[4].Symbol[0] = 'B';
  Element[4].Symbol[1] = 'e';
  Element[4].NumIsotopes = 1;
  Element[4].IsoMass = new float[2];
  Element[4].IsoMass[0] = 9.012182;
  Element[4].IsoProb = new float[2];
  Element[4].IsoProb[0] = 1.0;

  Element[5].Symbol[0] = 'B';
  Element[5].NumIsotopes = 2;
  Element[5].IsoMass = new float[3];
  Element[5].IsoMass[0] = 10.012937;
  Element[5].IsoMass[1] = 11.009305;
  Element[5].IsoProb = new float[3];
  Element[5].IsoProb[0] = 0.199;
  Element[5].IsoProb[1] = 0.801;

  Element[6].Symbol[0] = 'C';
  Element[6].NumIsotopes = 2;
  Element[6].IsoMass = new float[3];
  Element[6].IsoMass[0] = 12.0000000;
  Element[6].IsoMass[1] = 13.0033554;
  Element[6].IsoProb = new float[3];
  Element[6].IsoProb[0] = 0.98916;
  Element[6].IsoProb[1] = 0.01084;
  
  Element[7].Symbol[0] = 'N';
  Element[7].NumIsotopes = 2;
  Element[7].IsoMass = new float[3];
  Element[7].IsoMass[0] = 14.0030732;
  Element[7].IsoMass[1] = 15.0001088;
  Element[7].IsoProb = new float[3];
  Element[7].IsoProb[0] = 0.99633;
  Element[7].IsoProb[1] = 0.00366;
  
  Element[8].Symbol[0] = 'O';
  Element[8].NumIsotopes = 3;
  Element[8].IsoMass = new float[4];
  Element[8].IsoMass[0] = 15.9949141;
  Element[8].IsoMass[1] = 16.9991322;
  Element[8].IsoMass[2] = 17.9991616;
  Element[8].IsoProb = new float[4];
  Element[8].IsoProb[0] = 0.997576009706;
  Element[8].IsoProb[1] = 0.000378998479;
  Element[8].IsoProb[2] = 0.002044991815;

  Element[9].Symbol[0] = 'F';
  Element[9].NumIsotopes = 1;
  Element[9].IsoMass = new float[2];
  Element[9].IsoMass[0] = 18,9984032;
  Element[9].IsoProb = new float[2];
  Element[9].IsoProb[0] = 1.0;

  Element[10].Symbol[0] = 'N';
  Element[10].Symbol[1] = 'e';
  Element[10].NumIsotopes = 3;
  Element[10].IsoMass = new float[4];
  Element[10].IsoMass[0] = 19.992435;
  Element[10].IsoMass[1] = 20.993843;
  Element[10].IsoMass[2] = 21.991383;
  Element[10].IsoProb = new float[4];
  Element[10].IsoProb[0] = 0.9048;
  Element[10].IsoProb[1] = 0.0027;
  Element[10].IsoProb[2] = 0.0925;

  Element[11].Symbol[0] = 'N';
  Element[11].Symbol[1] = 'a';
  Element[11].NumIsotopes = 1;
  Element[11].IsoMass = new float[2];
  Element[11].IsoMass[0] = 22.989767;
  Element[11].IsoProb = new float[2];
  Element[11].IsoProb[0] = 1.0;

  Element[12].Symbol[0] = 'M';
  Element[12].Symbol[1] = 'g';
  Element[12].NumIsotopes = 3;
  Element[12].IsoMass = new float[4];
  Element[12].IsoMass[0] = 23.985042;
  Element[12].IsoMass[1] = 24.985837;
  Element[12].IsoMass[2] = 25.982593;
  Element[12].IsoProb = new float[4];
  Element[12].IsoProb[0] = 0.7899;
  Element[12].IsoProb[1] = 0.1000;
  Element[12].IsoProb[2] = 0.1101;
  
  Element[13].Symbol[0] = 'A';
  Element[13].Symbol[1] = 'l';
  Element[13].NumIsotopes = 1;
  Element[13].IsoMass = new float[2];
  Element[13].IsoMass[0] = 26.981539;
  Element[13].IsoProb = new float[2];
  Element[13].IsoProb[0] = 1.0;
  
  Element[14].Symbol[0] = 'S';
  Element[14].Symbol[1] = 'i';
  Element[14].NumIsotopes = 3;
  Element[14].IsoMass = new float[4];
  Element[14].IsoMass[0] = 27.976927;
  Element[14].IsoMass[1] = 28.976495;
  Element[14].IsoMass[2] = 29.973770;
  Element[14].IsoProb = new float[4];
  Element[14].IsoProb[0] = 0.9223;
  Element[14].IsoProb[1] = 0.0467;
  Element[14].IsoProb[2] = 0.0310;

  Element[15].Symbol[0] = 'P';
  Element[15].NumIsotopes = 1;
  Element[15].IsoMass = new float[2];
  Element[15].IsoMass[0] = 30.973762;
  Element[15].IsoProb = new float[2];
  Element[15].IsoProb[0] = 1.0;

  Element[14].Symbol[0] = 'S';
  Element[14].NumIsotopes = 4;
  Element[14].IsoMass = new float[5];
  Element[14].IsoMass[0] = 31.972070;
  Element[14].IsoMass[1] = 32.971456;
  Element[14].IsoMass[2] = 33.967866;
  Element[14].IsoMass[3] = 35.967080;
  Element[14].IsoProb = new float[5];
  Element[14].IsoProb[0] = 0.95021;
  Element[14].IsoProb[1] = 0.00745;
  Element[14].IsoProb[2] = 0.04221;
  Element[14].IsoProb[3] = 0.00013;

  Element[15].Symbol[0] = 'C';
  Element[15].Symbol[1] = 'l';
  Element[15].NumIsotopes = 2;
  Element[15].IsoMass = new float[3];
  Element[15].IsoMass[0] = 34.9688531;
  Element[15].IsoMass[1] = 36.9659034;
  Element[15].IsoProb = new float[3];
  Element[15].IsoProb[0] = 0.755290;
  Element[15].IsoProb[1] = 0.244710;
  
  Element[16].Symbol[0] = 'A';
  Element[16].Symbol[1] = 'r';
  Element[16].NumIsotopes = 3;
  Element[16].IsoMass = new float[4];
  Element[16].IsoMass[0] = 35.967545;
  Element[16].IsoMass[1] = 37.962732;
  Element[16].IsoMass[2] = 39.962384;
  Element[16].IsoProb = new float[4];
  Element[16].IsoProb[0] = 0.00337;
  Element[16].IsoProb[1] = 0.00063;
  Element[16].IsoProb[2] = 0.99600;

  for (int Z=0;Z<-MAXAtomNo;Z++) {
    strcpy(Orig[Z].Symbol,Element[Z].Symbol);
    Orig[Z].NumIsotopes = Element[Z].NumIsotopes;
    for (int i=0; i<Element[Z].NumIsotopes; i++) {
      Element[Z].IntMass[i] = (int)(Element[Z].IsoMass[i]+0.5);
      Orig[Z].IsoMass[i]=Element[Z].IsoMass[i];
      Orig[Z].IsoProb[i]=Element[Z].IsoProb[i];
      Orig[Z].IntMass[i]=Element[Z].IntMass[i];
    }
        
    Element[Z].NumAtoms = 0;
    Element[Z].IsoMass[Element[Z].NumIsotopes] = 0;
    Element[Z].IsoProb[Element[Z].NumIsotopes] = 0;
    Orig[Z].NumAtoms = 0;
    Orig[Z].IsoMass[Orig[Z].NumIsotopes] = 0;
    Orig[Z].IsoProb[Orig[Z].NumIsotopes] = 0;
  }

}


 
/*************************************************/
/* FUNCTION CalcVariances - called by main()     */
/*************************************************/
void CMercury8::CalcVariances(double *MolVar, double *IntMolVar, int NumElements){
  int i, j, Z;
  double Var, IntVar;
  double avemass, intavemass;
  
  *MolVar = *IntMolVar = 0;
  for (i=0; i<NumElements; i++) {
    Z = AtomicNum[i];
    avemass = intavemass = 0;
    for (j=0; j<Element[Z].NumIsotopes; j++){
      avemass += Element[Z].IsoMass[j] * Element[Z].IsoProb[j];
      intavemass += Element[Z].IntMass[j] * Element[Z].IsoProb[j];
    };
    Var = IntVar = 0;
    for (j=0; j<Element[Z].NumIsotopes; j++){
      Var += (Element[Z].IsoMass[j] - avemass) * (Element[Z].IsoMass[j] - avemass) * Element[Z].IsoProb[j];
      IntVar += (Element[Z].IntMass[j] - intavemass) * (Element[Z].IntMass[j] - intavemass) * Element[Z].IsoProb[j];
    };
    *MolVar += Element[Z].NumAtoms * Var;
    *IntMolVar += Element[Z].NumAtoms * IntVar;
  };

	//monoisotopic mass (zero mass, EXACT)
	zeroMass=0;
	for (i=0; i<NumElements; i++) {
    Z = AtomicNum[i];
		zeroMass+=(Element[Z].IsoMass[0] * Element[Z].NumAtoms);
	};
	monoMass=zeroMass;
  
};  /* End of CalcVariances() */

/*************************************************/
/* FUNCTION CalcMassRange - called by main()     */
/*************************************************/
void CMercury8::CalcMassRange(int *MassRange, double MolVar, int charge, int type) {
   int i;
   double dPoints;
 
   //This is insufficient without adding the one to the end
   if ((type == 1) || (charge == 0)) dPoints = (sqrt(1+MolVar)*10);
   else  dPoints = (sqrt(1+MolVar)*10/charge);  /* +/- 5 sd's : Multiply charged */

   /* Set to nearest (upper) power of 2 */
   for (i=1024; i>0; i/=2) {
     if (i < dPoints) {
       *MassRange = i * 2 * 2;   //MRH: Added extra power of 2 since this rule is often insufficient
       i = 0;
     };
   };
   
}  /* End of CalcMassRange() */
 
/*************************************************/
/* FUNCTION AddElement - called by ParseMF()     */
/*************************************************/

//Atom is the atomic abbreviation, Ecount is nth element in the formula
//Acount is the number of atoms of the element in the formula
void CMercury8::AddElement(char Atom[3], int Ecount, int Acount) {

  int Z, FOUND=0;
 
  for (Z=1; Z<=MAXAtomNo; Z++) {

    if (strcmp(Atom,Element[Z].Symbol) == 0) {

      if (Element[Z].NumAtoms != 0) { 

	printf("\nError - the element %s has been entered twice in molecular formula\n", Element[Z].Symbol);
	exit(-1);

      } else {

	AtomicNum[Ecount] = Z;
        //cerr << Atom << " " << AtomicNum[Ecount] << endl;
	Element[Z].NumAtoms = Acount;
	//Element[Z].WrapMass = new int[Element[Z].NumIsotopes+1];
	//Element[Z].WrapMass[Element[Z].NumIsotopes] = 0;
	FOUND=1;
	break;

      };

    };
  };

  if (!FOUND) {
    printf("\nError - Unknown element in Molecular Formula\n");
    exit(-1);
  };

}
 
/*************************************************/
/* FUNCTION ParseMF - called by main()           */
/*************************************************/
//Return codes
//	0: Success
//	-1: Invalid character
int CMercury8::ParseMF(char MF[], int *elementcount) {
  int COND, ERRFLAG;
  int atomcount;
  char Atom[3], errorch ='\0';

  atomcount=0; COND=0; ERRFLAG=0;
  Atom[0] = Atom[1] = Atom[2] = '\0';

  unsigned int pos=0;
  unsigned int peek=0;
  //unsigned int count=0;
  char digit[2];
  bool bFirst=true;
  atomcount=0;
  *elementcount=0;
  while(pos<strlen(MF) && ERRFLAG==0){
    if(isupper(MF[pos])){
      //Add the last atom
      if(!bFirst) {
        AddElement(Atom,(*elementcount)++,atomcount);
        atomcount=0;
      } else {
        bFirst=false;
      }
      
      peek=pos+1;
      if(peek==strlen(MF)){
        //reached end of string, add this single atom
        Atom[0]=MF[pos];
        Atom[1]=Atom[2]='\0';
        atomcount=1;
        pos++;
        continue;
      }
       
      if(isupper(MF[peek])){
        //This is a single atom
        Atom[0]=MF[pos];
        Atom[1]=Atom[2]='\0';
        atomcount=1;
        pos++;
      } else if(islower(MF[peek])){
        //Set the atom name
        Atom[0]=MF[pos];
        Atom[1]=MF[peek];
        Atom[2]='\0';
        atomcount=0;
        pos+=2;
      } else if(isdigit(MF[peek])){
        //Set the atom name
        Atom[0]=MF[pos];
        Atom[1]=Atom[2]='\0';
        atomcount=0;
        pos++;
      } else {
        errorch=MF[peek];
        ERRFLAG=1;
      }
         
    } else if(isdigit(MF[pos])){
      digit[0]=MF[pos];
      digit[1]='\0';
      atomcount*=10;
      atomcount+=atoi(digit);
      pos++;
     } else {
      errorch=MF[pos];
      ERRFLAG=1;
    }
  }

  if(ERRFLAG==0) AddElement(Atom,(*elementcount)++,atomcount);
  else printf("There was an error\n");

  if (ERRFLAG) {
    printf("\nError in format of input...  The character '%c' is invalid\n",errorch);
    return -1;
  } else {
    return 0;
  }

}
 
/*************************************************/
/* FUNCTION CalcFreq - called by main()          */
/*    Could be done with less code, but this     */
/*    saves a few operations.                    */
/*************************************************/
void CMercury8::CalcFreq(Hardklor::complex* FreqData, int Ecount, int NumPoints, int MassRange, int MassShift) {
  
  int    i, j, k, Z;
  double real, imag, freq, X, theta, r, tempr;
  double a, b, c, d;
 
  /* Calculate first half of Frequency Domain (+)masses */
  for (i=0; i<NumPoints/2; i++) {
    
    freq = (double)i/MassRange;
    r = 1;
    theta = 0;
    for (j=0; j<Ecount; j++) {
      Z = AtomicNum[j];
      real = imag = 0;
      for (k=0; k<Element[Z].NumIsotopes; k++) {
				X = TWOPI * Element[Z].IntMass[k] * freq;
				real += Element[Z].IsoProb[k] * cos(X);
				imag += Element[Z].IsoProb[k] * sin(X);
      }
      
      /* Convert to polar coordinates, r then theta */
      tempr = sqrt(real*real+imag*imag);
      r *= pow(tempr,Element[Z].NumAtoms);
      if (real > 0) theta += Element[Z].NumAtoms * atan(imag/real);
      else if (real < 0) theta += Element[Z].NumAtoms * (atan(imag/real) + PI);
      else if (imag > 0) theta += Element[Z].NumAtoms * HALFPI;
      else theta += Element[Z].NumAtoms * -HALFPI;
      
    }  /* end for(j) */
    
    /* Convert back to real:imag coordinates and store */
    a = r * cos(theta);
    b = r * sin(theta);
    c = cos(TWOPI*MassShift*freq);
    d = sin(TWOPI*MassShift*freq);
    FreqData[i].real = a*c - b*d;
    FreqData[i].imag = b*c + a*d;
    
  }  /* end for(i) */
  
  /* Calculate second half of Frequency Domain (-)masses */
  for (i=NumPoints/2; i<NumPoints; i++) {
    
    freq = (double)(i-NumPoints)/MassRange;
    r = 1;
    theta = 0;
    for (j=0; j<Ecount; j++) {
      Z = AtomicNum[j];
      real = imag = 0;
      for (k=0; k<Element[Z].NumIsotopes; k++) {
				X = TWOPI * Element[Z].IntMass[k] * freq;
				real += Element[Z].IsoProb[k] * cos(X);
				imag += Element[Z].IsoProb[k] * sin(X);
      }
      
      /* Convert to polar coordinates, r then theta */
      tempr = sqrt(real*real+imag*imag);
      r *= pow(tempr,Element[Z].NumAtoms);
      if (real > 0) theta += Element[Z].NumAtoms * atan(imag/real);
      else if (real < 0) theta += Element[Z].NumAtoms * (atan(imag/real) + PI);
      else if (imag > 0) theta += Element[Z].NumAtoms * HALFPI;
      else theta += Element[Z].NumAtoms * -HALFPI;
      
    }  /* end for(j) */
    
    /* Convert back to real:imag coordinates and store */
    a = r * cos(theta);
    b = r * sin(theta);
    c = cos(TWOPI*MassShift*freq);
    d = sin(TWOPI*MassShift*freq);
    FreqData[i].real = a*c - b*d;
    FreqData[i].imag = b*c + a*d;
    
  }  /* end of for(i) */
  
}  /* End of CalcFreq() */
  
 
/*************************************************/
/* FUNCTION main() - main block of FFTISO        */
/*************************************************/
//Return Codes:
//	0: Success
//	1: Invalid molecular formula
//	2: Cannot write to file
int CMercury8::GoMercury(char* MolForm, int Charge, const char* filename) {
  
  unsigned int i;
  int	 NumElements=0;			/* Number of elements in molecular formula */
  FILE	 *outfile;			/* output file pointer */
  
  //MolForm is the only required data
  if (strlen(MolForm) == 0) {
    //printf("\nNo molecular formula!\n");
    return 1;
  }
  
  //Parse the formula, check for validity
  if (ParseMF(MolForm,&NumElements) == -1)     {
    MolForm[0] = '\0';
    NumElements = 0;
    return 1;
  }

  //Run the user requested Mercury
  if(bAccMass) AccurateMass(NumElements,Charge);
  else Mercury(NumElements,Charge);
  
  //If the user requested relative abundance, convert data
  if(bRelAbun) RelativeAbundance(FixedData);

  //If the user requested the data to file, output it here.
  if (filename[0]!=0){
    if ((outfile = fopen(filename,"w")) == NULL) {
      printf("\nError - Cannot create file %s\n",filename);
      return 2;
    }
    for(i=0;i<FixedData.size();i++){
      fprintf(outfile,"%lf %lf\n",FixedData[i].mass,FixedData[i].data);
    }
    fclose(outfile);
    
  }
  
  //If the user wants the data on the screen, let it be so.
  if(showOutput){
    for(i=0;i<FixedData.size();i++){
      printf("%.4f\t%.5f\n",FixedData[i].mass,FixedData[i].data);
    }
  }
  
  //Clear all non-user input so object can be reused:
  Reset();
  
  //Output a final, useful message
  if(showOutput) printf("Mercury successful!\n");
  
  return 0;
  
}


//This function resets all atomic states to those at initialization.
//This is so the same object can be reused after performing a calculation
//or after user intervention, such as Enrich().
void CMercury8::Reset(){
  unsigned int i;
  int j;
  
  for (i=0;i<20;i++) AtomicNum[i]=0;
  for (i=0; i<=MAXAtomNo; i++)  Element[i].NumAtoms = 0;
  
  for (i=0;i<EnrichAtoms.size();i++){
    for (j=0;j<Element[EnrichAtoms[i]].NumIsotopes;j++){
      Element[EnrichAtoms[i]].IsoProb[j] = Orig[EnrichAtoms[i]].IsoProb[j];
    }
  }
  EnrichAtoms.clear();
  
}

//This function calculates the accurate mass. It's a bit slower, and it performs poorly
//on the fringes of large proteins because of the estimation of MassRange algorithm, and the
//fact that abundances would be so small as to be nearly indistinguishable from zero.
void CMercury8::AccurateMass(int NumElements, int Charge){
  int 	 i,j, k;
  int	  MassRange;
  int   PtsPerAmu;
  int   NumPoints;			/* Working # of datapoints (real:imag) */
  Hardklor::complex *FreqData;              /* Array of real:imaginary frequency values for FFT */
  double MW;
  double MIMW, tempMW, MolVar, IntMolVar;
  int intMW, MIintMW;
  int dummyLong;
  int dummyInt;


  Hardklor::complex *AltData;
  Hardklor::complex *AltData2;
  int MaxIntMW;
  int PMIintMW;
  int IsoShift;

  Result r;
  vector<Result> vParent;
  vector<Result> vProduct;
  
  //If we made it this far, we have valid input, so calculate molecular weight
  CalcWeights(MW, MIMW, tempMW, intMW, MIintMW, MaxIntMW, IsoShift, NumElements);

  //If the user specified an Echo, output the data now.
  if (showOutput){
    if (Charge != 0) {
      printf("Average Molecular Weight: %.3lf, at m/z: %.3f\n",MW,MW/fabs((double)Charge));
      printf("Average Integer MW: %ld, at m/z: %.3f\n\n",(long)intMW,(float)intMW/fabs((double)Charge));
    } else {
      printf("Average Molecular Weight: %.3lf\n",MW);
      printf("Average Integer MW: %ld\n\n",(long)intMW);
    }
  }

  //Set our parental molecular weight. This will be used as the lower bounds.
  //PMIintMW = (int)round(MIMW);
  PMIintMW = (int)(MIMW+0.5);
  
  //Calculate mass range to use based on molecular variance 
  CalcVariances(&MolVar,&IntMolVar,NumElements);
  CalcMassRange(&MassRange,MolVar,1,1);
  PtsPerAmu = 1;
  
  //Allocate memory for Axis arrays
  NumPoints = MassRange * PtsPerAmu;
  FreqData = new Hardklor::complex[NumPoints];
  
  //Start isotope distribution calculation
  //MH notes: How is this different from using -MW instead of -intMW?
  CalcFreq(FreqData,NumElements,NumPoints,MassRange,-intMW);
  FFT(FreqData,NumPoints,false);

  //Converts Hardklor::complex numbers back to masses
  ConvertMass(FreqData,NumPoints,PtsPerAmu,MW,tempMW,intMW,MIintMW,1,MolVar,IntMolVar);

  //Put our data in the global array
  FixedData.clear();
  for(j=NumPoints/2; j<NumPoints; j++){
    r.data=FreqData[j].real;
    r.mass=FreqData[j].imag;
    FixedData.push_back(r);
  };

  for(j=0; j<NumPoints/2; j++){
    r.data=FreqData[j].real;
    r.mass=FreqData[j].imag;
    FixedData.push_back(r);
  };

  //Convert the data to integers
  MassToInt(FreqData,NumPoints);

  //Shift the lower bound to reflect the range of points that will be in common
  //with all product variants. This is necessary for proteins in which the monoisotopic
  //mass is not visible, thus we must adjust for the first data point in each product
  //if(FixedData.at(0).mass > (PMIintMW+IsoShift)) PMIintMW=(int)round(FixedData.at(0).mass)+IsoShift;
  if(FixedData[0].mass > (PMIintMW+IsoShift)) PMIintMW=(int)(FixedData[0].mass+0.5)+IsoShift;

  //Reduce the parent set to this new boundary
  GetPeaks(FreqData,NumPoints,vParent,PMIintMW,MaxIntMW);

  //Set the upper bound to the max determined in the parent distribution
  //MaxIntMW = (int)round(vParent.at(vParent.size()-1).mass);
  MaxIntMW = (int)(vParent[vParent.size()-1].mass+0.5);

  //set parental masses to 0
  for(i=0;i<(int)vParent.size();i++) vParent[i].mass=0;

  //We have now completed the analysis on the actual protein.
  //Next we compute compositions for each element and each isotope
  for(i=0;i<NumElements;i++){
    
    //Subtract one atom
    Element[AtomicNum[i]].NumAtoms--;

    //Calculate the product masses and variances
    //However, use the same mass range as the parent (calculated above)
    CalcWeights(MW, MIMW, tempMW, intMW, MIintMW, dummyLong, dummyInt, NumElements);
    CalcVariances(&MolVar,&IntMolVar,NumElements);

    //Allocate memory for Axis arrays
    AltData = new Hardklor::complex[NumPoints];
    AltData2 = new Hardklor::complex[NumPoints];
    
    //Start isotope distribution calculation
    CalcFreq(AltData,NumElements,NumPoints,MassRange,-intMW);
    FFT(AltData,NumPoints,false);

    ConvertMass(AltData,NumPoints,PtsPerAmu,MW,tempMW,intMW,MIintMW,1,MolVar,IntMolVar);
    MassToInt(AltData,NumPoints);


    //Add the integer isotope mass
    for(j=0;j<Element[AtomicNum[i]].NumIsotopes;j++){

      //add mass to each point
      for(k=0;k<NumPoints;k++) {
				AltData2[k].imag=AltData[k].imag+Element[AtomicNum[i]].IntMass[j];
				AltData2[k].real=AltData[k].real;
      }

      GetPeaks(AltData2,NumPoints,vProduct,PMIintMW,MaxIntMW);

      //find ratio of abundances, multiply by abundance & number of atoms
      for(k=0;k<(int)vProduct.size();k++){
				vProduct[k].data/=vParent[k].data;
				vProduct[k].data*=Element[AtomicNum[i]].IsoProb[j];
				vProduct[k].data*=(Element[AtomicNum[i]].NumAtoms+1);

				//Add to the real mass for the parent
				vParent[k].mass += vProduct[k].data * Element[AtomicNum[i]].IsoMass[j];
			}
      
    }

		delete [] AltData;
		delete [] AltData2;
    
    //Add back the atom
    Element[AtomicNum[i]].NumAtoms++;
    
  }

  //output accurate masses:
  FixedData.clear();
  for(i=0;i<(int)vParent.size();i++){
    if(vParent[i].data < 0.000001) continue;
    vParent[i].mass=(vParent[i].mass+ProtonMass*Charge)/Charge;
    FixedData.push_back(vParent[i]);
  }

	delete [] FreqData;

}


void CMercury8::ConvertMass(Hardklor::complex* Data, int NumPoints, int PtsPerAmu,
		double MW, double tempMW, int intMW, int MIintMW, int charge,
		double MolVar, double IntMolVar) {

  int i;
  double mass, ratio, CorrIntMW;

  //These are here to prevent warnings...
  (unsigned int)charge;
  (unsigned int)MIintMW;

  if (IntMolVar == 0) ratio = 1;
  else ratio = sqrt(MolVar) / sqrt(IntMolVar);
  
  CorrIntMW = tempMW * ratio;
  for (i=NumPoints/2; i<NumPoints; i++) {
    mass = (double)(i-NumPoints)/PtsPerAmu + intMW;
    mass *= ratio;
    mass += MW - CorrIntMW;
    //mass /= charge;
    Data[i].imag = mass;
  };
  
  for (i=0; i<NumPoints/2; i++) {
    mass = (double)i/PtsPerAmu + intMW;
    mass *= ratio;
    mass += MW - CorrIntMW;
    //mass /= charge;
    Data[i].imag = mass;
  };

};

void CMercury8::MassToInt(Hardklor::complex* Data, int NumPoints) {

  int i, mass;

  for (i=NumPoints/2; i<NumPoints; i++) {

    //Since rounding poses problems, adjust when deviant values occur
    //This assumes that the average width >1
	mass = (int)(Data[i].imag+0.5);
    if(i>NumPoints/2){
      if(mass == Data[i-1].imag) mass++;

    };

    Data[i].imag = mass;
  };
  
  for (i=0; i<NumPoints/2; i++) {

    //Since rounding poses problems, adjust when deviant values occur
    //This assumes that the average width >1
	mass = (int)(Data[i].imag+0.5);	
    if(i>0){
      if(mass == Data[i-1].imag) mass++;
    };

    Data[i].imag = mass;
  };
  
};


void CMercury8::GetPeaks(Hardklor::complex* Data, int NumPoints, vector<Result>& v, 
			 int lower, int upper){
  int i;
  Result r;
  v.clear();

  for(i=NumPoints/2; i<NumPoints; i++){
    if(Data[i].imag>=lower && Data[i].imag<=upper){
      r.data=Data[i].real;
      r.mass=Data[i].imag;
      v.push_back(r);
    };
  };

  for(i=0; i<NumPoints/2; i++){
    if(Data[i].imag>=lower && Data[i].imag<=upper){
      r.data=Data[i].real;
      r.mass=Data[i].imag;
      v.push_back(r);
    };
  };

};

void CMercury8::RelativeAbundance(vector<Result>& v){

  unsigned int i;
  double max=0;
	double sum=0;

	FracAbunData.clear();
  
  /* Normalize intensity to 0%-100% scale */
  for (i=0; i<v.size(); i++) {
		if (v[i].data > max) max = v[i].data;
		//also, store fractional abundances anyway (useful for Hardklor, which uses both);
		FracAbunData.push_back(v[i]);
		sum+=v[i].data;
	}
  for (i=0; i<v.size(); i++) {
		v[i].data = 100 * v[i].data/max;
		FracAbunData[i].data/=sum;
	}

}

void CMercury8::CalcWeights(double& MW, double& MIMW, double& tempMW, int& intMW, 
			    int& MIintMW, int& MaxIntMW, int& IsoShift, int NumElements){

  int j, k, Z;

  MW = MIMW = tempMW = 0; intMW = MIintMW = MaxIntMW = 0;
  IsoShift=0;

  for (j=0; j<NumElements; j++) {
    Z = AtomicNum[j];
    for (k=0; k<Element[Z].NumIsotopes; k++) {
      MW += Element[Z].NumAtoms * Element[Z].IsoMass[k] * Element[Z].IsoProb[k];
      tempMW += Element[Z].NumAtoms * Element[Z].IntMass[k] * Element[Z].IsoProb[k];
      if (k==0) {
				MIMW += Element[Z].NumAtoms * Element[Z].IsoMass[k];
				MIintMW += Element[Z].NumAtoms * Element[Z].IntMass[k];
      }
      if(k==Element[Z].NumIsotopes-1){

				//IsoShift is only used for accurate masses as an indication of how
				//much each product distribution will differ in mass range from the
				//parent distribution.
				if((Element[Z].IntMass[k]-Element[Z].IntMass[0]) > IsoShift){
					IsoShift = Element[Z].IntMass[k]-Element[Z].IntMass[0];
				}

				MaxIntMW += Element[Z].NumAtoms * Element[Z].IntMass[k];
      }
    }
  }
  
  MW -= ElectronMass;
  tempMW -= ElectronMass;
  MIMW -= ElectronMass;
  intMW = (int)(tempMW+0.5);

};

//Calculates Mercury as done originally
void CMercury8::Mercury(int NumElements, int Charge) {
  
  int 	  j;
  int	  MassRange;
  int     PtsPerAmu;
  int    NumPoints;		       // Working # of datapoints (real:imag) 
  Hardklor::complex *FreqData;              // Array of real:imaginary frequency values for FFT 
  double   MW;
  double  MIMW, tempMW, MolVar, IntMolVar;
  int    intMW, MIintMW, MaxIntMW;
  clock_t start, end;
  float   timex, seconds;
  int	  minutes, dummyInt;
  Result r;

	//Calculate molecular weight
  CalcWeights(MW, MIMW, tempMW, intMW, MIintMW, MaxIntMW, dummyInt, NumElements);
   
  //If the user specified an Echo, output the data now.
  if (showOutput){
    if (Charge != 0) {
      printf("Average Molecular Weight: %.3lf, at m/z: %.3f\n",MW,MW/fabs((double)Charge));
      printf("Average Integer MW: %ld, at m/z: %.3f\n\n",(long)intMW,(float)intMW/fabs((double)Charge));
    } else {
      printf("Average Molecular Weight: %.3lf\n",MW);
      printf("Average Integer MW: %ld\n\n",(long)intMW);
    }
  }
 
  //Calculate mass range to use based on molecular variance 
  CalcVariances(&MolVar,&IntMolVar,NumElements);
  CalcMassRange(&MassRange,MolVar,Charge,1);
  PtsPerAmu = 1;

  monoMass+=(Charge*(ProtonMass));
  monoMass/=Charge;
  if(showOutput) printf("True MonoMass: %.8lf\n",monoMass);

  //Allocate memory for Axis arrays
  NumPoints = MassRange * PtsPerAmu;
  FreqData = new Hardklor::complex[NumPoints];
  
  //Start isotope distribution calculation
  //MH notes: How is this different from using -MW instead of -intMW?
  start = clock();
  CalcFreq(FreqData,NumElements,NumPoints,MassRange,-intMW);
  FFT(FreqData,NumPoints,false);
  end = clock();

  //Output the results if the user requested an Echo.
  if(showOutput){
    timex = (float)(end - start) / CLOCKS_PER_SEC;
    minutes = (int)(timex/60);
    seconds = timex - (60*minutes);
    printf("Calculation performed in %d min %.3f sec\n",minutes,seconds);
  }
	
  //Not sure why we do this...
  if (Charge == 0) Charge = 1;

  //Convert Hardklor::complex numbers to masses
  ConvertMass(FreqData,NumPoints,PtsPerAmu,MW,tempMW,intMW,MIintMW,(int)fabs((double)Charge),MolVar,IntMolVar);

  //Put our data in the global array
  //This performs a bit of computation to eliminate meaningless data.
  FixedData.clear();
  for(j=NumPoints/2; j<NumPoints; j++){
		if((int)(FreqData[j].imag+0.5)<MIintMW) continue;
    r.data=FreqData[j].real;
		r.mass=(FreqData[j].imag+(ProtonMass)*Charge)/Charge;
		if( (monoMass-r.mass)*Charge > 0.5 ) continue;
    r.data=FreqData[j].real;
		r.mass=(FreqData[j].imag+ProtonMass*Charge)/Charge;
    FixedData.push_back(r);
  };

  for(j=0; j<NumPoints/2; j++){
	if((int)(FreqData[j].imag+0.5)>MaxIntMW) continue;
    if(FreqData[j].real<0) break;
    r.data=FreqData[j].real;
		r.mass=(FreqData[j].imag+ProtonMass*Charge)/Charge;
    FixedData.push_back(r);
  };
  
  //Clean up the memory
  delete [] FreqData;
 
};


void CMercury8::AccMass(bool b) {
  bAccMass = b;
}

void CMercury8::RelAbun(bool b) {
  bRelAbun = b;
}

double CMercury8::getZeroMass(){
  return zeroMass;
}

double CMercury8::getMonoMass(){
  return monoMass;
}
