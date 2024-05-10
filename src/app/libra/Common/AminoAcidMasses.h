/*
 * elemental masses from http://www.unimod.org/masses.html
 * bMassType:  0=AVG, 1=MONO 
 */

void INITIALIZE_MASS(double *pdMassAA,
        int bMonoMass)
{
   double H, O, C, N, P, S, Se;

   /*
    * elemental masses from http://www.unimod.org/masses.html
    */
   if (bMonoMass) /* monoisotopic masses */
   {
      H = pdMassAA['h'] =  1.007825035; /* hydrogen */
      O = pdMassAA['o'] = 15.99491463;  /* oxygen */
      C = pdMassAA['c'] = 12.0000000;   /* carbon */
      N = pdMassAA['n'] = 14.0030740;   /* nitrogen */
      P = pdMassAA['p'] = 30.973762;    /* phosphorus */
      S = pdMassAA['s'] = 31.9720707;   /* sulphur */
     Se = pdMassAA['e'] = 79.9165196;   /* selenium */
   }
   else  /* average masses */
   {
      H = pdMassAA['h'] =  1.00794;     /* hydrogen */
      O = pdMassAA['o'] = 15.9994;      /* oxygen */
      C = pdMassAA['c'] = 12.0107;      /* carbon */
      N = pdMassAA['n'] = 14.0067;      /* nitrogen */
      P = pdMassAA['p'] = 30.973761;    /* phosporus */
      S = pdMassAA['s'] = 32.065;       /* sulphur */
     Se = pdMassAA['e'] = 78.96;        /* selenium */
   }

   pdMassAA['G'] = C*2  + H*3  + N   + O ;
   pdMassAA['A'] = C*3  + H*5  + N   + O ;
   pdMassAA['S'] = C*3  + H*5  + N   + O*2 ;
   pdMassAA['P'] = C*5  + H*7  + N   + O ;
   pdMassAA['V'] = C*5  + H*9  + N   + O ;
   pdMassAA['T'] = C*4  + H*7  + N   + O*2 ;
   pdMassAA['C'] = C*3  + H*5  + N   + O   + S ;
   pdMassAA['L'] = C*6  + H*11 + N   + O ;
   pdMassAA['I'] = C*6  + H*11 + N   + O ;
   pdMassAA['N'] = C*4  + H*6  + N*2 + O*2 ;
   pdMassAA['D'] = C*4  + H*5  + N   + O*3 ;
   pdMassAA['Q'] = C*5  + H*8  + N*2 + O*2 ;
   pdMassAA['K'] = C*6  + H*12 + N*2 + O ;
   pdMassAA['E'] = C*5  + H*7  + N   + O*3 ;
   pdMassAA['M'] = C*5  + H*9  + N   + O   + S ;
   pdMassAA['H'] = C*6  + H*7  + N*3 + O ;
   pdMassAA['F'] = C*9  + H*9  + N   + O ;
   pdMassAA['R'] = C*6  + H*12 + N*4 + O ;
   pdMassAA['Y'] = C*9  + H*9  + N   + O*2 ;
   pdMassAA['W'] = C*11 + H*10 + N*2 + O ;
   pdMassAA['U'] = C*3  + H*5  + N   + O   + Se;

   pdMassAA['O'] = C*12  + H*19 + N*3 + O*2 ;
   pdMassAA['B'] = 0.0;
   pdMassAA['J'] = 0.0;
   pdMassAA['X'] = 0.0;
   pdMassAA['Z'] = 0.0;
   
} /*ASSIGN_MASS*/
