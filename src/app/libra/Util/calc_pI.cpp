/***************************************************************************
*
*  FILENAME      :   pi.c
*
*  DESCRIPTION   :   Given a protein sequence, these functions return its pI.
* 
*  ARGUMENTS     :   char * seq; string containing the sequence, 
*                    all upper case letters, no leading or trailing blanks.
*        
*  AUTHORS       :   ROA from Amos' BASIC procedure
*
*  VERSION       :   1.6
*  DATE          :   1/25/95
*  
*  Copyright 1993 by Melanie/UIN/HCUG. All rights reserved.
*
***************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <math.h>

#define PH_MIN  0       /* minimum pH value */
#define PH_MAX  14      /* maximum pH value */
#define MAXLOOP 2000    /* maximum number of iterations */
#define EPSI    0.0001  /* desired precision */


/* the 7 amino acid which matter */
static int R = 'R' - 'A',
           H = 'H' - 'A',
           K = 'K' - 'A',
           D = 'D' - 'A',
           E = 'E' - 'A',
           C = 'C' - 'A',
           Y = 'Y' - 'A';

/*
 *  table of pk values : 
 *  Note: the current algorithm does not use the last two columns. Each 
 *  row corresponds to an amino acid starting with Ala. J, O and U are 
 *  inexistant, but here only in order to have the complete alphabet.
 *
 *          Ct    Nt   Sm     Sc     Sn
 */

static double  pk [26][5] = {
/* A */    3.55, 7.59, 0.0  , 0.0  , 0.0   ,
/* B */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* C */    3.55, 7.50, 9.00 , 9.00 , 9.00  ,
/* D */    4.55, 7.50, 4.05 , 4.05 , 4.05  ,
/* E */    4.75, 7.70, 4.45 , 4.45 , 4.45  ,
/* F */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* G */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* H */    3.55, 7.50, 5.98 , 5.98 , 5.98  ,
/* I */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* J */    0.00, 0.00, 0.0  , 0.0  , 0.0   ,
/* K */    3.55, 7.50, 10.00, 10.00, 10.00 ,
/* L */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* M */    3.55, 7.00, 0.0  , 0.0  , 0.0   ,
/* N */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* O */    0.00, 0.00, 0.0  , 0.0  , 0.0   ,
/* P */    3.55, 8.36, 0.0  , 0.0  , 0.0   ,
/* Q */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* R */    3.55, 7.50, 12.0 , 12.0 , 12.0  ,
/* S */    3.55, 6.93, 0.0  , 0.0  , 0.0   ,
/* T */    3.55, 6.82, 0.0  , 0.0  , 0.0   ,
/* U */    0.00, 0.00, 0.0  , 0.0  , 0.0   ,
/* V */    3.55, 7.44, 0.0  , 0.0  , 0.0   ,
/* W */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* X */    3.55, 7.50, 0.0  , 0.0  , 0.0   ,
/* Y */    3.55, 7.50, 10.00, 10.00, 10.00 ,
/* Z */    3.55, 7.50, 0.0  , 0.0  , 0.0   };

static double EXP10(double value);
double COMPUTE_PI(char *seq,unsigned long seq_length, int charge_increment);

static double EXP10(double value)
{
   return( pow(10.0,value) );
}


double COMPUTE_PI(char *seq,unsigned long seq_length, int charge_increment)
{
   int             comp[26];    /* Amino acid composition of the protein */
   register int    nterm_res,   /* N-terminal residue */
                   cterm_res;   /* C-terminal residue */
   register unsigned long i;
   register double charge,
                   ph_mid,
                   ph_min,
                   ph_max;
   double          cter,
                   nter;
   double          carg,
                   clys,
                   chis,
                   casp,
                   cglu,
                   ctyr,
                   ccys;
  
   memset((char *)comp, 0, 26 * sizeof(int));

   for (i = 0; i < seq_length; i++)        /* compute the amino acid composition */
   {
      comp[seq[i] - 'A']++;
   }
  
   nterm_res = seq[0] - 'A';               /* Look up N-terminal residue */
   cterm_res = seq[seq_length-1] - 'A';    /* Look up C-terminal residue */
  
   ph_min = PH_MIN;
   ph_max = PH_MAX;
  
   for (i = 0, charge = 1.0; i<MAXLOOP && (ph_max - ph_min)>EPSI; i++)
   {
      ph_mid = ph_min + (ph_max - ph_min) / 2.0;
    
      cter = EXP10(-pk[cterm_res][0]) / (EXP10(-pk[cterm_res][0]) + EXP10(-ph_mid));
      nter = EXP10(-ph_mid) / (EXP10(-pk[nterm_res][1]) + EXP10(-ph_mid));
    
      carg = comp[R] * EXP10(-ph_mid) / (EXP10(-pk[R][2]) + EXP10(-ph_mid));
      chis = comp[H] * EXP10(-ph_mid) / (EXP10(-pk[H][2]) + EXP10(-ph_mid));
      clys = comp[K] * EXP10(-ph_mid) / (EXP10(-pk[K][2]) + EXP10(-ph_mid));
    
      casp = comp[D] * EXP10(-pk[D][2]) / (EXP10(-pk[D][2]) + EXP10(-ph_mid));
      cglu = comp[E] * EXP10(-pk[E][2]) / (EXP10(-pk[E][2]) + EXP10(-ph_mid));
    
      ccys = comp[C] * EXP10(-pk[C][2]) / (EXP10(-pk[C][2]) + EXP10(-ph_mid));
      ctyr = comp[Y] * EXP10(-pk[Y][2]) / (EXP10(-pk[Y][2]) + EXP10(-ph_mid));
    
      charge = carg + clys + chis + nter + charge_increment 
         - (casp + cglu + ctyr + ccys + cter);
    
      if (charge > 0.0)
      {
         ph_min = ph_mid;
      }
      else
      {
         ph_max = ph_mid;
      }
   }
  
   return(ph_mid);
}

