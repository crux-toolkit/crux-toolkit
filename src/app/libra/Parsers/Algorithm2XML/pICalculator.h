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
#ifndef PICALC_H
#define PICALC_H

#include "Common/sysdepend.h"

#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "Common/Array.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Validation/MixtureDistribution/ICATMixtureDistr.h"

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

#define PK_X 26
#define PK_Y 5
const double  pk [PK_X][PK_Y] = {
/* A */   { 3.55, 7.59, 0.0  , 0.0  , 0.0   },
/* B */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* C */   { 3.55, 7.50, 9.00 , 9.00 , 9.00  },
/* D */   { 4.55, 7.50, 4.05 , 4.05 , 4.05  },
/* E */   { 4.75, 7.70, 4.45 , 4.45 , 4.45  },
/* F */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* G */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* H */   { 3.55, 7.50, 5.98 , 5.98 , 5.98  },
/* I */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* J */   { 0.00, 0.00, 0.0  , 0.0  , 0.0   },
/* K */   { 3.55, 7.50, 10.00, 10.00, 10.00 },
/* L */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* M */   { 3.55, 7.00, 0.0  , 0.0  , 0.0   },
/* N */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* O */   { 0.00, 0.00, 0.0  , 0.0  , 0.0   },
/* P */   { 3.55, 8.36, 0.0  , 0.0  , 0.0   },
/* Q */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* R */   { 3.55, 7.50, 12.0 , 12.0 , 12.0  },
/* S */   { 3.55, 6.93, 0.0  , 0.0  , 0.0   },
/* T */   { 3.55, 6.82, 0.0  , 0.0  , 0.0   },
/* U */   { 0.00, 0.00, 0.0  , 0.0  , 0.0   },
/* V */   { 3.55, 7.44, 0.0  , 0.0  , 0.0   },
/* W */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* X */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* Y */   { 3.55, 7.50, 10.00, 10.00, 10.00 },
/* Z */   { 3.55, 7.50, 0.0  , 0.0  , 0.0   }};
class pICalculator {
  friend class pIMixtureDistr;
  friend class VariableOffsetpIMixtureDistr;
  friend class KernelDensityPIMixtureDistr;

 public:

  pICalculator();
  ~pICalculator();
  pICalculator(int run_idx);
  // seq must be stripped of all modifications
  double Peptide_pI(const char* seq, ModificationInfo* modinfo, double pI);
  double Peptide_pI(const char* seq, ModificationInfo* modinfo);
  // can have modifications
  //double ModifiedPeptide_pI(char* seq, int offset, int len);
  double ModifiedPeptide_pI(const char* seq, int offset, int len);
  void calc_pIstats();
  Boolean recalc_pIstats(Array<double>* probs, double min_prob, Array<int>* ntts, int min_ntt);
  //Boolean recalc_pIstats(Array<double>* probs, double min_prob); 
  double getpIScore(double pI, const char* pep);
  void recalc_pIstats(Array<double>* probs);
  void write_pIstats(ostream& out);
  int getRunIdx();
  double getRunPIMean();
  double getRunPIStdDev();

 protected:
  void common_ctor_init(int run_idx);
  int run_idx_;
  double run_pI_sum_;
  int run_pI_count_;
  int run_pI_used_count_;
  double run_pI_mean_;
  double run_pI_stddev_;
  Array<double> pIs_;
  Array<char*>* peps_;

  double COMPUTE_PI(const char *seq,unsigned long seq_length, int charge_increment, ModificationInfo* modinfo);

};

#endif
