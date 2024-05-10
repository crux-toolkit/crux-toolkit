/*

Program       : spectStrct                                                       
Author        : Brian Pratt, Insilicos LLC                                                       
Date          : 7-6--6 

This is a consolidation and optimzation of the Savitsky-Golay smoothing 
code that was widely cut and pasted throughout the TPP project, which is
Copyright(C) 2003 ISB .

Optimizations are Copyright (C) 2006 Insilicos LLC


This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Brian Pratt
www.insilicos.com

*/

#ifndef SPECTSCRT_H_INCL_
#define SPECTSCRT_H_INCL_

/*
  structure of spectrum
*/
#include <math.h>
#include "constants.h" // for _ASAPRATIO_MXQ_ defn
#include "Array.h"

class spectStrct {
public:
	spectStrct():size(0),xval(NULL),yval(NULL) {
		initCacheVals();
	}
	spectStrct(int sz):size(0),xval(NULL),yval(NULL) {
		initCacheVals();
		setsize(sz);
	}

   ~spectStrct() {
      clear();
   }

  int size;
  double *xval;
  double *yval;
private:
   mutable double *root4;  // for cacheing expensive root calcs
   mutable Array<double *>pows; // to avoid beating on heap manager
   mutable double *sums;
   mutable double *vec;
   mutable double **mtrx;
   mutable int *indx;
   mutable int maxcount;
   void initCacheVals();
   void clear_cachevals(); // blow cached calc aids
   void freeDataFilterCache() const; // release memory for cached calcs
   void setDataFilterCacheSize(int order, int count) const; // resize memory for cached calcs
public:
   void setsize(int sz) {
      clear();
      size = sz;
      xval = (double *) calloc(size, sizeof(double));
      yval = (double *) calloc(size, sizeof(double));
   }
  double dataFilter(int dtIndx, int lower, int upper, int order) const;  // used in smoothing
  void set_yvals(double *src) {
      memmove(yval,src,size*sizeof(double));
      clear_cachevals(); // blow cached calc aids
  }
  void clear() {
     free(xval);
     free(yval);
     clear_cachevals(); // blow cached calc aids
     xval = NULL;
     yval = NULL;
  }
  spectStrct & operator = (const spectStrct &rhs) {
	  clear();
	  setsize(rhs.size);
	  memmove(xval,rhs.xval,size*sizeof(double));
	  set_yvals(rhs.yval);
	  return *this;
  }

	// This function smoothes rough spectrum at a specific point.
	double smoothDataPtFlx(int dtIndx, 
			       double range, double threshold, int smoothWindow=10) const;
	// This function smoothes rough spectrum within a x range.
	void smoothSpectFlx(double range, int repeats, int smoothWindow=10, bool wavelet=true);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  This function performs LU decomposition.
*/
void myLUDcmp(double **mtrx, int order, int *indx, double *d);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  This function performs LU backsubstition.
*/
void myLUBksb(double **mtrx, int order, int *indx, double *vec);


/*
  structure of LC spectra
*/
class lcSpectStrct {
public:
   lcSpectStrct():expScanNums(NULL) {
   }
   ~lcSpectStrct() {
      free(expScanNums);
   }
  spectStrct rawSpectrum[_ASAPRATIO_MXQ_][2]; 
  spectStrct fitSpectrum[_ASAPRATIO_MXQ_][2]; 
  long *expScanNums;
};


#endif // SPECTSCRT_H_INCL_
