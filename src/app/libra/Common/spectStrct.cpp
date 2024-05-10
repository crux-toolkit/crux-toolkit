/*

Program       : spectStrct                                                       
Author        : Brian Pratt, Insilicos LLC                                                       
Date          : 7-6-06 

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

#include <stdlib.h>
#include "spectStrct.h"

// this code gets called a LOT in ASAPRatio, so we're caching some calculations and avoiding heap thrash

void spectStrct::initCacheVals() {
	root4= NULL;
	sums=NULL;
	vec=NULL;
	mtrx=NULL;
	indx=NULL;
}

void spectStrct::clear_cachevals() {
	free(root4); // blow cached calc aids
	root4 = NULL;
	freeDataFilterCache();
}

void spectStrct::freeDataFilterCache() const {
	for (int pi=pows.size();pi--;) {
		delete[] pows[pi];
	}
	free(sums);
	sums = NULL;
	free(vec);
	vec = NULL;
	for (int di=pows.size()/2; di--;) {
		free(mtrx[di]);
	}
	free(mtrx);
	mtrx = NULL;
	pows.clear();
	free(indx);
	indx = NULL;
}

void spectStrct::setDataFilterCacheSize(int order, int count) const {

  if (!root4) {
     root4 = (double *)calloc(size+1,sizeof(double));
  }
	if (pows.size() < 2*(order+1)) { // resize it 
		for (int pi=pows.size();pi--;) {
			delete[] pows[pi];
			pows[pi] = NULL;
		}
		free(vec);
		vec = (double *) malloc((order+1) * sizeof(double));
		free(indx);
		indx = (int *) malloc((order+1) * sizeof(int));
		free(sums);
		sums = (double *)malloc(2*(order+1)*sizeof(double));
		for (int m = pows.size()/2;m--;) {
			free(mtrx[m]);
		}
		free(mtrx);
		mtrx = (double **) malloc((order+1)* sizeof(double *));
		for (int ii = order+1; ii--; ) {
			mtrx[ii] = (double *) malloc((order+1)*sizeof(double));
		}
		pows.setSize(2*(order+1));
		maxcount = 0;
	}
	memset(sums,0,2*(order+1)*sizeof(double));
	if (count > maxcount) {
		count *= 2; // over allocate to avoid lots of reallocs
		for (int i=pows.size();i--;) {
			if (maxcount) {
				delete[] pows[i];
			}
			pows[i] = new double[count];
		}
		maxcount = count;
	}
}

double spectStrct::dataFilter(int dtIndx, 
                              int lower, int upper, int order) const
{
  double x = xval[dtIndx];
  double sum, d;
  
  // do trivial orders 0 and 1
  const int count = 1+upper-lower; // range is lower...upper inclusive
  setDataFilterCacheSize(order,count);  // (re)size the pow() calc area
  for (int n = count, mi = upper; n--;) {
	  pows[0][n] = 1.0;
	  if (!root4[mi]) {
		  // cache those pow(yval,.25) calcs, they're expensive
		  if (yval[mi] > 1.) {
			  root4[mi] = sqrt(sqrt(yval[mi]));
		  } else {
			  root4[mi] = 1.0;
		  }
	  }
	  sums[1] += root4[mi]*(pows[1][n] = xval[mi]-x);
	  sums[0] += root4[mi--];
  }
  // with that seed, do higher orders
  for (int ord=2;ord<2*(order+1);ord++) {
     for (int k = count, m = upper;k--;) {
		sums[ord] += root4[m--]*(pows[ord][k] = pows[ord/2][k]*pows[(ord&1)+(ord/2)][k]);
     }
  }
  for (int oi = 0; oi <= order; ++oi) {
    for (int oj = oi; oj <= order; ++oj) {
      mtrx[oi][oj] = mtrx[oj][oi] = sums[oi+oj];
    }
  }

  // LU decomposition
  myLUDcmp(mtrx, order+1, indx, &d);

  // get "vec"
  for (int i = order+1; i--;) {
    sum = 0.0;
    for (int k = count, m = upper;k--;m--) {
       sum += pows[i][k]*yval[m]*root4[m];
    }
    vec[i] = sum;
  }

  // LU backsubstition
  myLUBksb(mtrx, order+1, indx, vec);
  x = vec[0];



  return x;

}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
  This function smoothes rough spectrum within a x range.
*/
void spectStrct::smoothSpectFlx(double range, int repeats, int smoothWindow, bool wavelet)
{

  int size = this->size;
  double *temp;
  double threshold;

  double *C0, *h, *g, *C1,*C2,*C4;                //sqrt(2)/2,-sqrt(2)/2,0,0
  double **d;
  int i,j,m,loop=0;
  double x,y,k=0;

  if (wavelet) {
    // use wavelet smoothing from WaveletQuant implementation

    // begin code from WaveletQuant

    //------------------------------wavelet-----------------------------

    if (size<16)
      {
	m=size;
	while(m>1) {
	  m=m/2;
	  loop++;
	}
      }
    else loop=4;
 
    d=(double **)malloc(loop*sizeof(double *));
    C0=(double *)calloc(size, sizeof(double));

    for(i=0;i<size;i++)
      {
	C0[i]=this->yval[i];
      }
    h=(double *)calloc(size, sizeof(double));
    g=(double *)calloc(size, sizeof(double));
    for(i=0;i<size;i++) 
      {
	h[i]=0;
      }
    for(i=0;i<size;i++)
      {
	g[i]=0;
      }	
    //db4
    h[0]=0.2304;h[1]=0.7148;h[2]=0.6309;h[3]=-0.0280;h[4]=-0.1870;h[5]=0.0308;h[6]=0.0329;h[7]=-0.0106;
    int t=8;
    int flag=-1;
    for(i=0;i<t;i++)
      {
	g[i]=flag*h[t-1-i];
	flag*=(-1);
      }
    //loop (unreadable)
    C2=(double *)calloc(size,sizeof(double));
    C1=(double *)calloc(size,sizeof(double));
    for(m=0;m<loop;m++)
      {
	for(i=0;i<size;i++)
	  {
	    C1[i]=0;
	    C2[i]=0;
	  }
	for(i=0; i<size; i++)
	  {
	    for(j=0; j<size; j++)
	      {
		if(j-2*i<0) 
		  x=0;
		else
		  x=C0[j]*h[j-2*i];
		C1[i]+=x;

		if(j-2*i<0) 
		  y=0;
		else
		  y=C0[j]*g[j-2*i];
		C2[i]+=y;
	      }
	  }

	for(i=0; i<size; i++)
	  {
	    C0[i]=C1[i];
	  }
	d[m]=(double *)malloc(size*sizeof(double));
        for(i=0; i<size; i++)
	  {
	    d[m][i]=C2[i];
	  }
      }
    // (unreadable)
    double MAD,T,P_temp,avaY;
    double sum=0.0;
    long N,P;
    N=this->size;
    temp=(double *)malloc(N*sizeof(double));

    for(i=0;i<N;i++) {
      temp[i]=0;
      sum+=d[0][i];    
    }
    avaY=sum/N;       
    
    for(i=0;i<N;i++)  {
      if(d[0][i]-avaY >=0 ) 
	temp[i]=d[0][i]-avaY;
      else
	temp[i]=d[0][i];
    }

    sum = 0.;

    for(i=0;i<N;i++) {
      sum+=temp[i];
    }
  

    MAD=sum/N;
    T=MAD/0.6745*sqrt(2*log((double)N));
    //free(temp);
    // (unreadable)
    for(m=0;m<loop;m++)
      for(i=0; i<size; i++)
	{
	  if(fabs(d[m][i])<=T)
	    d[m][i]=0;
	}
    
    C4=(double *)calloc(size,sizeof(double));

    // (unreadable)
    for(m=loop-1;m>=0;m--)
      {
	for(i=0;i<size;i++)
	  C4[i]=0;
	for(i=0; i<size; i++)
	  for(j=0; j<size; j++)
	    {
	      if(i-2*j<0) 
		{
		  x=0;
		  y=0;
		}
	      else	
		{
		  x=C0[j]*h[i-2*j];
		  y=d[m][j]*g[i-2*j];
		}
	      C4[i]+=x+y;
	    }
	for(i=0; i<size; i++)
	  C0[i]=C4[i];
      }
    for(i=0;i<size;i++)
      {
	this->yval[i]=C0[i];
      }
    free(C0);
    free(h);
    free(g);
    free(C1);
    free(C2);
    free(C4);
    for(i=0;i<loop;i++)
      free(d[i]);
    free(d);
    //------------------------------------------------------------------------------------------

    // end code from WaveletQuant


#ifndef DEBUG_NOFREE
    free(temp);
#endif

  }
  else {
    temp = (double *) calloc(size, sizeof(double));
    
    threshold = this->yval[0];
    for (i = 0; i < size; ++i) {
      if(this->yval[i] < threshold) {
	threshold = this->yval[i];
      }
    }
    
    for (j = 0; j < repeats; ++j) {
      for (i = 0; i < size; ++i) {
	temp[i] = this->smoothDataPtFlx(i, range, threshold, smoothWindow);
      }
      this->set_yvals(temp);
    }
    
    free(temp);
    
  }

    return;
}


/////////////////////////////////////////////////////////////////////////
/*
 This function smoothes rough spectrum at a specific point.
*/
double spectStrct::smoothDataPtFlx(int dtIndx, 
		       double range, double threshold, int smoothWindow) const
{

  int order = 4;
  int lower, upper;
  double value;

  // get boundary
  lower = dtIndx;
  while(lower > 0 
	&& this->xval[dtIndx]-this->xval[lower] < range) 
    --lower;
  upper = dtIndx;
  while(upper < this->size-1
	&& this->xval[upper]-this->xval[dtIndx] < range) 
    ++upper;

  // get filter value
  while(upper-lower < smoothWindow) {
    if(lower > 0)
      --lower;
    if(upper < this->size-1) 
      ++upper;
    if(lower <= 0 
       && upper >= this->size-1) 
      break;
  }
  order = order < upper-lower-1 ? order : upper-lower-1;
  if(order < 1)
    return this->yval[dtIndx];
  else
    value = this->dataFilter( dtIndx, lower, upper, order); 

  if(value > threshold) 
    return value;
  else
    return threshold;

}  



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  This function performs LU decomposition.
*/
void myLUDcmp(double **mtrx, int order, int *indx, double *d)
{
  double big, dum, sum, temp;
  double *vv;
  int i, j, k, imax=-1;

  vv = (double *) calloc(order, sizeof(double));

  *d = 1.0;
  for (i = 0; i < order; i++) {
    big = 0.0;
	for (j = 0; j < order; j++) {
		if ((temp=fabs(mtrx[i][j])) > big) {
			big = temp;
		}
	}
    vv[i] = 1.0/big;
  } 

  for (j = 0; j < order; j++) {
    for (i = 0; i < j; i++) {
      sum = mtrx[i][j];
	  for (k = 0; k < i; k++) {
		  sum -= mtrx[i][k]*mtrx[k][j];
	  }
      mtrx[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < order; i++) {
      sum = mtrx[i][j];
	  for (k = 0; k < j; k++) {
	sum -= mtrx[i][k]*mtrx[k][j];
	  }
      mtrx[i][j] = sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < order; k++) {
	dum = mtrx[imax][k];
	mtrx[imax][k] = mtrx[j][k];
	mtrx[j][k] = dum;
      }
      *d *= -1;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
	if (mtrx[j][j] == 0.0) {
		mtrx[j][j] = 1.e-20;
	}
    if (j != order) {
      dum = 1.0/(mtrx[j][j]);
	  for (i = j+1; i < order; i++) {
		  mtrx[i][j] *= dum;
	  }
    }
  }
  free(vv);
 
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  This function performs LU backsubstition.
*/
void myLUBksb(double **mtrx, int order, int *indx, double *vec)
{
  double sum;
  int ii, ip;
  int i, j;


  ii = -1;
  for (i = 0; i < order; i++) {
    ip = indx[i];
    sum = vec[ip];
    vec[ip] = vec[i];
	if (ii != -1) {
		for (j = ii; j < i; j++) {
			sum -= mtrx[i][j]*vec[j];
		}
	} else if (sum) {
		ii = i;
	}
    vec[i] = sum;
  }
  for (i = order-1; i >= 0; i--) {
    sum = vec[i];
	for (j = i+1; j < order; j++) {
		sum -= mtrx[i][j]*vec[j];
	}
    vec[i] = sum/mtrx[i][i];
  }

  return;
}

