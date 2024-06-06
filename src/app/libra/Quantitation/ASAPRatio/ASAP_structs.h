#ifndef ASAP_STRUCTS_H
#define ASAP_STRUCTS_H

#include "Common/constants.h"

// structure of data
class dataStrct {
 public:
  int indx; // 0 new, 1 analyzed, -1 invalid
  double ratio[2]; // ratio and error 
  double inv_ratio[2]; // inverted ratio and error 
  int dataNum; // data number
  int *dataIndx; // index for data
  int *dataCnts; // whether to count the data: 0 no, 1 yes
  double weight; // data weight
  int bofIndx;
  int msms_run_idx;

  void sort_for_output() {
    for (int i=dataNum;i--;) {
      for (int j=i;j--;) {
	if (dataIndx[i]<dataIndx[j]) {
	  int tmp = dataIndx[i];
	  dataIndx[i] = dataIndx[j];
	  dataIndx[j] = tmp;
	  tmp = dataCnts[i];
	  dataCnts[i] = dataCnts[j];
	  dataCnts[j] = tmp;
	}
      }
    }
  }
};

// structure of unique sequence
class seqDataStrct {
 public:
  int indx; // 0 new, 1 analyzed, -1 invalid
  double ratio[2]; // ratio and error 
  double inv_ratio[2]; // inverted ratio and error
  int dataNum; // data number
  dataStrct *peaks; // unique peaks
  int *dataCnts; // whether to count the data: 0 no, 1 yes
  double weight; // data weight from ProteinProphet
  char lightSeq[500]; // light labeled sequence

  void sort_for_output() {
    for (int ii=dataNum;ii--;) {
      peaks[ii].sort_for_output();
    }
    for (int i=dataNum;i--;) {
      for (int j=i;j--;) {
	if (peaks[i].dataIndx[0]<peaks[j].dataIndx[0]) {
	  dataStrct tmp = peaks[i];
	  peaks[i] = peaks[j];
	  peaks[j] = tmp;
	  int itmp = dataCnts[i];
	  dataCnts[i] = dataCnts[j];
	  dataCnts[j] = itmp;
	}
      }
    }
  }
};


// structure of protein 
typedef struct {
  int indx; // 0 new, 1 analyzed, 2 verified, -1 invalid STATUS
  double ratio[2]; // ratio and error
  double inv_ratio[2]; // inverted ratio and error
  int dataNum; // data number
  seqDataStrct *sequences; // unique sequences
  int *dataCnts; // whether to count the data: 0 no, 1 yes INCLUDE
} proDataStrct;


typedef struct {
  double mean;  // mean
  double merr; // fitting error on mean
  double stddev; // standard deviation of the data
} pValueStrct;


#endif
