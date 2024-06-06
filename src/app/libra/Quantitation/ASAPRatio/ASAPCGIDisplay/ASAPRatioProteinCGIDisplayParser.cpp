/*
Program       : ASAPRaRatioProteinCGIDisplayParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioProteinCGIDisplayParser.cpp 8486 2021-07-01 21:42:28Z dshteyn $


Copyright (C) 2003 Andrew Keller

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

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include "ASAPRatioProteinCGIDisplayParser.h"



ASAPRatioProteinCGIDisplayParser::ASAPRatioProteinCGIDisplayParser(proDataStrct* protein, Array<char*>* xmlfiles, Boolean heavy2light, char* colored_aas) : Parser(NULL) {

  protein_ = protein; // just borrowing it...
  xmlfiles_ = xmlfiles;
  heavy2light_ = heavy2light;
  colored_aas_ = NULL;
  if(colored_aas != NULL) {
    colored_aas_ = new char[strlen(colored_aas)+1];
    strcpy(colored_aas_, colored_aas);
    //cout << "c'tor with colored: " << colored_aas_ << endl;
  }

  databases_ = new Array<char*>;
  basenames_ = new Array<char*>;
  pepproph_timestamps_ = new Array<char*>;
  iproph_timestamps_ = new Array<char*>;
  asap_timestamps_ = new Array<char*>;
  asap_quantHighBGs_ = new Array<Boolean>;
  asap_zeroBGs_ = new Array<Boolean>;
  asap_mzBounds_ = new Array<double>;
  asap_wavelets_ = new Array<bool>;
  elutions_ = new Array<int>;
  aa_modifications_ = new Array<char*>;
  term_modifications_ = new Array<char*>;
  misc_run_conditions_ = new Array<char*>;

  result_index_ = 0;
  pepdata_index_ = 0;
  current_index_ = 0;

  //cout << "ready in c'tor" << endl;

  components_ = new Array<ProDataComponent*>;

  for(int seq = 0; seq < protein->dataNum; seq++)
    for(int pk = 0; pk < protein->sequences[seq].dataNum; pk++)
      for(int data = 0; data < protein->sequences[seq].peaks[pk].dataNum; data++)
	components_->insertAtEnd(new ProDataComponent(seq, pk, data, protein->sequences[seq].peaks[pk].bofIndx, 
						      protein->sequences[seq].peaks[pk].dataIndx[data], heavy2light_, protein->sequences[seq].peaks[pk].msms_run_idx ));


  // now sort by xml_index and data_index
  ProDataComponent** sorted = new ProDataComponent*[components_->length()];
  int k;
  for(k = 0; k < components_->length(); k++)
    sorted[k] = (*components_)[k];

  //qsort(sorted, components_->length(), sizeof(ProDataComponent*), OrderBySeqPkDataInds);
  qsort(sorted, components_->length(), sizeof(ProDataComponent*), OrderByXmlAndDataInds);
  for(k = 0; k < components_->length(); k++)
    components_->replace(k, sorted[k]);


  //DDS:
  //  initMSMSRunIdx();
  //current_index_ = 0;
  //for(int k = 0; k < components_->length(); k++)
  // sorted[k] = (*components_)[k];
  //qsort(sorted, components_->length(), sizeof(ProDataComponent*), OrderBySeqPkDataInds);
  //qsort(sorted, components_->length(), sizeof(ProDataComponent*), OrderByXmlAndDataInds);
  //for(int k = 0; k < components_->length(); k++)
  //components_->replace(k, sorted[k]);


  init(NULL);


  // now sort by seq, pk, and data inds
  sorted = new ProDataComponent*[components_->length()];
  for(k = 0; k < components_->length(); k++)
    sorted[k] = (*components_)[k];

  qsort(sorted, components_->length(), sizeof(ProDataComponent*), OrderBySeqPkDataInds);
  //qsort(sorted, components_->length(), sizeof(ProDataComponent*), OrderByXmlAndDataInds);
  for(k = 0; k < components_->length(); k++) 
    components_->replace(k, sorted[k]);

  
  // verify ordered correctly
  //for(int k = 0; k < components_->length(); k++)
   //  cout << "seq: " << (*components_)[k]->seq_ << " peak: " << (*components_)[k]->peak_ << " data: " << (*components_)[k]->data_ << " DataIdx: " << (*components_)[k]->data_index_<< endl; //<< " msms_runI:" << (*components_)[k]->msms_run_idx_ << " xmlI: " << (*components_)[k]->xml_index_ << " Idx: " << (*components_)[k]->data_index_ << endl;
  

  update(); // compute new ratios based on pepdatastrcts obtained from xml files

  // this is just for test purposes....
  /*
  cout << "<HTML><table cellpadding=\"2\" bgcolor=\"white\">" << endl;
  write(cout);
  cout << "</table></HTML>" << endl;
  */

}


ASAPRatioProteinCGIDisplayParser::~ASAPRatioProteinCGIDisplayParser() {
  cleanup(databases_);
  cleanup(basenames_);
  cleanup(pepproph_timestamps_);  
  cleanup(iproph_timestamps_);
  cleanup(asap_timestamps_);
  delete asap_quantHighBGs_;
  delete asap_zeroBGs_;
  delete asap_wavelets_;
  delete asap_mzBounds_;
  if(components_ != NULL) {
    for(int k = 0; k < components_->length(); k++)
      if((*components_)[k] != NULL)
	delete (*components_)[k];
    delete components_;
  }
  
}

void ASAPRatioProteinCGIDisplayParser::update() {
  // here's where X-J must combine together pepdatastrct info with protein_ include info to 
  // compute new Pk, Seq, and protein_ ratios
  double *ratios, *errors, *inv_ratios, *inv_errors, *weights;
  int *outliers, *seqIndx;
  int seqNum, tmpNum;
  double cmnErr, ratio, error, inv_ratio, inv_error;
  int num[2] = {0, 0};

  int i, j;

  if (protein_->indx == -1) {
    protein_->ratio[0] = -2.;
    protein_->ratio[1] = 0.;
    protein_->inv_ratio[0] = -2.;
    protein_->inv_ratio[1] = 0.;
    return;
  }

  // collect all sequences
  seqIndx = (int *) calloc(protein_->dataNum, sizeof(int));
  seqNum = 0;
  for (i = 0; i < protein_->dataNum; ++i) {
    updateSeqStrctRatio(i);
    //ASAPRatio_getSeqDataStrctRatio(&(data->sequences[i]), pepBofFiles);
    // double check dataCnts
    if(protein_->dataCnts[i] == 1) {
      protein_->dataCnts[i] = 0;
      for (j = 0; j < protein_->sequences[i].dataNum; ++j){
	if(protein_->sequences[i].dataCnts[j] == 1){
	  protein_->dataCnts[i] = 1;
	  break;
	}
      }
    }
    if(protein_->dataCnts[i] == 1 && protein_->sequences[i].indx != -1) {
      seqIndx[i] = 1;
      ++seqNum;
    }
    else
      seqIndx[i] = 0;
  }

  // only 0 or 1 valid sequences
  if(seqNum < 2) {
    if(seqNum < 1) {
      protein_->ratio[0] = -2.;
      protein_->ratio[1] = 0.;
      protein_->inv_ratio[0] = -2.;
      protein_->inv_ratio[1] = 0.;
    } // if(seqNum < 1) {
    else {
      for (i = 0; i < protein_->dataNum; ++i) {
	if(seqIndx[i] == 1) {
	  protein_->ratio[0] = protein_->sequences[i].ratio[0];
	  protein_->ratio[1] = protein_->sequences[i].ratio[1];
	  protein_->inv_ratio[0] = protein_->sequences[i].inv_ratio[0];
	  protein_->inv_ratio[1] = protein_->sequences[i].inv_ratio[1];
	  break;
	}
      }
    } // else {

    if(protein_->indx == 0)
      protein_->indx = 1;
    
    free(seqIndx);

    return;
  } //if(seqNum < 2) {


  // allocate memory
  ratios = (double *) calloc(seqNum, sizeof(double));
  errors = (double *) calloc(seqNum, sizeof(double));
  inv_ratios = (double *) calloc(seqNum, sizeof(double));
  inv_errors = (double *) calloc(seqNum, sizeof(double));
  weights = (double *) calloc(seqNum, sizeof(double));
  outliers = (int *) calloc(seqNum, sizeof(int));

  // get data
  tmpNum = 0;
  for(i = 0; i < protein_->dataNum; ++i) {
    if(seqIndx[i] == 1) {
      ratios[tmpNum] = protein_->sequences[i].ratio[0];
      errors[tmpNum] = protein_->sequences[i].ratio[1];
      inv_ratios[tmpNum] = protein_->sequences[i].inv_ratio[0];
      inv_errors[tmpNum] = protein_->sequences[i].inv_ratio[1];
      ++tmpNum;
    }
  }

  // for(int k = 0; k < tmpNum; k++)
  //  cout << k << ": " << ratios[k] << " +- " << errors[k] << endl;

    
  // check valid error
  cmnErr = 0.;
  tmpNum = 0;
  for(i = 0; i < seqNum; ++i) {
    if(ratios[i] == 0.) // 0
      ++num[0];
    else if(ratios[i] == -1.) // -1
      ++num[1];
    if(errors[i] > 0.) { // valid error
      ++tmpNum;
      cmnErr += errors[i];
    }
  }// for(i = 0; i < seqNum; ++i) {
  
  if(tmpNum < 1) { // all ratios are 0 or -1 or -2
    if(num[0] > num[1]) {// more 0 than -1
      ratio = 0.;
      inv_ratio = 999.;
    }
    else if(num[0] < num[1]) {// more -1 than 0
      ratio = -1.;
      inv_ratio = -1.;
    }
    else { // same 0 and -1
      ratio = -2.;
      inv_ratio = -2.;
    }
    error = 0.;
    inv_error = 0.;

    // get outliers
    for (i = 0; i < seqNum; ++i) {
      if(ratios[i] == ratio)
	outliers[i] = 0;
      else
	outliers[i] = 1;
    }
    
    // store ratios, etc
    protein_->ratio[0] = ratio;
    protein_->ratio[1] = error;

    protein_->inv_ratio[0] = inv_ratio;
    protein_->inv_ratio[1] = inv_error;

    if(protein_->indx == 0) {
      tmpNum = 0;
      for(i = 0; i < protein_->dataNum; ++i) {
	if(seqIndx[i] == 1) {
	  protein_->dataCnts[i] = 1 - outliers[tmpNum];
	  ++tmpNum;
	}
      }
      protein_->indx = 1;
    }
    
    free(seqIndx);
    free(ratios);
    free(errors);    
    free(inv_ratios);
    free(inv_errors);
    free(weights);
    free(outliers);

    return;
  }
  else // get ave err
    cmnErr /= tmpNum;


  // which one is correct????????

  /*
original here
  // get weight
  for(i = 0; i < seqNum; ++i) {
    if(errors[i] <= 0.)
      weights[i] = 1./cmnErr/cmnErr;
    else if(errors[i] < cmnErr/4.)
      weights[i] = 16./cmnErr/cmnErr;
    else
      weights[i] = 1./errors[i]/errors[i];
  }// for(i = 0; i < currPepNum; ++i) {
  */
  
  // get weight
  for(i = 0; i < seqNum; ++i) {
    if(errors[i] < cmnErr)
      weights[i] = 1./cmnErr/cmnErr;
    else
      weights[i] = 1./errors[i]/errors[i];
  }


  //for(int k = 0; k < seqNum; k++)
  //  cout << "wts etc: " << k << ": " << weights[k] << " outliers:  " << outliers[k] << endl;


  // calculate ratio and error
  if(protein_->indx == 0) {
    getDataRatio(&(protein_->ratio[0]), &(protein_->ratio[1]), &(protein_->inv_ratio[0]), &(protein_->inv_ratio[1]),
		 _ASAPRATIO_CONFL_, ratios, errors, inv_ratios, inv_errors, weights, outliers, seqNum, 1);
    tmpNum = 0;
    for(i = 0; i < protein_->dataNum; ++i) {
      if(seqIndx[i] == 1) {
	protein_->dataCnts[i] = 1 - outliers[tmpNum];
	++tmpNum;
      }
    }      
  }
  else
    getDataRatio(&(protein_->ratio[0]), &(protein_->ratio[1]), &(protein_->inv_ratio[0]), &(protein_->inv_ratio[1]),
		 _ASAPRATIO_CONFL_, ratios, errors, inv_ratios, inv_errors, weights, outliers, seqNum, 0);

  // reset indx
  if(protein_->indx == 0)
    protein_->indx = 1;

  // free memory
  free(seqIndx);
  free(ratios);
  free(errors);
  free(inv_ratios);
  free(inv_errors);
  free(weights);
  free(outliers);

  return;

}

void ASAPRatioProteinCGIDisplayParser::updateSeqStrctRatio(int seq) {
  double *ratios, *errors, *inv_ratios, *inv_errors, *weights;
  int *outliers, *peakIndx;
  int peakNum, tmpNum;

  int i, j;

  seqDataStrct* data = &protein_->sequences[seq];

  if (data->indx == -1) {
    data->ratio[0] = -2.;
    data->ratio[0] = 0.;
    data->inv_ratio[0] = -2.;
    data->inv_ratio[1] = 0.;
    return;
  }

  // collect all peaks
  peakIndx = (int *) calloc(data->dataNum, sizeof(int));
  peakNum = 0;
  for (i = 0; i < data->dataNum; ++i) {
    updatePeakStrctRatio(seq, i);
    //ASAPRatio_getDataStrctRatio(&(data->peaks[i]), pepBofFiles[data->peaks[i].bofIndx]);
    // double check dataCount
    if(data->dataCnts[i] == 1) {
      data->dataCnts[i] = 0;
      for (j = 0; j < data->peaks[i].dataNum; ++j){
	if(data->peaks[i].dataCnts[j] == 1){
	  data->dataCnts[i] = 1;
	  break;
	}
      }
    }
    if(data->dataCnts[i] == 1 && data->peaks[i].indx != -1) {
      peakIndx[i] = 1;
      ++peakNum;
    }
    else
      peakIndx[i] = 0;
  }

  // only 0 or 1 valid peaks
  if(peakNum < 2) {
    if(peakNum < 1) {
      data->ratio[0] = -2.;
      data->ratio[1] = 0.;
      data->inv_ratio[0] = -2.;
      data->inv_ratio[1] = 0.;
    } // if(peakNum < 1) {
    else {
      for (i = 0; i < data->dataNum; ++i) {
	if(peakIndx[i] == 1) {
	  data->ratio[0] = data->peaks[i].ratio[0];
	  data->ratio[1] = data->peaks[i].ratio[1];
	  data->inv_ratio[0] = data->peaks[i].inv_ratio[0];
	  data->inv_ratio[1] = data->peaks[i].inv_ratio[1];
	  break;
	}
      }
    } // else {

    if(data->indx == 0)
      data->indx = 1;
    
    free(peakIndx);

    return;
  } //if(peakNum < 2) {


  // allocate memory
  ratios     = (double *) calloc(peakNum, sizeof(double));
  errors     = (double *) calloc(peakNum, sizeof(double));
  inv_ratios = (double *) calloc(peakNum, sizeof(double));
  inv_errors = (double *) calloc(peakNum, sizeof(double));
  weights    = (double *) calloc(peakNum, sizeof(double));
  outliers   = (int *)    calloc(peakNum, sizeof(int));

  // get data
  tmpNum = 0;
  for(i = 0; i < data->dataNum; ++i) {
    if(peakIndx[i] == 1) {
      ratios[tmpNum] = data->peaks[i].ratio[0];
      errors[tmpNum] = data->peaks[i].ratio[1];
      inv_ratios[tmpNum] = data->peaks[i].inv_ratio[0];
      inv_errors[tmpNum] = data->peaks[i].inv_ratio[1];
      weights[tmpNum] = data->peaks[i].weight;
      ++tmpNum;
    }
  }
    
  // calculate ratio and error
  if(data->indx == 0) {
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]),
		 _ASAPRATIO_CONFL_, ratios, errors, inv_ratios, inv_errors, weights, outliers, peakNum, 1);
    tmpNum = 0;
    for(i = 0; i < data->dataNum; ++i) {
      if(peakIndx[i] == 1) {
	data->dataCnts[i] = 1 - outliers[tmpNum];
	++tmpNum;
      }
    }      
  }
  else
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]),
		 _ASAPRATIO_CONFL_, ratios, errors, inv_ratios, inv_errors, weights, outliers, peakNum, 0);
  
  // reset indx
  if(data->indx == 0)
    data->indx = 1;

  // free memory
  free(peakIndx);
  free(ratios);
  free(errors);
  free(inv_ratios);
  free(inv_errors);
  free(weights);
  free(outliers);

  return;
}

void ASAPRatioProteinCGIDisplayParser::updatePeakStrctRatio(int seq, int pk) {
  dataStrct *data = &(protein_->sequences[seq].peaks[pk]);
  double *ratios, *errors, *inv_ratios, *inv_errors, *weights;
  int *outliers, *pepIndx;
  int pepNum, tmpNum;
  double tol = 1.e-2;
  int i;
  pepDataStrct *peptides;

  if (data->indx == -1) {
    data->ratio[0] = -2.;
    data->ratio[1] = 0.;
    data->inv_ratio[0] = -2.;
    data->inv_ratio[1] = 0.;
    return;
  }

  // get peptides
  peptides = (pepDataStrct *) calloc(data->dataNum, sizeof(pepDataStrct));
  pepIndx = (int *) calloc(data->dataNum, sizeof(int));
  pepNum = 0;
  for(i = 0; i < data->dataNum; i++) {
    if(data->dataCnts[i] == 1) { 
      peptides[pepNum] = getPepDataStrct(seq, pk, i);
      if(peptides[pepNum].indx != -1) {
	pepIndx[i] = 1;
	++pepNum;
      }
      else 
	pepIndx[i] = 0;
    }
    else {
      pepIndx[i] = 0;
    }

  } // next one

  // calculate ratio

  // only 0 or 1 valid peptides
  if(pepNum < 2) {
    if(pepNum < 1) {
      data->ratio[0] = -2.;
      data->ratio[1] = 0.;
      data->inv_ratio[0] = -2.;
      data->inv_ratio[1] = 0.;
      data->weight = 0.;
    }
    else {
      data->ratio[0] = peptides[0].pepRatio[0];
      data->ratio[1] = peptides[0].pepRatio[1];
      data->inv_ratio[0] = peptides[0].pepH2LRatio[0];
      data->inv_ratio[1] = peptides[0].pepH2LRatio[1];
      data->weight = peptides[0].pepArea;
    } // else {
    if(data->indx == 0) {
      for (i = 0; i < data->dataNum; ++i) {
	if(data->dataCnts[i] == 1 && pepIndx[i] == 0)  
	  data->dataCnts[i] = 0;
      }
      data->indx = 1;
    }

    free(pepIndx);
    free(peptides);

    return;
  }

  // allocate memory
  ratios = (double *) calloc(pepNum, sizeof(double));
  errors = (double *) calloc(pepNum, sizeof(double));
  inv_ratios = (double *) calloc(pepNum, sizeof(double));
  inv_errors = (double *) calloc(pepNum, sizeof(double));
  weights = (double *) calloc(pepNum, sizeof(double));
  outliers = (int *) calloc(pepNum, sizeof(int));

  // get data
  for(i = 0; i < pepNum; ++i) {
    ratios[i] = peptides[i].pepRatio[0];
    errors[i] = peptides[i].pepRatio[1];    
    inv_ratios[i] = peptides[i].pepH2LRatio[0];
    inv_errors[i] = peptides[i].pepH2LRatio[1];
    weights[i] = peptides[i].pepArea;
  }
    
  // calculate ratio and error
  if(data->indx == 0) {
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]),
		 _ASAPRATIO_CONFL_, ratios, errors, inv_ratios, inv_errors, weights, outliers, pepNum, 1);
    tmpNum = 0;
    for (i = 0; i < data->dataNum; ++i) {
      if(pepIndx[i] == 1) {
	data->dataCnts[i] = 1 - outliers[tmpNum];
	++tmpNum;
      }
    }      
  }
  else
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]),
		 _ASAPRATIO_CONFL_, ratios, errors, inv_ratios, inv_errors, weights, outliers, pepNum, 0);

  // weight
  data->weight = 0.;
  for (i = 0; i < pepNum; ++i) {
    if(outliers[i] == 0 && weights[i] > data->weight)
      data->weight = weights[i];
  }

  // if identical ratios
  tmpNum = 0;
  for (i = 0; i < pepNum; ++i) {
    if(outliers[i] == 0 
       && fabs(ratios[i]-data->ratio[0]) > tol*data->ratio[0]) {
      ++tmpNum;
    }
  }
  if(tmpNum == 0) {
    data->ratio[1] = 0.;
    data->inv_ratio[1] = 0.;
    tmpNum = 0;
    for (i = 0; i < pepNum; ++i) {
      if(outliers[i] == 0) {
	data->ratio[1] += errors[i];
	data->inv_ratio[1] += inv_errors[i];
	++tmpNum;
      }
    }
    if(tmpNum > 0) {
      data->ratio[1] /= tmpNum;
      data->inv_ratio[1] /= tmpNum;
    }
    else {
      for (i = 0; i < pepNum; ++i) {
	if(outliers[i] == 0) {
	  data->ratio[1] = errors[i];
	  data->inv_ratio[1] = inv_errors[i];
	  break;
	}
      }
    }
  } // if(tmpNum == 0) {
  
  // reset indx
  if(data->indx == 0)
    data->indx = 1;

  // free memory
  free(peptides);
  free(pepIndx);
  free(ratios);
  free(errors); 
  free(inv_ratios);
  free(inv_errors);
  free(weights);
  free(outliers);

  return;
}


// For a given set of data, dataErrs, and dataWghs, this function identifies 
// any outliers and gets the mean and confidence interval of ratio.
void ASAPRatioProteinCGIDisplayParser::getDataRatio(double *ratio, double *error,  double *inv_ratio, double *inv_error, double confL, 
		  double *data, double *dataErrs, double *inv_data, double *inv_dataErrs,
		  double *dataWghs, int *dataIndx, int dataSize,
		  int testType)
{
  
  //  void findMeanAndStdDevWeight(double *mean, double *error,
  //			       double *data, double *weight, int size);
  //  void DixonTest(double *data, int *outliers, int size);
  
  int counts[4] = {0, 0, 0, 0};
  double *dtRatios, *dtErrors, *dt_inv_Ratios, *dt_inv_Errors, *dtWeights;
  int *dtIndx;
  int pass;
  int count, vldNum;
  double sum, inv_sum;
  double tmpError, tmp_inv_Error;
  double acc = 0.01;
  int i, j;
  int exclude_olutliers = 1;


  // check whether there are valid data
  for(i = 0; i < dataSize; ++i) {
    if(data[i] == 0.) 
      ++counts[0];
    else if(data[i] == -1.) 
      ++counts[1];
    else if(data[i] == -2.) 
      ++counts[2];
    else
      ++counts[3];
  }

  // easy output for invalid data set
  if(counts[3] < 1) {
    if(counts[0] > counts[1]) {
      *ratio = 0.;
      *error = 0.;
      *inv_ratio = 999.;
      *inv_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != 0.)
	    dataIndx[i] = exclude_olutliers;
	  else
	    dataIndx[i] = 0;
      }
    }
    else if(counts[0] < counts[1]) {
      *ratio = -1.;
      *error = 0.;      
      *inv_ratio = -1.;
      *inv_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != -1.)
	    dataIndx[i] = exclude_olutliers;
	  else
	    dataIndx[i] = 0;
      }
    }
    else {
      *ratio = -2.;
      *error = 0.;
      *inv_ratio = -2.;
      *inv_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != -2.)
	    dataIndx[i] = exclude_olutliers;
	  else
	    dataIndx[i] = 0;
      }
    }
    return;
  } // if(count[3] < 1) {

  // allocate memory
  dtRatios = (double *) calloc(dataSize, sizeof(double));
  dtErrors = (double *) calloc(dataSize, sizeof(double));
  dt_inv_Ratios = (double *) calloc(dataSize, sizeof(double));
  dt_inv_Errors = (double *) calloc(dataSize, sizeof(double));
  dtWeights = (double *) calloc(dataSize, sizeof(double));
  dtIndx = (int *) calloc(dataSize, sizeof(int));


  // identify outliers

  // collect valid data, transform into log(ratio)
  for (i = 0; i < dataSize; ++i) 
    if(data[i] > 0.)
      dataIndx[i] = 0;
    else
      dataIndx[i] = exclude_olutliers;

  count = 0;
  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) {
      pass = 1;
      for (j = 0; pass == 1 && j < count; ++j) {
	if(fabs(log(data[i])-dtRatios[j]) < acc*dtRatios[j]){
	  pass = 0;
	}
      }
      if(pass == 1) {
	dtRatios[count] = log(data[i]);
	dt_inv_Ratios[count] = log(inv_data[i]);
	++count;
      }
    }
  }

  // identify any outliers
  if(testType != 0)
    DixonTest(dtRatios, dtIndx, count);
  else
    for (i = 0; i < count; ++i)
      dtIndx[i] = 0;

  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) {
      pass = 1;
      for (j = 0; pass == 1 && j < count; ++j) {
	if(dtIndx[j] == 1 
	   && fabs(log(data[i])-dtRatios[j]) < acc*dtRatios[j]){
	  pass = 0;
	}
      }
      if(pass == 0) 
	dataIndx[i] = 1;
    }
  }

  // get ratio and error

  // collect valid date
  count = 0; 
  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) { // not an outlier
      dtRatios[count] = data[i];
      dtErrors[count] = dataErrs[i];
      dt_inv_Ratios[count] = inv_data[i];
      dt_inv_Errors[count] = inv_dataErrs[i];
      dtWeights[count] = dataWghs[i];
      ++count;
    }
  }

  // calculate ratio and error
  if(count < 1) { // no valid data
    *ratio = -2.;
    *error = 0.;    
    *inv_ratio = -2.;
    *inv_error = 0.;
  } //if(count < 1) { // no valid data
  else if(count == 1) { // only one valid data
    *ratio = dtRatios[0];
    *error = dtErrors[0];
    *inv_ratio = dt_inv_Ratios[0];
    *inv_error = dt_inv_Errors[0];
  }//  else if(count == 1) { // only one valid data
  else {
    // transform into log(ratio)
    for (i = 0; i < count; ++i) {
      dtErrors[i] /= dtRatios[i];
      dtRatios[i] = log(dtRatios[i]);
      dt_inv_Errors[i] /= dt_inv_Ratios[i];
      dt_inv_Ratios[i] = log(dt_inv_Ratios[i]);
    }
    // calculate the light:heavy ratio by weight 
    findMeanAndStdDevWeight(ratio, error, dtRatios, inv_ratio, inv_error, dt_inv_Ratios, dtWeights, count);

    sum = 0.;
    inv_sum = 0.;
    vldNum = 0;
    for (i = 0; i < count; ++i) {
      if(dtErrors[i] > 0.) {
	sum += 1./dtErrors[i]/dtErrors[i];
	inv_sum += 1./dt_inv_Errors[i]/dt_inv_Errors[i];
	++vldNum;
      }
    }
    if(vldNum > 0) {
      tmpError = 1./sqrt(sum);
      tmp_inv_Error = 1./sqrt(inv_sum);
    } // if(vldNum > 0) {
    else {
      tmpError = 0.;
      tmp_inv_Error = 0.; 
    }

    *error = sqrt((*error)*(*error)+tmpError*tmpError);
    *inv_error = sqrt((*inv_error)*(*inv_error)+tmp_inv_Error*tmp_inv_Error);
    
    *inv_ratio = exp(*inv_ratio);
    *inv_error *= (*inv_ratio);
    // transform back to ratio
    *ratio = exp(*ratio);
    *error *= (*ratio);
  }//else {
  
  // free memory
  free(dtRatios);
  free(dtErrors);
  free(dt_inv_Ratios);
  free(dt_inv_Errors);
  free(dtWeights);
  free(dtIndx);

  return;
}


// For a set of data and weight, this function finds the mean and standard deviation.
void ASAPRatioProteinCGIDisplayParser::findMeanAndStdDevWeight(double *mean, double *error,
			     double *data, double *inv_mean, double *inv_error,
			     double *inv_data, double *weight, int size)
{
  double sum0, sum1, sum2, inv_sum1, inv_sum2, sumW;
  double nEff;
  double mnValue, mxValue, mn_inv_Value, mx_inv_Value;
  int count, i;
  
  if (size < 2) {
    *mean = *data;
    *error = 0.;
    *inv_mean = *inv_data;
    *inv_error = 0.;
    return;
  }

  // ensure weight is valid
  count = 0;
  sum0 = 0.;
  mnValue = data[0];
  mxValue = data[0];
  mn_inv_Value = inv_data[0];
  mx_inv_Value = inv_data[0];
  for (i = 0; i < size; ++i) {
    if(weight[i] >= 0.) {
      ++count;
      sum0 += weight[i];
    }
    mnValue = mnValue < data[i] ? mnValue : data[i];
    mxValue = mxValue > data[i] ? mxValue : data[i];

    mn_inv_Value = mn_inv_Value < inv_data[i] ? mn_inv_Value : inv_data[i];
    mx_inv_Value = mx_inv_Value > inv_data[i] ? mx_inv_Value : inv_data[i];
  }

  if(mnValue >= mxValue) {
    *mean = mnValue;
    *error = 0.;
    *inv_mean = mn_inv_Value;
    *inv_error = 0.;
    return;
  }

  if(count < size || sum0 == 0.) {   // no all have valid weight
    if(count < 1 || sum0 == 0.) {   // if no data has weight
      for (i = 0; i < size; ++i) {
	weight[i] = 1.;
      }
    } // if(count < 1 || sum0 == 0.) {   // if no data has weight
    else {
      sum0 /= count;
      for (i = 0; i < size; ++i) {
	if(weight[i] < 0.) {
	  weight[i] = sum0;
	}
      }
    }
  }

  // get mean and std. dev.
  sum0 = 0.;
  sum1 = 0.;
  sum2 = 0.;
  inv_sum1 = 0.;
  inv_sum2 = 0.;
  sumW = 0.;
  for (i = 0; i < size; ++i) {
    sum0 += weight[i];
    sum1 += data[i]*weight[i];
    sum2 += data[i]*data[i]*weight[i];
    inv_sum1 += inv_data[i]*weight[i];
    inv_sum2 += inv_data[i]*inv_data[i]*weight[i];
    sumW += weight[i]*weight[i];
  }

  // get mean
  *mean = sum1/sum0;
  *inv_mean = inv_sum1/sum0;
  
  // get std. dev.
  if(sum2*sum0-sum1*sum1 > 0.) {
    nEff = sum0*sum0/sumW;
    if(nEff > 2.) {
      *error = sqrt((sum2*sum0-sum1*sum1)*nEff/(nEff-1.))/sum0;
      *inv_error = sqrt((inv_sum2*sum0-inv_sum1*inv_sum1)*nEff/(nEff-1.))/sum0;
    }
    else {
      *error = sqrt(2.*(sum2*sum0-sum1*sum1))/sum0;
      *inv_error = sqrt(2.*(inv_sum2*sum0-inv_sum1*inv_sum1))/sum0;
    }
  }
  else {
     *error = 0.;
     *inv_error = 0.;
   }
    
  return;
}


// This function uses Dixon's test with alpha = 0.05 to identify any outliers.
void ASAPRatioProteinCGIDisplayParser::DixonTest(double *data, int *outliers, int size)
{
  //    double PadeApprx(double x, double *xa, double *ya, int size);

  // cutoff values in Dixon's test: n = 3, ..., 30, INF.
  double ya[29] = {0.941, 0.765, 0.642, 0.560, 0.507, 0.554,
		   0.512, 0.477, 0.576, 0.546, 0.521, 0.546,
		   0.525, 0.507, 0.490, 0.475, 0.462, 0.450,
		   0.440, 0.430, 0.421, 0.413, 0.406, 0.399,
		   0.393, 0.387, 0.381, 0.376, 0.};

  // values of 1/n: 1/3, ..., 1/30, 1/INF.
  double xa[29] = {0.333333, 0.250000, 0.200000, 0.166667, 0.142857, 
		   0.125000, 0.111111, 0.100000, 0.090909, 0.083333, 
		   0.076923, 0.071429, 0.066667, 0.062500, 0.058824, 
		   0.055556, 0.052632, 0.050000, 0.047619, 0.045455, 
		   0.043478, 0.041667, 0.040000, 0.038462, 0.037037, 
		   0.035714, 0.034483, 0.033333, 0.};
  
  int cnstSize = 29;
  int *dataIndx;
  int startIndx, endIndx;
  int count;
  double ratio1, ratio2;
  double cutoff, x;
  int i, j;

  // assume none is an outlier
  for (i = 0; i < size; ++i)
    outliers[i] = 0;
  if (size < 3) // not enough data for checking
    return;
  
  // get dataIndx for ordered data
  dataIndx = (int *) calloc(size, sizeof(int));
  for (i = 0; i < size; ++i)
    dataIndx[i] = i;
  for(i = 0; i < size; ++i) {
    for(j = 0; j < size-i-1; ++j) {
      if(data[dataIndx[j]] > data[dataIndx[j+1]]) {
	count = dataIndx[j];
	dataIndx[j] = dataIndx[j+1];
	dataIndx[j+1] = count;
      }
    }
  }
  
  // check for outliers
  count = 0;
  startIndx = 0;
  endIndx = size;
  while(size > 2  // enough data for checking
	&& count != size // look for more when an outlier is identified
	&& data[dataIndx[startIndx]] != data[dataIndx[endIndx-1]]) { 

    // restore size
    count = size;

    // get cutoff
    if (size < 3)
      cutoff = 1.;
    else if(size <= cnstSize+1)
      cutoff = ya[size-3];
    else {
      x = 1./((double) size);
      cutoff = PadeApprx(x, xa, ya, cnstSize);
    }

    // get ratio
    if(size < 8) {
      ratio1 = (data[dataIndx[startIndx+1]]-data[dataIndx[startIndx]])
	/(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx]]);
      ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-2]])
	/(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx]]);
    }
    else if(size < 11) {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-2]]) { 
	ratio1 = (data[dataIndx[startIndx+1]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-2]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+1]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-2]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+1]]);
      }
      else
	ratio2 = 0.;  
    }
    else if(size < 14) {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-2]]) { 
	ratio1 = (data[dataIndx[startIndx+2]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-2]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+1]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-3]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+1]]);
      }
      else
	ratio2 = 0.;  
    }
    else {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-3]]) { 
	ratio1 = (data[dataIndx[startIndx+2]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-3]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+2]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-3]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+2]]);
      }
      else
	ratio2 = 0.;  
    }

    // check ratio
    if(ratio1 > ratio2) {
      if(ratio1 > cutoff) { // an outlier
	outliers[dataIndx[startIndx]] = 1;
	--size;
	++startIndx;
      }
    } // if(ratio1 > ratio2) {
    else {
      if(ratio2 > cutoff) { // an outlier
	outliers[dataIndx[endIndx-1]] = 1;
	--size;
	--endIndx;
      }
    } //else {
  } // while(size > 2  // enough data for checking  

  free(dataIndx);
  return;
}


// This function returns the value of Pade Approximation.
double ASAPRatioProteinCGIDisplayParser::PadeApprx(double x, double *xa, double *ya, int size)
{
  double y, dy;
  double tiny = 1.e-25;
  int m,i,ns=1;
  double w,t,hh,h,dd,*c,*d;
  double *xb, *yb;
  int n = size;

  // convert into 1 ... n
  xb = (double *) calloc(n+1, sizeof(double));
  yb = (double *) calloc(n+1, sizeof(double));
  for (i = 0; i < size; ++i) {
    xb[i+1] = xa[i];
    yb[i+1] = ya[i];
  }

  // use ratint
  c = (double *) calloc(n+1, sizeof(double));
  d = (double *) calloc(n+1, sizeof(double));

  hh=fabs(x-xb[1]);
  for (i=1;i<=n;i++) {
    h=fabs(x-xb[i]);
    if (h == 0.0) {
      y=yb[i];
      dy=0.0;
      free(c);
      free(d);
      free(xb);
      free(yb);
      return y;
    } else if (h < hh) {
      ns=i;
      hh=h;
    }
    c[i]=yb[i];
    d[i]=yb[i]+tiny;
  }
  y=yb[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      w=c[i+1]-d[i];
      h=xb[i+m]-x;
      t=(xb[i]-x)*d[i]/h;
      dd=t-c[i+1];
      if (dd == 0.0) {
	printf("Error in routine PadeApprx\n");
	free(c);
	free(d);
	free(xb);
	free(yb);
	return y;
      }
      dd=w/dd;
      d[i]=c[i+1]*dd;
      c[i]=t*dd;
    }
    y += (dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free(c);
  free(d);
  free(xb);
  free(yb);

  return y;
}


void ASAPRatioProteinCGIDisplayParser::write(ostream& os) {
  const int max_base_len = 127; // truncate 'experiment' titles over this length
  char jscript_tag[50000];
  char radio_tag[5000];
  char asap_tag[5000];
  char tmp_buffer[200];

  os << "<span class='sidebuttons'>" << endl;
  os << "<input type='button' onclick='showallpeps(true,\"allpepstable\");' value='Expand All'><br/><br/>" << endl;
  os << "<input type='button' onclick='showallpeps(false,\"allpepstable\");' value='Collapse All'>" << endl;
  os << "</span>" << endl;

  os << "<table id='allpepstable' class='tppsimpletable'>" << endl;
  os << "<tr>" << endl;
  os << "<th colspan='4' class='head'>Peptide</th>" << endl;
  os << "<th colspan='1' class='head'>Weight</th>" << endl;
  os << "<th colspan='2' class='head'>Ratio(s)</th>" << endl;
  os << "<th colspan='4' class='head'>Acceptance</th>" << endl;
  os << "</tr>" << endl;

  for(int seq = 0; seq < protein_->dataNum; seq++) {
    setSeqRadioTag(radio_tag, seq, protein_->dataCnts[seq]);
    setASAPRatioTag(asap_tag, protein_->sequences[seq].ratio[0], protein_->sequences[seq].ratio[1]);

    os << "<tr class='upep tpphov' onclick='show(\"peptide"<<(seq+1)<<"\",this)'>" << endl;
    os << "<th style='text-align:left;' colspan='4' class='subhead'>&nbsp;" << (seq+1) << ". "
       << "<span class='tpporange'>" << protein_->sequences[seq].lightSeq << "</span>"
       << "<span title='ntt=" << (*components_)[getResultIndex(seq, 0, 0, False)]->result_->num_tol_term_ << "' style='margin-left:20px;float:right;'>";

    if ( (*components_)[getResultIndex(seq, 0, 0, False)]->result_->num_tol_term_ == 2)
      os << "<span class='ntt on'></span>";
    else
      os << "<span class='ntt off'></span>";

    if ( (*components_)[getResultIndex(seq, 0, 0, False)]->result_->num_tol_term_ == 0)
      os << "<span class='ntt off'></span>";
    else
      os << "<span class='ntt on'></span>";

    os << "</span>"
       << "</th>" << endl;

    sprintf(tmp_buffer, "%0.2f", protein_->sequences[seq].weight);
    os << "<th style='text-align:left;' class='subhead'>"
       << "<span style='max-width:100px;' class='tpp_loadprogress'>"
       << "<span class='tpp_loadprogressbar' style='font-size:small;font-weight:normal;width:" << (protein_->sequences[seq].weight*100) << "%'>" << tmp_buffer << "</span>"
       << "</span>" << endl;
    if ( (*components_)[getResultIndex(seq, 0, 0, False)]->result_->proteins_->size() == 1) // unique mapping?
      os << "&#x25b2;";
    os << "</th>" << endl;

    os << "<th class='subhead'>&nbsp;</th>" << endl;

    os << "<th class='subhead asapratio'>" << asap_tag << "</th>" << endl;

    // borrowed from PepXMLViewer:
    int intratio = int(50+33*log(0.00000001 + protein_->sequences[seq].ratio[0])); // log scale for -4.5 < log(r) < +4.5
    string qclass =
      (protein_->sequences[seq].ratio[0] < 0) ? "quantNC" :
      (intratio   < 0) ? "quant0"   :
      (intratio > 100) ? "quant100" :
      "quant" + std::to_string(intratio);

    string style =
      (protein_->sequences[seq].ratio[0] < 0) ? "display:none;" :
      (intratio   < 0) ? "position:relative;left:-26px;" :
      (intratio > 100) ? "position:relative;left:26px;"  :
      "position:relative;left:" + std::to_string(intratio/2 -25) + "px;";

    os << "<th colspan='2' class='" << qclass << " subhead'><span style='" << style << "'>&#x2022;</span></th>" << endl;
    //os << "<th class='subhead'>&nbsp;</th>" << endl;
    //os << "<th class='subhead'>&nbsp;</th>" << endl;

    os << "<th class='subhead' class='ratioynm' style='background-color:";
    if (protein_->dataCnts[seq] == 1)
      os << "#207f0e;'>&check;";
    else if (protein_->dataCnts[seq] == 0)
      os << "#f00;'>&cross;";
    else
      os << "#5e6a71;'>&cross;";
    os << "</th>" << endl;

    os << "<th class='subhead'>" << radio_tag << "</th>\n" << endl;
    os << "</tr>" << endl;


    string disp = (protein_->dataNum < 3) ? "table-row" : "none";

    for(int pk = 0; pk < protein_->sequences[seq].dataNum; pk++) {
      for(int data = 0; data < protein_->sequences[seq].peaks[pk].dataNum; data++) {
	if(data == 0) {
	  setPeakRadioTag(radio_tag, seq, pk, protein_->sequences[seq].dataCnts[pk]);
	  setASAPRatioTag(asap_tag, protein_->sequences[seq].peaks[pk].ratio[0], protein_->sequences[seq].peaks[pk].ratio[1]);

	  // shorten name if too long
	  if(strlen((*basenames_)[(*components_)[getResultIndex(seq, pk, data, False)]->basename_index_]) > max_base_len) {
	    strcpy(tmp_buffer, "(...)" );
	    strcat(tmp_buffer,
		   (*basenames_)[(*components_)[getResultIndex(seq, pk, data, False)]->basename_index_] +
		   (strlen((*basenames_)[(*components_)[getResultIndex(seq, pk, data, False)]->basename_index_])-max_base_len)
		   );
 	  }
	  else {
	    strcpy(tmp_buffer, (*basenames_)[(*components_)[getResultIndex(seq, pk, data, False)]->basename_index_]);
	  }

	  os << "<tr style='display:" << disp << ";border-top:2px solid #002664;' data-peptide='peptide" << (seq+1) << "'>" << endl;

	  os << "<td style='font-weight:bold;' class='name'>&#x25b6; ";
	  os << "<a title='view peptide results in PepXMLViewer' target='pepxml' href='"<< getCgiUrl() << "PepXMLViewer.cgi" 
	     << "?xmlFileName=" << (*basenames_)[(*components_)[getResultIndex(seq, pk, data, False)]->basename_index_]
	     << "&requiredPeptideText=^" << protein_->sequences[seq].lightSeq << "$"
	     << "&sortField=Gcalc_neutral_pep_mass"
	     << "'>" << (*components_)[getResultIndex(seq, pk, data, False)]->scan_ << "</a></td>" << endl;

	  os << "<td class='name' colspan='5'>Experiment: <span style='font-weight:bold;'>" << tmp_buffer  << "</span></td>" << endl;
	  os << "<td class='name' style='text-align:right'>" << asap_tag << "</td>" << endl;
	  os << "<td class='name'>&nbsp;</td>" << endl;

	  os << "<td class='ratioynm' style='background-color:";
	  if (protein_->sequences[seq].dataCnts[pk] == 1)
	    os << "#207f0e;'>&check;";
	  else if (protein_->sequences[seq].dataCnts[pk] == 0)
	    os << "#f00;'>&cross;";
	  else
	    os << "#5e6a71;'>&cross;";
	  os << "</td>" << endl;

	  os << "<td class='name'>&nbsp;</td>" << endl;
	  os << "<td class='name'>" << radio_tag << "</td>" << endl;
	  os << "</tr>" << endl;
	}

	setDataRadioTag(radio_tag, seq, pk, data, protein_->sequences[seq].peaks[pk].dataCnts[data]);

	os << "<tr style='display:" << disp << ";' class='tpphov' data-peptide='peptide" << (seq+1) << "'>" << endl;
	//os << "<td style='background-color:#fff;' width=\"10\">&#x25b6;</td>" << endl;


#ifdef USE_STD_MODS
	(*components_)[getResultIndex(seq, pk, data, True)]->write(os, xmlfiles_, basenames_, pepproph_timestamps_, iproph_timestamps_, databases_, asap_timestamps_, asap_quantHighBGs_, asap_zeroBGs_, asap_mzBounds_, asap_wavelets_, radio_tag, misc_run_conditions_, colored_aas_);
#endif
#ifndef USE_STD_MODS
	(*components_)[getResultIndex(seq, pk, data, True)]->write(os, xmlfiles_, basenames_, pepproph_timestamps_, iproph_timestamps_, databases_, asap_timestamps_, asap_quantHighBGs_, asap_zeroBGs_, asap_mzBounds_, asap_wavelets_, radio_tag, aa_modifications_, term_modifications_, misc_run_conditions_);
#endif

        os << "<td class='ratioynm' style='background-color:";
        if (protein_->sequences[seq].peaks[pk].dataCnts[data] == 1)
          os << "#207f0e;'>&check;";
        else if (protein_->sequences[seq].peaks[pk].dataCnts[data] == 0)
          os << "#f00;'>&cross;";
        else
          os << "#5e6a71;'>&cross;";
        os << "</td>" << endl;

	os << "<td>&nbsp;</td>" << endl;
	os << "<td>&nbsp;</td>" << endl;

	os << "<td>" << radio_tag << "</td>" << endl;
	os << "</tr>" << endl;
      } // next data

    } // next pk

  } // next seq

  os << "</table>\n<br>\n\n";
}


void ASAPRatioProteinCGIDisplayParser::setASAPRatioTag(char* text, double mean, double error) {
  // check for heavy2light
  if(heavy2light_) {
    if(mean == 0.0) {
      mean = 999.;
      error = 0.0;
    }
    else if(mean >= 999.) {
      mean = 0.0;
      error = 0.0;
    }
    else {
      error /= mean * mean;
      mean = 1.0 / mean;
    }
  }

  double percent_error;
  if(mean == 0.0)
    percent_error = 0.0;
  else 
    percent_error = error * 100 / mean;

  if(percent_error >= 50.0)
    sprintf(text, "%0.2f &plusmn; %0.2f <font color=\"red\">(%0.0f%%)</font>", mean, error, percent_error);
  else if(percent_error < 10.0)
    sprintf(text, "%0.2f &plusmn; %0.2f (%0.1f%%)", mean, error, percent_error);
  else
    sprintf(text, "%0.2f &plusmn; %0.2f <span class='tpporange'>(%0.0f%%)</span>", mean, error, percent_error);
}


void ASAPRatioProteinCGIDisplayParser::makeRadioTag(char* tag, int value, char* html_name) {
  // these are 'select' options
  if(value == 1)
    sprintf(tag, "\n<select onchange='this.style.backgroundColor=\"#ff5f00\";' onclick='event.stopPropagation();' name='%s'>\n<option value='1' selected='yes'/>YES</option>\n<option value='0'/>No</option>\n<option value='-1'/>NA</option>\n</select>\n", html_name);
  else if(value == 0)
    sprintf(tag, "\n<select onchange='this.style.backgroundColor=\"#ff5f00\";' onclick='event.stopPropagation();' name='%s'>\n<option value='1'/>YES</option>\n<option value='0' selected='yes'/>No</option>\n<option value='-1'/>NA</option>\n</select>\n", html_name);
  else
    sprintf(tag, "\n<select onchange='this.style.backgroundColor=\"#ff5f00\";' onclick='event.stopPropagation();' name='%s'>\n<option value='1'/>YES</option>\n<option value='0'/>No</option>\n<option value='-1' selected='yes'/>NA</option>\n</select>\n", html_name);

}

void ASAPRatioProteinCGIDisplayParser::setDataRadioTag(char* tag, int seq, int pk, int data, int value) {
  char* radio_name = new char[200];
  sprintf(radio_name, "peak_%d_%d_dataCnts_%d", seq, pk, data);
  makeRadioTag(tag, value, radio_name);
}

void ASAPRatioProteinCGIDisplayParser::setPeakRadioTag(char* tag, int seq, int pk, int value) {
  char* radio_name = new char[200];
  sprintf(radio_name, "sequence_%d_dataCnts_%d", seq, pk);
  makeRadioTag(tag, value, radio_name);
}

void ASAPRatioProteinCGIDisplayParser::setSeqRadioTag(char* tag, int seq, int value) {
  char* radio_name = new char[200];
  sprintf(radio_name, "protein_dataCnts_%d", seq);
  makeRadioTag(tag, value, radio_name);
}


int OrderByXmlAndDataInds(void const *a, void const *b) {
  ProDataComponent** dd1 = (ProDataComponent**)a;
  ProDataComponent** dd2 = (ProDataComponent**)b;
  if((*dd1)->xml_index_ < (*dd2)->xml_index_)
    return -1;
  else if((*dd1)->xml_index_ > (*dd2)->xml_index_)
    return 1;
  
  //DDS:
  //if((*dd1)->msms_run_idx_ < (*dd2)->msms_run_idx_)
  //  return -1;
  //else if((*dd1)->msms_run_idx_ > (*dd2)->msms_run_idx_)
  //  return 1;
  
  if((*dd1)->data_index_ < (*dd2)->data_index_)
    return -1;
  else if((*dd1)->data_index_ > (*dd2)->data_index_)
    return 1;
  return 0;
}

int OrderBySeqPkDataInds(void const *a, void const *b) {
  ProDataComponent** dd1 = (ProDataComponent**)a;
  ProDataComponent** dd2 = (ProDataComponent**)b;
  if((*dd1)->seq_ < (*dd2)->seq_)
    return -1;
  else if((*dd1)->seq_ > (*dd2)->seq_)
    return 1;
  
  if((*dd1)->peak_ < (*dd2)->peak_)
    return -1;
  else if((*dd1)->peak_ > (*dd2)->peak_)
    return 1;
  
  if((*dd1)->data_ < (*dd2)->data_)
    return -1;
  else if((*dd1)->data_ > (*dd2)->data_)
    return 1;

  return 0;
}

int ASAPRatioProteinCGIDisplayParser::getTimestampIndex(const char* timestamp) {
  for(int k = 0; k < asap_timestamps_->length(); k++)
    if(! strcmp((*asap_timestamps_)[k], timestamp))
      return k;

  return -1;
}


void ASAPRatioProteinCGIDisplayParser::cleanup(Array<char*>* data) {
  if(data != NULL) {
    for(int k = 0; k < data->length(); k++)
      if((*data)[k] != NULL)
	delete (*data)[k];
    delete data;
  }
}

void ASAPRatioProteinCGIDisplayParser::initMSMSRunIdx() {
  char* data = NULL;
  Tag* tag = NULL;
  int index = -1;
  char *nextline = new char[line_width_];
  cout << "Enter" << endl;
  for(int k = 0; k < xmlfiles_->length(); k++) {
    RACI fin((*xmlfiles_)[k]); // can read gzipped xml
    if(! fin) {
      cout << "ASAPRatioProteinCGIDisplayParser: error opening " << (*xmlfiles_)[k] << endl;
      exit(1);
    }
    int msms_run_idx = -1;
    while(current_index_ < components_->length() && 
	  (*components_)[current_index_]->xml_index_ == k && 
	  fin.getline(nextline, line_width_)) {
     
      data = strstr(nextline, "<");
      if (data != NULL) {
	tag = new Tag(data);
	if(tag != NULL) {
	  if(tag->isStart() && ! strcmp(tag->getName(), "msms_run_summary")) {
	      msms_run_idx++;
	      delete tag;
	  }
	  else if(tag->isStart() && ! strcmp(tag->getName(), "spectrum_query")) {
	    index = atoi(tag->getAttributeValue("index"));
	    delete tag;
	  }
	}
      }
      while (current_index_ < components_->length() &&
	     (*components_)[current_index_]->xml_index_ == k && 
	     (*components_)[current_index_]->data_index_ == index) {
	(*components_)[current_index_++]->msms_run_idx_ = msms_run_idx;
      } 
    }
    fin.close();
  }
  cout << "Exit" << endl;
  delete[] nextline;
}

void ASAPRatioProteinCGIDisplayParser::parse(const char* xmlfile) {
  Tag* tag = NULL;
  //  int line_width = 10000;
  char *nextline = new char[line_width_];
  char* data = NULL;

  Boolean analyze = False;
  Boolean collection = False;

  double probability = -4.0;
  Boolean adjusted_prob = False;
  Boolean incomplete_prob = False;

  Array<Tag*>* tags = NULL;
  long scan=-1;
  int precursor_charge=0;
  //char peptide[200];
  int index = -1;
  Array<Tag*>* search_result_tags = NULL;

  SearchResult* nextresult = NULL;

  prodatacomponent_struct comp_data;

  int asap_index = -1;

  //comp_data.xpressratio = -1.0;
  comp_data.xpressratio[0] = 0;

  comp_data.lightfirstscan = -1;
  comp_data.lightlastscan = -1;
  comp_data.heavyfirstscan = -1;
  comp_data.heavylastscan = -1;
  comp_data.lightmass = -1.0;
  comp_data.heavymass = -1.0;
  comp_data.masstol = -1.0;
  comp_data.xpresslight = -1;


  comp_data.asap_mean = -3.0;
  comp_data.asap_err = -1.0;
  comp_data.asap_inv_mean = -3.0;
  comp_data.asap_inv_err = -1.0;
  comp_data.asapratio_index = -1;
  comp_data.score_summary[0] = 0;
  char search_engine[300];
  search_engine[0] = 0;

  Boolean search_score_found = False;

  char impossible_match[] = "^^^^^^^^^^^^^^^^^";
  char comet_md5_match[] = "parameter name=\"md5_check_sum\"";
  char misc_run_match[1000];
  strcpy(misc_run_match, impossible_match); // impossible to match
  char modification_match[] = "_modification";

  char aa_modification[4000];
  char term_modification[4000];

  for(int k = 0; k < xmlfiles_->length(); k++) {
    // first = True;

    //cout << "ready to open inputfile: " << (*xmlfiles_)[k] << endl; exit(1);
  
    RACI fin((*xmlfiles_)[k]); // can read gzipped xml
    if(! fin) {
      cout << "ASAPRatioProteinCGIDisplayParser: error opening " << (*xmlfiles_)[k] << endl;
      exit(1);
    }

    Boolean header = True;
    Boolean filterParams = False;
    Boolean quantHighBG = False;
    Boolean zeroBG = False;
    double mzBound = 0.5;
    bool wavelet = false;
    char index_match[100];
    if(current_index_ < components_->length())
      sprintf(index_match, "index=\"%d\"", (*components_)[current_index_]->data_index_);
    else
      strcpy(index_match, impossible_match); // something impossible to match

    char asap_match[] = "asapratio_timestamp";
    //char xpress_match[] = "xpressratio_timestamp";

    int msms_run_idx = -1;
    //cout << (*xmlfiles_)[k] << endl;
    while(current_index_ < components_->length() && 
	  (*components_)[current_index_]->xml_index_ == k && 
	  fin.getline(nextline, line_width_)) {
     
      if(header || analyze || strstr(nextline, index_match) != NULL || strstr(nextline, asap_match) != NULL ||
	 strstr(nextline, misc_run_match) != NULL || strstr(nextline, modification_match) != NULL ||
	 strstr(nextline, "search_summary") != NULL || strstr(nextline, "analysis_timestamp") ||
	 strstr(nextline, "search_database") || strstr(nextline, "analysis_summary") ||
	 /*	 strstr(nextline, xpress_match) != NULL ||*/ strstr(nextline, "msms_run_summary") != NULL) {// ||
	//strstr(nextline, "search_result") != NULL) { 

	data = strstr(nextline, "<");
	while(data != NULL) {
	  tag = new Tag(data);

	  //tag->write(cout);
	  //cout << tag->getName() << endl;
	  if(tag != NULL) {

	    if(header && tag->isStart() && ! strcmp(tag->getName(), "asapratio_summary")) {
	      elutions_->insertAtEnd(atoi(tag->getAttributeValue("elution")));
	      //cout << "elutions: " << elutions_->length() << ", asaptime: " << asap_timestamps_->length() << endl;
	    }
	    else if(header && tag->isStart() && // ! strcmp(tag->getName(), "asapratio_timestamp") &&
		    ! strcmp(tag->getName(), "analysis_summary") &&
		    ! strcmp(tag->getAttributeValue("analysis"), "asapratio") &&
		    elutions_->length() == asap_timestamps_->length()) {
	      char* next = new char[strlen(tag->getAttributeValue("time"))+1];
	      strcpy(next, tag->getAttributeValue("time"));
	      asap_timestamps_->insertAtEnd(next);
	      quantHighBG = False;
	      zeroBG = False;
	      wavelet = false;
	      mzBound = 0.5;
	      filterParams = True;
	    }
	    else if(header && tag->isStart() && //! strcmp(tag->getName(), "peptideprophet_timestamp")) {
		    ! strcmp(tag->getName(), "analysis_summary") &&
		    ! strcmp(tag->getAttributeValue("analysis"), "peptideprophet")) {
	      char* next = new char[strlen(tag->getAttributeValue("time"))+1];
	      strcpy(next, tag->getAttributeValue("time"));
	      pepproph_timestamps_->insertAtEnd(next);
	    }
	    else if(header && tag->isStart() && //! strcmp(tag->getName(), "interprophet_timestamp")) {
		    ! strcmp(tag->getName(), "analysis_summary") &&
		    ! strcmp(tag->getAttributeValue("analysis"), "interprophet")) {
	      char* next = new char[strlen(tag->getAttributeValue("time"))+1];
	      strcpy(next, tag->getAttributeValue("time"));
	      iproph_timestamps_->insertAtEnd(next);
	    }
	    
	    if (filterParams && ! strcmp(tag->getName(), "parameter")) {
	      if (! strcmp(tag->getAttributeValue("name"), "quantHighBG") && 
		  ! strcmp(tag->getAttributeValue("value"), "True")) {
		quantHighBG = True;
	      }
	      else if (! strcmp(tag->getAttributeValue("name"), "zeroBG") && 
		       ! strcmp(tag->getAttributeValue("value"), "True")) {
		zeroBG = True;
	      }
	      else if (! strcmp(tag->getAttributeValue("name"), "wavelet") && 
		       ! strcmp(tag->getAttributeValue("value"), "True")) {
		wavelet = true;
	      }
	      else if (! strcmp(tag->getAttributeValue("name"), "mzBound")) {
		mzBound = atof(tag->getAttributeValue("value"));
		if (mzBound <= 0 || mzBound >= 1) {
		  mzBound = 0.5;
		}
	      }
	    }

	    if (filterParams && tag->isEnd() && ! strcmp(tag->getName(), "analysis_summary")) {
	      filterParams = False;
	      asap_quantHighBGs_->insertAtEnd(quantHighBG);
	      asap_zeroBGs_->insertAtEnd(zeroBG);
	      asap_wavelets_->insertAtEnd(wavelet);
	      asap_mzBounds_->insertAtEnd(mzBound);
	    }

	    if(tag->isStart() && ! strcmp(tag->getName(), "msms_run_summary")) {
	      header = False; // done
	      char* next = new char[strlen(tag->getAttributeValue("base_name"))+1];
	      strcpy(next, tag->getAttributeValue("base_name"));
	      basenames_->insertAtEnd(next);
	      msms_run_idx++;
	      //cout << "basenames: " << basenames_->length() << " with " << next << endl;

	      //first = False;
	    }
	    else if(tag->isStart() && ! strcmp(tag->getName(), "search_database")) {
	      char* next = new char[strlen(tag->getAttributeValue("local_path"))+1];
	      strcpy(next, tag->getAttributeValue("local_path"));	      
	      databases_->insertAtEnd(next);
	    }

	    else if(tag->isStart() && ! strcmp(tag->getName(), "search_summary")) {
	      strcpy(search_engine, tag->getAttributeValue("search_engine"));

	      if(! strcasecmp(search_engine, "COMET"))
		strcpy(misc_run_match, comet_md5_match);
	      else
		strcpy(misc_run_match, impossible_match); // impossible to match
	      int masstype = ! strcmp(tag->getAttributeValue("precursor_mass_type"), "monoisotopic") ? 1 : 0;
	      int fragmasstype = ! strcmp(tag->getAttributeValue("fragment_mass_type"), "monoisotopic") ? 1 : 0;
	      if(masstype || fragmasstype)
		(*components_)[current_index_]->setMassType(masstype, fragmasstype);
	      // reset
	      aa_modification[0] = 0;
	      term_modification[0] = 0;
	    }
	    else if(tag->isStart() && ! strcmp(tag->getName(), "parameter") && ! strcmp(tag->getAttributeValue("name"), "md5_check_sum")) { // comet case
	      char* nextmisc = new char[strlen(tag->getAttributeValue("value"))+1];
	      strcpy(nextmisc, tag->getAttributeValue("value"));
	      misc_run_conditions_->insertAtEnd(nextmisc);
	    }
#ifndef USE_STD_MODS
	    else if(tag->isStart() && ! strcmp(tag->getName(), "aminoacid_modification")) {
	      strcat(aa_modification, tag->getAttributeValue("aminoacid"));
	      if(tag->getAttributeValue("symbol") != NULL)
		strcat(aa_modification, tag->getAttributeValue("symbol"));
	      strcat(aa_modification, "-");
	      strcat(aa_modification, tag->getAttributeValue("mass"));
	      strcat(aa_modification, ":");
	    }
	    else if(tag->isStart() && ! strcmp(tag->getName(), "terminal_modification")) {
	      strcat(term_modification, tag->getAttributeValue("terminus"));
	      if(tag->getAttributeValue("symbol") != NULL)
		strcat(term_modification, tag->getAttributeValue("symbol"));
	      strcat(term_modification, "-");
	      strcat(term_modification, tag->getAttributeValue("mass"));
	      strcat(term_modification, ":");
	    }
	    else if(tag->isEnd() && ! strcmp(tag->getName(), "search_summary")) {
	      char* next_aamod = new char[strlen(aa_modification)+1];
	      strcpy(next_aamod, aa_modification);
	      aa_modifications_->insertAtEnd(next_aamod);
	      char* next_termmod = new char[strlen(term_modification)+1];
	      strcpy(next_termmod, term_modification);
	      term_modifications_->insertAtEnd(next_termmod);
	    }
#endif
	    //else if(! header && tag->isStart() && ! strcmp(tag->getName(), "xpressratio_timestamp")) {
	    //  comp_data.xpresslight = atoi(tag->getAttributeValue("display_ref"));
	    //}
	    else if(! header && tag->isStart() && // ! strcmp(tag->getName(), "asapratio_timestamp")) {
		    ! strcmp(tag->getName(), "analysis_timestamp") &&
		    ! strcmp(tag->getAttributeValue("analysis"), "asapratio")) {

	      asap_index = getTimestampIndex(tag->getAttributeValue("time"));
	      // cout << "asap time index for " << tag->getAttributeValue("time") << " is " << asap_index << endl;
	    }
	    else if(tag->isStart() && ! strcmp(tag->getName(), "spectrum_query")) {
	      precursor_charge = atoi(tag->getAttributeValue("assumed_charge"));
	      scan = (long)(atoi(tag->getAttributeValue("start_scan")));
	      comp_data.asapratio_index = atoi(tag->getAttributeValue("index"));

	      //cerr << "scan: " << scan << endl;
	      //tag->write(cerr);
	      index = atoi(tag->getAttributeValue("index"));
	      analyze = index == (*components_)[current_index_]->data_index_;
	      if(analyze) {
		search_result_tags = new Array<Tag*>;
		search_result_tags->insertAtEnd(tag);
		//cout << "entering search result for index " << index << endl;
	      }
	    }
	    else if(analyze) {
	      if(tag->isStart() && ! strcmp(tag->getName(), "asapratio_result")) {
		if(tags == NULL)
		  tags = new Array<Tag*>;
		tags->insertAtEnd(tag);
		collection = True;
		comp_data.asap_inv_mean = atof(tag->getAttributeValue("heavy2light_mean"));
		comp_data.asap_inv_err = atof(tag->getAttributeValue("heavy2light_error"));
		comp_data.asap_mean = atof(tag->getAttributeValue("mean"));
		comp_data.asap_err = atof(tag->getAttributeValue("error"));
		//comp_data.asapratio_index = atoi(tag->getAttributeValue("index"));
	      }
	      //else if(tag->isStart() && ! strcmp(tag->getName(), "search_hit") && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) {
	      //	tags = new Array<Tag*>;
	      // }
	      else if(tag->isStart() && ! strcmp(tag->getName(), "peptideprophet_result")) {
		probability = atof(tag->getAttributeValue("probability"));
		if(tag->getAttributeValue("analysis") != NULL) {
		  adjusted_prob = ! strcmp(tag->getAttributeValue("analysis"), "adjusted");
		  incomplete_prob = ! strcmp(tag->getAttributeValue("analysis"), "incomplete");
		}
		else {
		  adjusted_prob = False;
		  incomplete_prob = False;
		}
		delete tag;
	      }
	      else if(tag->isStart() && ! strcmp(tag->getName(), "interprophet_result")) {
		probability = atof(tag->getAttributeValue("probability"));
		if(tag->getAttributeValue("analysis") != NULL) {
		  adjusted_prob = ! strcmp(tag->getAttributeValue("analysis"), "adjusted");
		  incomplete_prob = ! strcmp(tag->getAttributeValue("analysis"), "incomplete");
		}
		else {
		  adjusted_prob = False;
		  incomplete_prob = False;
		}
		delete tag;
	      }
	      else if(tag->isStart() && ! strcmp(tag->getName(), "search_score_summary")) {
		search_score_found = True;
	      }
	      else if(search_score_found && tag->isStart() && ! strcmp(tag->getName(), "parameter")) {
		if(strlen(comp_data.score_summary) == 0) { // first one
		  strcpy(comp_data.score_summary, tag->getAttributeValue("name"));
		}
		else {
		  strcat(comp_data.score_summary, " ");
		  strcat(comp_data.score_summary, tag->getAttributeValue("name"));
		}
		strcat(comp_data.score_summary, ":");
		strcat(comp_data.score_summary, tag->getAttributeValue("value"));
	      }
	      else if(search_score_found && tag->isEnd() && ! strcmp(tag->getName(), "search_score_summary")) 
		search_score_found = False;
	      else if(tag->isEnd() && ! strcmp(tag->getName(), "asapratio_result")) {
		// process
		tags->insertAtEnd(tag);
		//cout << "setting pepdatastructs with elution: " << (*elutions_)[asap_index] << endl;
		//cout << "setting with elution index " << asap_index << " of " << elutions_->length() << endl; exit(1);

		setPepDataStruct(tags, (*elutions_)[asap_index], scan, precursor_charge);
		collection = False;
	      }

	      else if(tag->isStart() && ! strcmp(tag->getName(), "xpressratio_result")) {
		if(heavy2light_)
		  strcpy(comp_data.xpressratio, tag->getAttributeValue("heavy2light_ratio"));
		else
		  strcpy(comp_data.xpressratio, tag->getAttributeValue("ratio"));

		comp_data.lightfirstscan = atoi(tag->getAttributeValue("light_firstscan"));
		comp_data.lightlastscan = atoi(tag->getAttributeValue("light_lastscan"));
		comp_data.heavyfirstscan = atoi(tag->getAttributeValue("heavy_firstscan"));
		comp_data.heavylastscan = atoi(tag->getAttributeValue("heavy_lastscan"));
		comp_data.lightmass = atof(tag->getAttributeValue("light_mass"));
		comp_data.heavymass = atof(tag->getAttributeValue("heavy_mass"));
		comp_data.masstol = atof(tag->getAttributeValue("mass_tol"));
		//cout << comp_data.xpressratio << endl;

		delete tag;
	      }
	      else if(collection)
		tags->insertAtEnd(tag);
	      else {
		if(tag->isEnd() && ! strcmp(tag->getName(), "search_hit")) {

		  // here must process everything and reset
		  nextresult = getSearchResult(search_result_tags, search_engine);

		  if(nextresult != NULL && asap_index >= 0) {
		    //cout << "ready to process xml file " << k << " index " << index << endl;

		    nextresult->probability_ = probability;
		    nextresult->adjusted_prob_ = adjusted_prob;
		    nextresult->incomplete_prob_ = incomplete_prob;
		   
		    while(current_index_ < components_->length() &&
			  //DDS:
			  //(*components_)[current_index_]->msms_run_idx_  == msms_run_idx &&
			  (*components_)[current_index_]->xml_index_ == k &&
			  (*components_)[current_index_]->data_index_ == index) {
		      (*components_)[current_index_++]->enter(data_, nextresult, basenames_->length()-1, pepproph_timestamps_->length()-1, iproph_timestamps_->length()-1, databases_->length()-1, asap_index, comp_data, scan);
		    }

		    //  ; //   (*components_)[current_index_++]->enter(data_, nextresult, basenames_->length()-1, pepproph_timestamps_->length()-1, databases_->length()-1, asap_index, xpressratio, asap_mean, asap_err, asapratio_index, score_summary);

		    if(current_index_ < components_->length())
		      sprintf(index_match, "index=\"%d\"", (*components_)[current_index_]->data_index_);
		    //cout << "next match: " << index_match << endl;
		  }
		  else {
		    cout << "error" << endl;
		    exit(1);
		  }

		  // reset
		  analyze = False;
		  index = -1;
		  nextresult = NULL;
		  //comp_data.xpressratio = -1.0;
		  probability = -4.0;
		  comp_data.asap_mean = -3.0;
		  comp_data.asap_err = -1.0;
		  comp_data.asap_inv_mean = -3.0;
		  comp_data.asap_inv_err = -1.0;
		  comp_data.asapratio_index = -1; 
		  comp_data.score_summary[0] = 0;
		  comp_data.xpressratio[0] = 0;
		  // get rid of tags and search_result_tags
		  if(tags != NULL) {
		    for(int k = 0; k < tags->length(); k++)
		      if((*tags)[k] != NULL) {
			delete (*tags)[k];
			(*tags)[k] = NULL;
		      }
		    delete tags;
		    tags = NULL;
		  }

		  if(search_result_tags != NULL) {
		    for(int k = 0; k < search_result_tags->length(); k++)
		      if((*search_result_tags)[k] != NULL)
			delete (*search_result_tags)[k];
		    delete search_result_tags;
		    search_result_tags = NULL;
		  }
		}
		else
		  search_result_tags->insertAtEnd(tag);
	      } // if process (but not collection)
	    } // if analyze
	    else
	      delete tag;
	  } //  if not null
	  data = strstr(data+1, "<");
	} // next tag
      } // if possible peptide present

    } // next line
    //cout << "closing " << (*xmlfiles_)[k] << endl;
    fin.close();
  } // next inputfile

  //cout << "done with group" << endl;
  //pro_ = ratio_->getProDataStruct();
  //fout.close();
  delete[] nextline;
}


pepDataStrct ASAPRatioProteinCGIDisplayParser::getPepDataStrct(int seq, int pk, int data) {
  if(pepdata_index_ < components_->length() &&
     (*components_)[pepdata_index_]->seq_ == seq &&
     (*components_)[pepdata_index_]->peak_ == pk &&
     (*components_)[pepdata_index_]->data_ == data) 
    return (*components_)[pepdata_index_++]->pepdata_;
  int k;
  for(k = pepdata_index_; k < components_->length(); k++)
    if((*components_)[k]->seq_ == seq &&
       (*components_)[k]->peak_ == pk &&
       (*components_)[k]->data_ == data) {
      pepdata_index_ = k;
      return (*components_)[pepdata_index_++]->pepdata_;
    }

  for(k = 0; k < components_->length() && pepdata_index_; k++)
    if((*components_)[k]->seq_ == seq &&
       (*components_)[k]->peak_ == pk &&
       (*components_)[k]->data_ == data) {
      pepdata_index_ = k;
      return (*components_)[pepdata_index_++]->pepdata_;
    }
  // error
  cout << "error: not component found for seq: " << seq << " pk: " << pk << " data: " << data << endl;
  exit(1);
}

int ASAPRatioProteinCGIDisplayParser::getResultIndex(int seq, int pk, int data, Boolean advance) {
  if(result_index_ < components_->length() &&
     (*components_)[result_index_]->seq_ == seq &&
     (*components_)[result_index_]->peak_ == pk &&
     (*components_)[result_index_]->data_ == data){ 
    if(advance){
      return result_index_++;
    }
    else{ 
      return result_index_;
    }
  }
  int k;
  for(k = result_index_; k < components_->length(); k++)
    if((*components_)[k]->seq_ == seq &&
       (*components_)[k]->peak_ == pk &&
       (*components_)[k]->data_ == data) {
      result_index_ = k;
      if(advance)
	return result_index_++;
      else 
	return result_index_;
    }

  for(k = 0; k < components_->length() && result_index_; k++)
    if((*components_)[k]->seq_ == seq &&
       (*components_)[k]->peak_ == pk &&
       (*components_)[k]->data_ == data) {
      result_index_ = k;
      if(advance)
	return result_index_++;
      else 
	return result_index_;
    }

  return -1;
}

char* ASAPRatioProteinCGIDisplayParser::display(int seq, int pk, int data) {
  if(result_index_ < components_->length() &&
     (*components_)[result_index_]->seq_ == seq &&
     (*components_)[result_index_]->peak_ == pk &&
     (*components_)[result_index_]->data_ == data) 
    return (*components_)[result_index_++]->display(basenames_, pepproph_timestamps_, iproph_timestamps_, databases_, asap_timestamps_);
  int k;
  for(k = result_index_; k < components_->length(); k++)
    if((*components_)[k]->seq_ == seq &&
       (*components_)[k]->peak_ == pk &&
       (*components_)[k]->data_ == data) {
      result_index_ = k;
      return (*components_)[result_index_++]->display(basenames_, pepproph_timestamps_, iproph_timestamps_, databases_, asap_timestamps_);
    }

  for(k = 0; k < result_index_; k++)
    if((*components_)[k]->seq_ == seq &&
       (*components_)[k]->peak_ == pk &&
       (*components_)[k]->data_ == data) {
      result_index_ = k;
      return (*components_)[result_index_++]->display(basenames_, pepproph_timestamps_, iproph_timestamps_, databases_, asap_timestamps_);
    }
  return NULL;
}


SearchResult* ASAPRatioProteinCGIDisplayParser::getSearchResult(Array<Tag*>* tags, char* engine) {
  MixtureDistrFactory fac;

  return fac.getSearchResult(tags, engine);
}

void ASAPRatioProteinCGIDisplayParser::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query")) {
    if(tag->isStart() && ! strcmp(tag->getAttributeValue("hit_rank"), "1")) {
      //tag->print();
      filter_ = True;
    }
    else {
      if(filter_ && tag->isEnd())
        filter_memory_ = True;
    }
  }
}


proDataStrct* ASAPRatioProteinCGIDisplayParser::getProDataStruct() {
  return protein_;
}

void ASAPRatioProteinCGIDisplayParser::setPepDataStruct(Array<Tag*>* tags, int elution, long scan, int precursor_charge) {
  Tag* next;
  int charge=0;

  //DDS: calculate pepArea variable
  bool use_this = false;
  double maxArea = 0;
  double tmpArea = 0;

  double tmp_ltime = 0;
  double tmp_ltime_wd = 0;
  double tmp_htime = 0;
  double tmp_htime_wd = 0;

  double ltime = 0;
  double ltime_wd = 0;
  double htime = 0;
  double htime_wd = 0;

  data_.scan = scan;
  data_.chrg = precursor_charge;
  data_.eltn = elution;

  for(int k = 0; k < tags->length(); k++) {
    next = (*tags)[k];
    if(next->isStart()) {
      if(! strcmp(next->getName(), "asapratio_result")) {
	data_.pepRatio[0] = atof(next->getAttributeValue("mean"));
	data_.pepRatio[1] = atof(next->getAttributeValue("error"));
	data_.pepH2LRatio[0] = atof(next->getAttributeValue("heavy2light_mean"));
	data_.pepH2LRatio[1] = atof(next->getAttributeValue("heavy2light_error"));
	// make change here
	if(data_.pepRatio[0] == -1) {
	  data_.pepRatio[0] = -2;
	  data_.pepH2LRatio[0] = -2;
	}
	else if(data_.pepRatio[0] >= 999.0) {
	  data_.pepRatio[0] = -1;
	  data_.pepH2LRatio[0] = -1;
	}
	
      }
      else if(! strcmp(next->getName(), "asapratio_peptide_data")) {
	data_.indx = atoi(next->getAttributeValue("status"));
	data_.cidIndx = atoi(next->getAttributeValue("cidIndex"));
	data_.msLight = atof(next->getAttributeValue("light_mass"));
	data_.msHeavy = atof(next->getAttributeValue("heavy_mass"));
	data_.areaFlag = atoi(next->getAttributeValue("area_flag"));
      }
      else if(! strcmp(next->getName(), "asapratio_contribution")) {
	charge = atoi(next->getAttributeValue("charge"));
	data_.pkRatio[charge-1] = atof(next->getAttributeValue("ratio"));
	data_.pkError[charge-1] = atof(next->getAttributeValue("error"));
	data_.pkCount[charge-1] = atoi(next->getAttributeValue("use"));
	
	if (use_this && tmpArea > maxArea) {
	  maxArea = tmpArea; 
	  ltime = tmp_ltime;
	  htime = tmp_htime;
	  ltime_wd = tmp_ltime_wd;
	  htime_wd = tmp_htime_wd;
	}
	
	if (data_.pkCount[charge-1] == 1) {
	  use_this = true;
	  tmpArea = 0; 
	  tmp_ltime = 0;
	  tmp_htime = 0;
	  tmp_ltime_wd = 0;
	  tmp_htime_wd = 0;
	}
	else {
	  use_this = false;
	  tmpArea = 0;
	  tmp_ltime = 0;
	  tmp_htime = 0;
	  tmp_ltime_wd = 0;
	  tmp_htime_wd = 0;
	}

      }
      else if(! strcmp(next->getName(), "asapratio_lc_lightpeak")) {
	int label = 0;
	data_.peaks[charge-1][label].indx = atoi(next->getAttributeValue("status"));
	data_.peaks[charge-1][label].valley[0] = atoi(next->getAttributeValue("left_valley"));
	data_.peaks[charge-1][label].valley[1] = atoi(next->getAttributeValue("right_valley"));
	data_.peaks[charge-1][label].bckgrnd = atof(next->getAttributeValue("background"));
	data_.peaks[charge-1][label].area[0] = atof(next->getAttributeValue("area"));
	data_.peaks[charge-1][label].area[1] = atof(next->getAttributeValue("area_error"));
	data_.peaks[charge-1][label].time[0] = atof(next->getAttributeValue("time"));
	data_.peaks[charge-1][label].time[1] = atof(next->getAttributeValue("time_width"));
	data_.peaks[charge-1][label].peak = atoi(next->getAttributeValue("is_heavy"));
	if (use_this) {
	  tmpArea += data_.peaks[charge-1][label].area[0];
	  tmp_ltime = data_.peaks[charge-1][label].time[0];
	  tmp_ltime_wd = data_.peaks[charge-1][label].time[1];
	}
      }
      else if(! strcmp(next->getName(), "asapratio_lc_heavypeak")) {
	int label = 1;
	data_.peaks[charge-1][label].indx = atoi(next->getAttributeValue("status"));
	data_.peaks[charge-1][label].valley[0] = atoi(next->getAttributeValue("left_valley"));
	data_.peaks[charge-1][label].valley[1] = atoi(next->getAttributeValue("right_valley"));
	data_.peaks[charge-1][label].bckgrnd = atof(next->getAttributeValue("background"));
	data_.peaks[charge-1][label].area[0] = atof(next->getAttributeValue("area"));
	data_.peaks[charge-1][label].area[1] = atof(next->getAttributeValue("area_error"));
	data_.peaks[charge-1][label].time[0] = atof(next->getAttributeValue("time"));
	data_.peaks[charge-1][label].time[1] = atof(next->getAttributeValue("time_width"));
	data_.peaks[charge-1][label].peak = atoi(next->getAttributeValue("is_heavy"));
	if (use_this) {
	  tmpArea += data_.peaks[charge-1][label].area[0];
	  tmp_htime = data_.peaks[charge-1][label].time[0];
	  tmp_htime_wd = data_.peaks[charge-1][label].time[1];
	}
      }
    } // if start
  } // next tag
  if (use_this && tmpArea > maxArea) {
    maxArea = tmpArea;
    ltime = tmp_ltime;
    htime = tmp_htime;
    ltime_wd = tmp_ltime_wd;
    htime_wd = tmp_htime_wd;
  }
  data_.pepArea = maxArea;
  data_.pepTime[0][0] = ltime;
  data_.pepTime[0][1] = ltime_wd;
  data_.pepTime[1][0] = htime;
  data_.pepTime[1][1] = htime_wd;
}
