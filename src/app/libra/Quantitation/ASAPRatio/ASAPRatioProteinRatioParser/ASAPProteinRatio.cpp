/*
Program       : ASAPRatioPeptideParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPProteinRatio.cpp 7996 2019-12-25 00:16:42Z real_procopio $

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

#include "ASAPProteinRatio.h"

ASAPProteinRatio::ASAPProteinRatio(double minprob, double minwt) {
  runseqdata_ = new Array<ASAPRatioMSMSSeqData*>;
  set_ = False;
  start_ = 0;
  pro_ = NULL;
  lightlabels_ = NULL;
  heavylabels_ = NULL;
  min_probability_ = minprob;
  min_weight_ = minwt;
}

ASAPProteinRatio::~ASAPProteinRatio() {
  if(pro_ != NULL) {
    for (int p = 0; p < pro_->dataNum; ++p) {
      for (int s = 0; s < pro_->sequences[p].dataNum; ++s) {
	free(pro_->sequences[p].peaks[s].dataIndx);
	free(pro_->sequences[p].peaks[s].dataCnts);
      }
      free(pro_->sequences[p].peaks);
      free(pro_->sequences[p].dataCnts);
    }
    free(pro_->sequences);
    free(pro_->dataCnts);

    delete pro_;
  }

  if(runseqdata_ != NULL) {
    for(int k = 0; k < runseqdata_->length(); k++)
      if((*runseqdata_)[k] != NULL)
	delete (*runseqdata_)[k];
    delete runseqdata_;
  }

  if(lightlabels_ != NULL)
    delete lightlabels_;
  if(heavylabels_ != NULL)
    delete heavylabels_;
}

void ASAPProteinRatio::setRunInfo(char* lightlabels, char* heavylabels) {
  if(lightlabels_ == NULL || heavylabels_ == NULL || 
     strcmp(lightlabels_, lightlabels) || strcmp(heavylabels_, heavylabels)) {
    if(lightlabels_ != NULL)
      delete lightlabels_;
    lightlabels_ = new char[strlen(lightlabels)+1];
    strcpy(lightlabels_, lightlabels);
    if(heavylabels_ != NULL)
      delete heavylabels_;
    heavylabels_ = new char[strlen(heavylabels)+1];
    strcpy(heavylabels_, heavylabels);

    label_partners_ = collectPrtnAAStrct(&prtnAANum_, lightlabels_, heavylabels_);
  }

  start_ = runseqdata_->length();
}


#ifdef USE_STD_MODS
void ASAPProteinRatio::enter(const char* peptide, double lightmass, const ModificationInfo* modinfo, const char* quant_labels, 
			     const pepDataStrct &data, int index, int xml_index, double wt, double prob, int msms_run_idx) {
  for(int k = start_; k < runseqdata_->length(); k++) {
    Boolean found = False;
    // first check for match of stripped pep, quantlabels, and lightmass
    if(! strcmp((*runseqdata_)[k]->lightsequence_, peptide) && 
       ! strcmp((*runseqdata_)[k]->quant_labels_, quant_labels) &&
       (*runseqdata_)[k]->light_mass_ - lightmass <= MOD_ERROR &&
       lightmass - (*runseqdata_)[k]->light_mass_ <= MOD_ERROR) { // now must check modifications more closely

      if((*runseqdata_)[k]->light_mass_ == 0.0) { // not valid quant, must match peptide and mods exactly
	if((*runseqdata_)[k]->mod_info_->equivalentModification(modinfo, MOD_ERROR, peptide, NULL))
	  found = True;
      } // non-quant pep
      else { // check without regard to quant label positions
	if((*runseqdata_)[k]->mod_info_->equivalentModification(modinfo, MOD_ERROR, peptide, quant_labels))
	  found = True;
      }
      if(found && 
	 //DDS:
	 (*runseqdata_)[k]->xml_index_ == xml_index &&
	 (*runseqdata_)[k]->msms_run_idx_ == msms_run_idx
	 ) {
	(*runseqdata_)[k]->data_->insertAtEnd(data);
	(*runseqdata_)[k]->indices_->insertAtEnd(index);
	(*runseqdata_)[k]->probs_->insertAtEnd(prob);
	//cout << "adding index " << index << " to peptide " << light <<  " (" << (*runseqdata_)[k]->lightsequence_ << ")" << endl;
	return;
      }
    } // if potential match
  } // next entry
  // still here
  runseqdata_->insertAtEnd(new ASAPRatioMSMSSeqData(peptide, lightmass, modinfo, quant_labels, data, index, xml_index, wt, prob, msms_run_idx));
}
#else
void ASAPProteinRatio::enter(const char* peptide, const pepDataStrct & data, int index, int xml_index, double wt, double prob, int msms_run_idx) {
  char* light = getLightSequence(peptide, label_partners_, prtnAANum_);
  for(int k = start_; k < runseqdata_->length(); k++)
    if(! strcmp(light, (*runseqdata_)[k]->lightsequence_) 
       //DDS:
       (*runseqdata_)[k]->xml_index_ == xml_index &&
       (*runseqdata_)[k]->msms_run_idx_ == msms_run_idx) {
      (*runseqdata_)[k]->data_->insertAtEnd(data);
      (*runseqdata_)[k]->indices_->insertAtEnd(index);
      (*runseqdata_)[k]->probs_->insertAtEnd(prob);
      //cout << "adding index " << index << " to peptide " << light <<  " (" << (*runseqdata_)[k]->lightsequence_ << ")" << endl;
      delete light;
      return;
    }

  // still here
  //cout << "making new entry with index " << index << " for peptide " << light << endl;
  runseqdata_->insertAtEnd(new ASAPRatioMSMSSeqData(light, data, index, /* basename_,*/ xml_index, wt, prob, msms_run_idx));
  delete light;
}
#endif


// the guts
void ASAPProteinRatio::computeProteinRatio() {
  sortMSMSSeqData();
  initProDataStrct();
  getProDataStrct(pro_);

  /*
  // display
  cout << "MSMSSEQDATA with MINPEPPROB " << min_probability_ << " and MINPEPWT " << min_weight_ << endl;
  for(int k = 0; k < runseqdata_->length(); k++) {
    cout << (*runseqdata_)[k]->lightsequence_ << " (" << (*runseqdata_)[k]->xml_index_ << ", " << (*runseqdata_)[k]->weight_ << "): ";
    for(int j = 0; j < (*runseqdata_)[k]->indices_->length(); j++)
      cout << (*(*runseqdata_)[k]->indices_)[j] << "[" << (*(*runseqdata_)[k]->probs_)[j] << "] ";
    cout << endl;
  }
  cout << "----------------------------------------------" << endl;
  */

  set_ = True;
}


// This function initializes a proDataStrct *protein.
void ASAPProteinRatio::initProDataStrct() {
  int seqNum;
  int startIndx, endIndx;
  int i;

  // initialize protein
  pro_ = new proDataStrct();
  pro_->indx = -1;
  pro_->sequences = (seqDataStrct *) calloc(runseqdata_->length(), sizeof(seqDataStrct));
  pro_->dataCnts = (int *) calloc(runseqdata_->length(), sizeof(int));

  seqNum = 0;
  endIndx = 0;
  while(endIndx < runseqdata_->length()){
    startIndx = endIndx;
    while(endIndx < runseqdata_->length() 
	  && strcmp((*runseqdata_)[endIndx]->lightsequence_, (*runseqdata_)[startIndx]->lightsequence_) == 0)
      ++endIndx;

    initSeqDataStrct(pro_->sequences+seqNum, startIndx, endIndx, min_probability_);
    
    if (pro_->sequences[seqNum].indx != -1)
      pro_->indx = 0;
    
    pro_->dataCnts[seqNum] = 0;
    for (i = 0; i < pro_->sequences[seqNum].dataNum; ++i){
      if (pro_->sequences[seqNum].dataCnts[i] == 1)
	pro_->dataCnts[seqNum] = 1;
    }
    if(pro_->sequences[seqNum].weight <= min_weight_)
      pro_->dataCnts[seqNum] = 0;

    ++seqNum;
  }
  pro_->dataNum = seqNum;
  pro_->sequences = (seqDataStrct *) realloc(pro_->sequences, seqNum*sizeof(seqDataStrct));
  pro_->dataCnts = (int *) realloc(pro_->dataCnts, seqNum*sizeof(int));

}


// This function get the initial seqDataStrct of a unique peptide sequence.
void ASAPProteinRatio::initSeqDataStrct(seqDataStrct *seqData, int start_index, int end_index, double mnProb)
{
  int totPepNum, totPepCount;
  int **pepIndx;
  int peakNum, peakIndx;
  double timeRange[2][2];
  double *weight;
  int i, j;
  
  // total peptide number
  totPepNum = 0;
  for (i = start_index; i < end_index; ++i){
    totPepNum += (*runseqdata_)[i]->data_->length();
  }  

  int dreamNum = end_index - start_index;

  // get seqData->peaks

  // memory
  seqData->peaks = (dataStrct *) calloc(totPepNum, sizeof(dataStrct));
  weight = (double *) calloc(totPepNum, sizeof(double));

  // pepIndx
  pepIndx = (int **) calloc(dreamNum, sizeof(int *));
  for (i = 0; i < dreamNum; ++i)
    pepIndx[i] = (int *) calloc((*runseqdata_)[start_index + i]->data_->length(), sizeof(int));
  for (i = 0; i < dreamNum; ++i){
    for (j = 0; j < (*runseqdata_)[start_index + i]->data_->length(); ++j){
      pepIndx[i][j] = 0;
    }
  }

  // collect data
  peakNum = 0;
  totPepCount = 0;
  while(totPepCount < totPepNum){
    seqData->peaks[peakNum].dataNum = 0;
    seqData->peaks[peakNum].dataIndx = (int *) calloc(totPepNum, sizeof(int));
    seqData->peaks[peakNum].dataCnts = (int *) calloc(totPepNum, sizeof(int));
    seqData->peaks[peakNum].weight = 0.;
    seqData->peaks[peakNum].bofIndx = -1;
    
    //DDS:
    seqData->peaks[peakNum].msms_run_idx = -1;
    
    weight[peakNum] = 0.;
    peakIndx = 1;
    for (i = 0; peakIndx == 1 && i < dreamNum; ++i){
      if(i == 0)
	strcpy(seqData->lightSeq, (*runseqdata_)[start_index + i]->lightsequence_);
      for (j = 0; j < (*runseqdata_)[start_index + i]->data_->length(); ++j){
	if(pepIndx[i][j] == 1)
	  continue;
	if(seqData->peaks[peakNum].bofIndx == -1){ // new peak
	  if(    (*((*runseqdata_)[start_index + i]->data_))[j].indx == -1){ // invalid data
	    seqData->peaks[peakNum].indx = -1;
	    seqData->peaks[peakNum].ratio[0] = -2.;
	    seqData->peaks[peakNum].ratio[1] = 0.;
	    seqData->peaks[peakNum].dataIndx[0] = (*((*runseqdata_)[start_index + i]->indices_))[j];
	    seqData->peaks[peakNum].dataCnts[0] = 0;
	    seqData->peaks[peakNum].dataNum = 1;
	    seqData->peaks[peakNum].weight = 0.;
	    seqData->peaks[peakNum].bofIndx = (*runseqdata_)[start_index + i]->xml_index_; 
	    
	    //DDS:
	    seqData->peaks[peakNum].msms_run_idx = (*runseqdata_)[start_index + i]->msms_run_idx_;
	    
	    weight[peakNum] = (*runseqdata_)[start_index + i]->weight_;
	    pepIndx[i][j] = 1;
	    peakIndx = 0;
	    ++totPepCount;
	    break;
	  }
	  else {
	    seqData->peaks[peakNum].bofIndx = (*runseqdata_)[start_index + i]->xml_index_;
	    seqData->peaks[peakNum].msms_run_idx = (*runseqdata_)[start_index + i]->msms_run_idx_;
	    timeRange[0][0] = (*((*runseqdata_)[start_index + i]->data_))[j].pepTime[0][0];
	    timeRange[0][1] = (*((*runseqdata_)[start_index + i]->data_))[j].pepTime[0][1];
	    timeRange[1][0] = (*((*runseqdata_)[start_index + i]->data_))[j].pepTime[1][0];
	    timeRange[1][1] = (*((*runseqdata_)[start_index + i]->data_))[j].pepTime[1][1];
	  }
	}

	if(pepIndx[i][j] == 0 // not yet analyzed
	   && (*runseqdata_)[start_index + i]->xml_index_ == seqData->peaks[peakNum].bofIndx // from the same run  
	   && (*runseqdata_)[start_index + i]->msms_run_idx_ == seqData->peaks[peakNum].msms_run_idx // from the same run 
	   && // eluted nearby: light
	   fabs((*((*runseqdata_)[start_index + i]->data_))[j].pepTime[0][0]-timeRange[0][0])
	   <= fabs((*((*runseqdata_)[start_index + i]->data_))[j].pepTime[0][1]+timeRange[0][1])
	   && // eluted nearby: heavy
	   fabs((*((*runseqdata_)[start_index + i]->data_))[j].pepTime[1][0]-timeRange[1][0])
	   <= fabs((*((*runseqdata_)[start_index + i]->data_))[j].pepTime[1][1]+timeRange[1][1])) {
	  seqData->peaks[peakNum].indx = 0;
	  seqData->peaks[peakNum].dataIndx[seqData->peaks[peakNum].dataNum] = 
	    (*((*runseqdata_)[start_index + i]->indices_))[j];
	  if((*((*runseqdata_)[start_index + i]->probs_))[j] >= mnProb)
	    seqData->peaks[peakNum].dataCnts[seqData->peaks[peakNum].dataNum] = 1;
	  else
	    seqData->peaks[peakNum].dataCnts[seqData->peaks[peakNum].dataNum] = 0;
	  ++(seqData->peaks[peakNum].dataNum);
	  seqData->peaks[peakNum].weight 
	    = seqData->peaks[peakNum].weight > 
	    (*((*runseqdata_)[start_index + i]->data_))[j].pepArea ? 
	    seqData->peaks[peakNum].weight : (*((*runseqdata_)[start_index + i]->data_))[j].pepArea;	  
	  weight[peakNum] = weight[peakNum] > (*runseqdata_)[start_index + i]->weight_ ?
	    weight[peakNum] : (*runseqdata_)[start_index + i]->weight_;
	  pepIndx[i][j] = 1;
	  ++totPepCount;
	}
      }
    }
    
    seqData->peaks[peakNum].dataIndx = (int *) realloc(seqData->peaks[peakNum].dataIndx, 
						       seqData->peaks[peakNum].dataNum*sizeof(int));
    seqData->peaks[peakNum].dataCnts = (int *) realloc(seqData->peaks[peakNum].dataCnts, 
						       seqData->peaks[peakNum].dataNum*sizeof(int));
    ++peakNum;
  }
  seqData->peaks = (dataStrct *) realloc(seqData->peaks, peakNum*sizeof(dataStrct));
  for (i = 0; i < dreamNum; ++i)
    free(pepIndx[i]);
  free(pepIndx);

  // get seqdata
  seqData->indx = -1;
  seqData->dataNum = peakNum;
  seqData->dataCnts = (int *) calloc(peakNum, sizeof(int));
  seqData->weight = 0.;
  for (i = 0; i < peakNum; ++i) {
    if(seqData->peaks[i].indx != -1)
      seqData->indx = 0;
    seqData->dataCnts[i] = 0;
    for (j = 0; j < seqData->peaks[i].dataNum; ++j)
      if(seqData->peaks[i].dataCnts[j] == 1)
	seqData->dataCnts[i] = 1;
    seqData->weight = seqData->weight > weight[i] ? seqData->weight : weight[i];
  } //   for (i = 0; i < peakNum; ++i) {
  free(weight);

}

// This function evaluates the ratio of a protein.
void ASAPProteinRatio::getProDataStrct(proDataStrct *data)
{
  double *ratios, *errors, *inv_ratios, *inv_errors, *weights;
  int *outliers, *seqIndx;
  int seqNum, tmpNum;
  double cmnErr, ratio, error, inv_ratio, inv_error;
  int num[2] = {0, 0};

  int i, j;

  if (data->indx == -1) {
    data->ratio[0] = -2.;
    data->ratio[1] = 0.;
    data->inv_ratio[0] = -2.;
    data->inv_ratio[1] = 0.;
    return;
  }

  // collect all sequences
  seqIndx = (int *) calloc(data->dataNum, sizeof(int));
  seqNum = 0;
  for (i = 0; i < data->dataNum; ++i) {
    getSeqDataStrctRatio(&(data->sequences[i]));
    // double check dataCnts
    if(data->dataCnts[i] == 1) {
      data->dataCnts[i] = 0;
      for (j = 0; j < data->sequences[i].dataNum; ++j){
	if(data->sequences[i].dataCnts[j] == 1){
	  data->dataCnts[i] = 1;
	  break;
	}
      }
    }
    if(data->dataCnts[i] == 1 && data->sequences[i].indx != -1) {
      seqIndx[i] = 1;
      ++seqNum;
    }
    else
      seqIndx[i] = 0;
  }

  // only 0 or 1 valid sequences
  if(seqNum < 2) {
    if(seqNum < 1) {
      data->ratio[0] = -2.;
      data->ratio[1] = 0.;
      data->inv_ratio[0] = -2.;
      data->inv_ratio[1] = 0.;
    } // if(seqNum < 1) {
    else {
      for (i = 0; i < data->dataNum; ++i) {
	if(seqIndx[i] == 1) {
	  data->ratio[0] = data->sequences[i].ratio[0];
	  data->ratio[1] = data->sequences[i].ratio[1];
	  data->inv_ratio[0] = data->sequences[i].inv_ratio[0];
	  data->inv_ratio[1] = data->sequences[i].inv_ratio[1];
	  break;
	}
      }
    }

    if(data->indx == 0)
      data->indx = 1;
    
    free(seqIndx);

    return;
  }


  // allocate memory
  ratios = (double *) calloc(seqNum, sizeof(double));
  errors = (double *) calloc(seqNum, sizeof(double));
  inv_ratios = (double *) calloc(seqNum, sizeof(double));
  inv_errors = (double *) calloc(seqNum, sizeof(double));
  weights = (double *) calloc(seqNum, sizeof(double));
  outliers = (int *) calloc(seqNum, sizeof(int));

  // get data
  tmpNum = 0;
  for(i = 0; i < data->dataNum; ++i) {
    if(seqIndx[i] == 1) {
      ratios[tmpNum] = data->sequences[i].ratio[0];
      errors[tmpNum] = data->sequences[i].ratio[1];
      inv_ratios[tmpNum] = data->sequences[i].inv_ratio[0];
      inv_errors[tmpNum] = data->sequences[i].inv_ratio[1];
      ++tmpNum;
    }
  }

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
  }
  
  if(tmpNum < 1) { // all ratios are 0 or -1 or -2
    if(num[0] > num[1]) {// more 0 than -1
      ratio = 0.;
      inv_ratio = 999.;
    }
    else if(num[0] < num[1]) { // more -1 than 0
      ratio = -1.;
      inv_ratio = -1.;
    }
    else {// same 0 and -1
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
    data->ratio[0] = ratio;
    data->ratio[1] = error;

    data->inv_ratio[0] = inv_ratio;
    data->inv_ratio[1] = inv_error;

    if(data->indx == 0) {
      tmpNum = 0;
      for(i = 0; i < data->dataNum; ++i) {
	if(seqIndx[i] == 1) {
	  data->dataCnts[i] = 1 - outliers[tmpNum];
	  ++tmpNum;
	}
      }
      data->indx = 1;
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


  // get weight
  for(i = 0; i < seqNum; ++i) {
    if(errors[i] < cmnErr)
      weights[i] = 1./cmnErr/cmnErr;
    else
      weights[i] = 1./errors[i]/errors[i];
  }


  //calculate ratio and error
  if(data->indx == 0) {
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]),_ASAPRATIO_CONFL_,
		 ratios, errors, inv_ratios, inv_errors, weights, outliers, seqNum, 1);
    tmpNum = 0;
    for(i = 0; i < data->dataNum; ++i) {
      if(seqIndx[i] == 1) {
	data->dataCnts[i] = 1 - outliers[tmpNum];
	++tmpNum;
      }
    }      
  }
  else
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]),_ASAPRATIO_CONFL_,
		 ratios, errors, inv_ratios, inv_errors, weights, outliers, seqNum, 0);

  // reset indx
  if(data->indx == 0)
    data->indx = 1;

  // free memory
  free(seqIndx);
  free(ratios);
  free(errors);
  free(inv_ratios);
  free(inv_errors);
  free(weights);
  free(outliers);
}

// This function calculates the ratio of a unique sequence.
void ASAPProteinRatio::getSeqDataStrctRatio(seqDataStrct *data)
{
  double *ratios, *errors, *inv_ratios, *inv_errors, *weights;
  int *outliers, *peakIndx;
  int peakNum, tmpNum;

  int i, j;

  if (data->indx == -1) {
    data->ratio[0] = -2.;
    data->ratio[1] = 0.;
    data->inv_ratio[0] = -2.;
    data->inv_ratio[1] = 0.;
    return;
  }

  // collect all peaks
  peakIndx = (int *) calloc(data->dataNum, sizeof(int));
  peakNum = 0;
  for (i = 0; i < data->dataNum; ++i) {
    getDataStrctRatio(&(data->peaks[i])); 

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
    }
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
    }

    if(data->indx == 0)
      data->indx = 1;

    free(peakIndx);
    return;
  }


  // allocate memory
  ratios = (double *) calloc(peakNum, sizeof(double));
  errors = (double *) calloc(peakNum, sizeof(double));
  inv_ratios = (double *) calloc(peakNum, sizeof(double));
  inv_errors = (double *) calloc(peakNum, sizeof(double));
  weights = (double *) calloc(peakNum, sizeof(double));
  outliers = (int *) calloc(peakNum, sizeof(int));

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
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]), _ASAPRATIO_CONFL_,
		 ratios, errors, inv_ratios, inv_errors, weights, outliers, peakNum, 1);
    tmpNum = 0;
    for(i = 0; i < data->dataNum; ++i) {
      if(peakIndx[i] == 1) {
	data->dataCnts[i] = 1 - outliers[tmpNum];
	++tmpNum;
      }
    }      
  }
  else
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]), _ASAPRATIO_CONFL_,
		 ratios, errors, inv_ratios, inv_errors, weights, outliers, peakNum, 0);
  
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
}


// This function calculates the ratio of a unique peak.
void ASAPProteinRatio::getDataStrctRatio(dataStrct *data) 
{
  pepDataStrct *peptides;
  double *ratios, *errors, *inv_ratios, *inv_errors, *weights;
  int *outliers, *pepIndx;
  int pepNum, tmpNum;
  double tol = 1.e-2;

  int i;

  if (data->indx == -1) {
    data->ratio[0] = -2.;
    data->ratio[1] = 0.;
    data->inv_ratio[0] = -2.;
    data->inv_ratio[1] = 0.;
    return;
  }

  // collect all valid peptides

  // read peptides
  peptides = (pepDataStrct *) calloc(data->dataNum, sizeof(pepDataStrct));
  pepIndx = (int *) calloc(data->dataNum, sizeof(int));
  pepNum = 0;
  for (i = 0; i < data->dataNum; ++i) {
    if(data->dataCnts[i] == 1) { 
      searchPepDataStrct(peptides+pepNum, data->dataIndx[i], data->bofIndx);
      if(peptides[pepNum].indx != -1) {
	pepIndx[i] = 1;
	//cout << "got " << pepNum << ": " << peptides[pepNum].pepRatio[0] << endl;
	++pepNum;
      }
      else 
	pepIndx[i] = 0;
    }
    else
      pepIndx[i] = 0;
  }

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
    }

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
    //cout << i << ": here with " << peptides[i].pepRatio[0] << endl;
    errors[i] = peptides[i].pepRatio[1];
    inv_ratios[i] = peptides[i].pepH2LRatio[0];
    //cout << i << ": here with " << peptides[i].pepRatio[0] << endl;
    inv_errors[i] = peptides[i].pepH2LRatio[1];
    weights[i] = peptides[i].pepArea;
  }
    
  // calculate ratio and error
  if(data->indx == 0) {
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]), _ASAPRATIO_CONFL_,
		 ratios, errors, inv_ratios, inv_errors, weights, outliers, pepNum, 1);
    tmpNum = 0;
    for (i = 0; i < data->dataNum; ++i) {
      if(pepIndx[i] == 1) {
	data->dataCnts[i] = 1 - outliers[tmpNum];
	++tmpNum;
      }
    }      
  }
  else
    getDataRatio(&(data->ratio[0]), &(data->ratio[1]), &(data->inv_ratio[0]), &(data->inv_ratio[1]), _ASAPRATIO_CONFL_,
		 ratios, errors, inv_ratios, inv_errors, weights, outliers, pepNum, 0);

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
  }
  
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


// This function searches pepDataStrct from dreamStrct.
void ASAPProteinRatio::searchPepDataStrct(pepDataStrct *pepData, int pepIndx, int xmlIndx)
{
  int indx = 0;
  int i, j;

  for (i = 0; indx == 0 && i < runseqdata_->length(); ++i){
    if((*runseqdata_)[i]->xml_index_ != xmlIndx)
      continue;
    for (j = 0; indx == 0 && j < (*runseqdata_)[i]->data_->length(); ++j){
      if((*((*runseqdata_)[i]->indices_))[j] == pepIndx){
	//pepData = &(*((*runseqdata_)[i]->data_))[j];
	*pepData = (*((*runseqdata_)[i]->data_))[j];
	indx = 1;
	//cout << "found it with ratio: " << pepData->pepRatio[0] << " from " << (*((*runseqdata_)[i]->data_))[j].pepRatio[0] << endl;
      }
    }
  } // for (i = 0; i < dreamNum; ++i){

  if(! indx) {
    cerr << "error, could not find pepdatastrct for pepindex " << pepIndx << " and xmlindex " << xmlIndx << endl;
    exit(1);
  }

}


void ASAPProteinRatio::getDataRatio(double *ratio, double *error, double *inv_ratio, double *inv_error, double confL, 
				    double *data, double *dataErrs, double *inv_data, double *inv_dataErrs, 
				    double *dataWghs, int *dataIndx, int dataSize, int testType)
{

  int counts[4] = {0, 0, 0, 0};
  double *dtRatios, *dtErrors, *dt_inv_Ratios, *dt_inv_Errors, *dtWeights;
  int *dtIndx;
  int pass;
  int count, vldNum;
  double sum, inv_sum;
  double tmpError, tmp_inv_Error;
  double acc = 0.01;
  int i, j;

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
	    dataIndx[i] = 1;
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
	    dataIndx[i] = 1;
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
	    dataIndx[i] = 1;
	  else
	    dataIndx[i] = 0;
      }
    }
    //cout << "bailing early" << endl;
    return;
  }

  // allocate memory
  dtRatios = (double *) calloc(dataSize, sizeof(double));
  dtErrors = (double *) calloc(dataSize, sizeof(double));
  dt_inv_Ratios = (double *) calloc(dataSize, sizeof(double));
  dt_inv_Errors = (double *) calloc(dataSize, sizeof(double));
  dtWeights = (double *) calloc(dataSize, sizeof(double));
  dtIndx = (int *) calloc(dataSize, sizeof(int));


  // collect valid data, transform into log(ratio)
  for (i = 0; i < dataSize; ++i) 
    if(data[i] > 0.)
      dataIndx[i] = 0;
    else
      dataIndx[i] = 1;

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
  }
  else if(count == 1) { // only one valid data
    *ratio = dtRatios[0];
    *error = dtErrors[0];
    *inv_ratio = dt_inv_Ratios[0];
    *inv_error = dt_inv_Errors[0];
  }
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
    }
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
  }
  
  // free memory
  free(dtRatios);
  free(dtErrors);
  free(dt_inv_Ratios);
  free(dt_inv_Errors);
  free(dtWeights);
  free(dtIndx);
  //cout << "returning " << ratio[0] << " +- " << ratio[1] << endl;
}


// For a set of data and weight, this function finds the mean and standard deviation.
void ASAPProteinRatio::findMeanAndStdDevWeight(double *mean, double *error,
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
    }
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
void ASAPProteinRatio::DixonTest(double *data, int *outliers, int size)
{
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
    }
    else {
      if(ratio2 > cutoff) { // an outlier
	outliers[dataIndx[endIndx-1]] = 1;
	--size;
	--endIndx;
      }
    }
  }

  free(dataIndx);
  return;
}


// This function returns the value of Pade Approximation.
double ASAPProteinRatio::PadeApprx(double x, double *xa, double *ya, int size)
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


proDataStrct* ASAPProteinRatio::getProDataStruct() {
  if(! set_)
    computeProteinRatio();
  //  cout << pro_->ratio[0] << " +- " << pro_->ratio[1] << endl;

  return pro_;
}

int MSMSSeqDataCmp(void const *a, void const *b) {
  ASAPRatioMSMSSeqData** dd1 = (ASAPRatioMSMSSeqData**)a;
  ASAPRatioMSMSSeqData** dd2 = (ASAPRatioMSMSSeqData**)b;

  int result = strcmp((*dd1)->lightsequence_, (*dd2)->lightsequence_);
  return result?result:((char *)a-(char *)b); // stabilize sort
}


ASAPRatioMSMSSeqData** ASAPProteinRatio::getSortedMSMSSeqData() {
  ASAPRatioMSMSSeqData** output = new ASAPRatioMSMSSeqData* [runseqdata_->length()];
  for(int k = 0; k < runseqdata_->length(); k++)
    output[k] = (*runseqdata_)[k];

  //qsort(prtnAAs, *prtnAANum, sizeof(pairStrct), pairStrctCmp);

  qsort(output, runseqdata_->length(), sizeof(ASAPRatioMSMSSeqData*), MSMSSeqDataCmp);
  return output;
}

void ASAPProteinRatio::sortMSMSSeqData() {
  ASAPRatioMSMSSeqData** sorted = getSortedMSMSSeqData();
  //return;
  for(int k = 0; k < runseqdata_->length(); k++)
    runseqdata_->replace(k, sorted[k]);
 
  delete[] sorted;
}

// This function coverts two strings into a list of amino acid partners. It returns NULL on failure.
// lightString: L1, L2, L3, ...
// heavyString: H1, H2, H3, ...
// Pairs (L1, H1), (L2, H2), (L3, H3), ..., must appear in order.
pairStrct* ASAPProteinRatio::collectPrtnAAStrct(int *prtnAANum, char *lightString, char *heavyString)
{
  pairStrct *prtnAAs;
  char **tmpInputs;
  char tmpString[100];
  int tmpNum;
  int i;

  // light
  cnvtUpper(lightString);
  tmpInputs = getStrSects(prtnAANum, lightString, ',');
  if(*prtnAANum < 1) {
    printf("Invalid input for LIGHT isotopes: \"%s\".\n", lightString);
    fflush(stdout);
    freeMtrx((void **)tmpInputs, *prtnAANum);
    return NULL;
  }

  prtnAAs = (pairStrct *) calloc(*prtnAANum, sizeof(pairStrct));
  for (i = 0; i < *prtnAANum; ++i) 
    strcpy(prtnAAs[i].prtnA, tmpInputs[i]);
  freeMtrx((void **)tmpInputs, *prtnAANum);

  // heavy
  cnvtUpper(heavyString);
  tmpInputs = getStrSects(&tmpNum, heavyString, ',');
  if(*prtnAANum != tmpNum) {
    printf("Number of LIGHT isotopes (%d) doesn't match that of HEAVY (%d). \n", 
	   *prtnAANum, tmpNum);
    fflush(stdout);
    free(prtnAAs);
    freeMtrx((void **)tmpInputs, tmpNum);
    return NULL;
  }
  for (i = 0; i < *prtnAANum; ++i) 
    strcpy(prtnAAs[i].prtnB, tmpInputs[i]);
  freeMtrx((void **)tmpInputs, *prtnAANum);

  // sort prtnAAs
  for (i = 0; i < *prtnAANum; ++i) {
    if(strlen(prtnAAs[i].prtnA) < strlen(prtnAAs[i].prtnB)) {
      strcpy(tmpString, prtnAAs[i].prtnA);
      strcpy(prtnAAs[i].prtnA, prtnAAs[i].prtnB);
      strcpy(prtnAAs[i].prtnB, tmpString);
    }
  }
  qsort(prtnAAs, *prtnAANum, sizeof(pairStrct), pairStrctCmp);

  return prtnAAs;
}

// This function gets consecutive sections of a string, separated by "sep".
char** ASAPProteinRatio::getStrSects(int *sectNum, char *string, char sep)
{
  int lngth = (int)strlen(string);
  char **sects;
  int startIndx, endIndx;
  int i;

  // sectNum
  *sectNum = 1;
  for (i = 0; i < lngth; ++i) {
    if (string[i] == sep) {
      ++(*sectNum);
    }
  }

  // sects
  sects = (char **) calloc(*sectNum, sizeof(char *));
  *sectNum = 0;
  startIndx = 0;
  for (i = 0; i < lngth; ++i) {
    if (string[i] == sep) {
      endIndx = i;
      sects[*sectNum] = (char *) calloc(endIndx-startIndx+1, sizeof(char));
      strncpy(sects[*sectNum], string+startIndx, endIndx-startIndx);
      getRidOfSpace(sects[*sectNum]);
      ++(*sectNum);
      startIndx = endIndx + 1;
    }
  }
  endIndx = i;
  sects[*sectNum] = (char *) calloc(endIndx-startIndx+1, sizeof(char));
  strncpy(sects[*sectNum], string+startIndx, endIndx-startIndx);
  getRidOfSpace(sects[*sectNum]);
  ++(*sectNum);
  
  return sects;
}

// This function frees a matrix.
void ASAPProteinRatio::freeMtrx(void **mtrx, int size)
{
  int i;
  
  for (i = 0; i < size; ++i)
    free(mtrx[i]);
  free(mtrx);

  return;
}

// This function gets rid of any space at the beginning or end of a string.
void ASAPProteinRatio::getRidOfSpace(char *string) 
{
  int lngth = (int)strlen(string);
  int i;

  // start
  for (i = 0; i < lngth; ++i) {
    if (isspace(string[i]) == 0) // not a space
      break;
  }
  strcpy(string, string+i);

  // end
  lngth = (int)strlen(string);
  for (i = lngth-1; i >= 0; --i) {
    if (isspace(string[i]) == 0) // not a space
      break;
  }
  string[i+1] = '\0';

  return;
}


char* ASAPProteinRatio::getLightSequence(const char* sequence, pairStrct *prtnAAs, int prtnAANum) {
  char* output = new char[strlen(sequence)+1]; // max length
  int index = 0;
  for(int k = 0; sequence[k]; k++) {
    Boolean found = False; // until proven
    for(int j = 0; j < prtnAANum; j++) {
      if(k + strlen(prtnAAs[j].prtnA) <= strlen(sequence)) {
	Boolean match = True; // until proven
	for(int i = 0; prtnAAs[j].prtnA[i]; i++)
	  if(sequence[k+i] != prtnAAs[j].prtnA[i]) {
	    match = False;
	    break;
	  }
	if(match) {
	  for(int i = 0; prtnAAs[j].prtnB[i]; i++)
	    output[index++] = prtnAAs[j].prtnB[i];
	  found = True;
	  // increment
	  k += (int)strlen(prtnAAs[j].prtnA) - 1;
	  j = prtnAANum;
	}
      } // next position of this havy label
    } // next pair heavy j
    if(! found)
      output[index++] = sequence[k];
  } // next position
  output[index] = 0;
  return output;
}
