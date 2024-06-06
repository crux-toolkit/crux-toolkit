/*
Program       : QuanticProteinRatio
Author        : David Shteynberg
Date          : 11.03.2020


Copyright (C) 2020 David Shteynberg

Based on ASAPRatioProteinRatio
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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
*/

#include "QuanticProteinRatio.h"

QuanticProteinRatio::QuanticProteinRatio(double minprob, double minwt) {
  runseqdata_ = new vector<QuanticMSMSSeqData*>;
  set_ = False;
  start_ = 0;
  pro_ = NULL;
 
  min_probability_ = minprob;
  min_weight_ = minwt;
}

QuanticProteinRatio::~QuanticProteinRatio() {
  runseqdata_->clear();
  delete runseqdata_;
  delete pro_;

}

void QuanticProteinRatio::setRunInfo() {
   start_ = runseqdata_->size();
}

void QuanticProteinRatio::enter(const char* mod_pep, const ModificationInfo* modinfo, QuanticPSMData * data, int index, int xml_index, double wt, double prob, int msms_run_idx, int chg) {
  for(int k = start_; k < runseqdata_->size(); k++) {
    bool found = false;

    if(!strcmp(mod_pep, (*runseqdata_)[k]->modseq_) && (*runseqdata_)[k]->mod_info_->equivalentModification(modinfo, MOD_ERROR, mod_pep, NULL) ) {
      found = True;
    }
    
    
    if(found && (*runseqdata_)[k]->chg_ == chg &&
       (*runseqdata_)[k]->xml_index_ == xml_index &&
       (*runseqdata_)[k]->msms_run_idx_ == msms_run_idx) {
      
      (*runseqdata_)[k]->data_->push_back(data);
      (*runseqdata_)[k]->indices_->push_back(index);
      (*runseqdata_)[k]->probs_->push_back(prob);
      //cout << "adding index " << index << " to peptide " << light <<  " (" << (*runseqdata_)[k]->lightsequence_ << ")" << endl;
      return;
    }
  }

  // still here
  //cout << "making new entry with index " << index << " for peptide " << light << endl;
  runseqdata_->push_back(new QuanticMSMSSeqData(mod_pep, modinfo, data, index, xml_index, wt, prob, msms_run_idx, chg));

}



// the guts
void QuanticProteinRatio::computeProteinQuants() {
  sortMSMSSeqData();
  //initProDataStrct();
  
  initProQuantStrct(); 
  getProQuantStrct(pro_);


  set_ = True;
}


// This function initializes a proDataStrct *protein.
//void QuanticProteinRatio::initProDataStrct() {

void QuanticProteinRatio::initProQuantStrct() {
  int seqNum;
  int startIndx, endIndx;
  int i;

  // initialize protein
  
  pro_ = new proQuantStrct();//  pro_ = new proDataStrct();
  
  pro_->indx_ = -1;

  


  
  seqNum = 0;
  endIndx = 0;
  while(endIndx < runseqdata_->size()){
    startIndx = endIndx;
    while(endIndx < runseqdata_->size() 
	  && strcmp((*runseqdata_)[endIndx]->modseq_, (*runseqdata_)[startIndx]->modseq_) == 0 &&
	  (*runseqdata_)[endIndx]->chg_  == (*runseqdata_)[startIndx]->chg_ ) 
      ++endIndx;

    pro_->peps_.push_back(new pepQuantStrct());
    pro_->peps_[pro_->peps_.size()-1]->dataInc_ = true;    
    //    initSeqDataStrct(pro_->sequences+seqNum, startIndx, endIndx, min_probability_);

    
    if (pro_->peps_[seqNum]->indx_ != -1)
      pro_->indx_ = 0;
      

    pro_->peps_[seqNum]->chg_ = (*runseqdata_)[startIndx]->chg_;
    pro_->peps_[seqNum]->weight_ = 0; //(*runseqdata_)[startIndx]->weight_;

    initPepQuantStrct(pro_->peps_[seqNum], startIndx, endIndx, min_probability_);

    if(pro_->peps_[seqNum]->weight_ <= min_weight_) {
      pro_->peps_[seqNum]->dataInc_ = false;
    }
    else {
      pro_->peps_[seqNum]->weight_ = 0;
      for (i = 0; i < pro_->peps_[seqNum]->psms_.size(); ++i){
	if (pro_->peps_[seqNum]->psms_[i]->dataInc_) {
	  pro_->peps_[seqNum]->dataInc_ = true;
	  //break;
	}
	if (pro_->peps_[seqNum]->psms_[i]->weight_ >= 0.) {
	  pro_->peps_[seqNum]->weight_ += pro_->peps_[seqNum]->psms_[i]->weight_;
	}
      }
    }


    ++seqNum;
  }

}


// This function get the initial seqDataStrct of a unique peptide sequence.
//void QuanticProteinRatio::initSeqDataStrct(seqDataStrct *seqData, int start_index, int end_index, double mnProb)
void QuanticProteinRatio::initPepQuantStrct(pepQuantStrct* seqData, int start_index, int end_index, double mnProb)
{
  int totPepNum, totPepCount;
  int **pepIndx;
  int peakNum, peakIndx;
  double timeRange[2][2];
  double *weight;
  int i, j;


  int lastChg = -1;
  
  // total peptide number
  totPepNum = 0;
  for (i = start_index; i < end_index; ++i){
    totPepNum += (*runseqdata_)[i]->data_->size();
  }  

  int dreamNum = end_index - start_index;

  // get seqData->peaks

  // memory
  //  seqData->peaks = (dataStrct *) calloc(totPepNum, sizeof(dataStrct));
  for (i = 0; i < dreamNum; i++) {
    seqData->psms_.push_back(new psmQuantStrct());
  }

  weight = (double *) calloc(totPepNum, sizeof(double));

  // pepIndx
  pepIndx = (int **) calloc(dreamNum, sizeof(int *));
  for (i = 0; i < dreamNum; ++i)
    pepIndx[i] = (int *) calloc((*runseqdata_)[start_index + i]->data_->size(), sizeof(int));
  for (i = 0; i < dreamNum; ++i){
    for (j = 0; j < (*runseqdata_)[start_index + i]->data_->size(); ++j){
      pepIndx[i][j] = 0;
    }
  }

  // collect data
  peakNum = 0;
  totPepCount = 0;
  while(totPepCount < totPepNum){

    if (peakNum >= seqData->psms_.size()) {
      seqData->psms_.push_back(new psmQuantStrct());
    }

    seqData->psms_[peakNum]->indx_ = 0.;
    seqData->psms_[peakNum]->weight_ = 0.;
    seqData->psms_[peakNum]->msms_run_idx_ = -1;
    
    weight[peakNum] = 0.;
    peakIndx = 1;
    for (i = 0;  i < dreamNum; ++i){
      if(i == 0) {
	seqData->pep_ = string((*runseqdata_)[start_index + i]->modseq_);
      }

   
      
      for (j = 0; j < (*runseqdata_)[start_index + i]->data_->size(); ++j){
	if(pepIndx[i][j] == 1)
	  continue;


	if (seqData->chg_ != (*((*runseqdata_)[start_index + i]->data_))[j]->chg_) {
	  peakNum++;
	}

	if(    (*((*runseqdata_)[start_index + i]->data_))[j]->indx_ == -1){ // invalid data
	  
	  seqData->psms_[peakNum]->indx_ = -1;
	  seqData->psms_[peakNum]->weight_ = 0.;
	  seqData->psms_[peakNum]->msms_run_idx_ = -1;


	  weight[peakNum] = (*runseqdata_)[start_index + i]->weight_;
	  pepIndx[i][j] = 1;
	  peakIndx = 0;
	  ++totPepCount;
	  break;
	}
	else {
	  seqData->psms_[peakNum]->msms_run_idx_ = (*runseqdata_)[start_index + i]->msms_run_idx_;
	  
	}
	
	
	if (// not yet analyzed
	    pepIndx[i][j] == 0 
	    // from the same run
	    && (*runseqdata_)[start_index + i]->msms_run_idx_ == seqData->psms_[peakNum]->msms_run_idx_ ) {
	  
	  seqData->psms_[peakNum]->indx_ = 0;
	  
	  
	  if((*((*runseqdata_)[start_index + i]->probs_))[j] >= mnProb)
	    seqData->psms_[peakNum]->dataInc_ = true;
	  else
	    seqData->psms_[peakNum]->dataInc_ = false;


	  seqData->psms_[peakNum]->specName_ = (*((*runseqdata_)[start_index + i]->data_))[j]->specName_;
	  seqData->psms_[peakNum]->chg_ = (*((*runseqdata_)[start_index + i]->data_))[j]->chg_;
	  seqData->psms_[peakNum]->scan_ = (*((*runseqdata_)[start_index + i]->data_))[j]->scan_;
	
	  seqData->psms_[peakNum]->weight_ =
	    seqData->psms_[peakNum]->weight_ > (*((*runseqdata_)[start_index + i]->data_))[j]->getTotalQuant() ? 
	    seqData->psms_[peakNum]->weight_ : (*((*runseqdata_)[start_index + i]->data_))[j]->getTotalQuant();	  

	  weight[peakNum] += seqData->psms_[peakNum]->weight_;
	  //weight[peakNum] > (*runseqdata_)[start_index + i]->weight_ ?
	  //  weight[peakNum] : (*runseqdata_)[start_index + i]->weight_;

	  pepIndx[i][j] = 1;
	  ++totPepCount;
	  

	  
	}

      }

    }
    
    
    ++peakNum;
  
  }
  
  
  for (i = 0; i < dreamNum; ++i)
    free(pepIndx[i]);
  free(pepIndx);
  
  // get seqdata
  seqData->indx_ = -1;

  //  seqData->dataNum_ = peakNums;
  //seqData->dataCnts = (int *) calloc(peakNum, sizeof(int));
  seqData->weight_ = 0.;
  for (i = 0; i < peakNum; ++i) {
   if(seqData->psms_[i]->indx_ != -1)
      seqData->indx_ = 0;

      if (seqData->psms_[i]->dataInc_)
      seqData->dataInc_ = true;
    else
      seqData->dataInc_ = false;

     seqData->weight_ = seqData->weight_ > weight[i] ? seqData->weight_ : weight[i];

   seqData->weight_ =  seqData->weight_ > weight[i] ? seqData->weight_ : weight[i];
  } //   for (i = 0; i < peakNum; ++i) {
  free(weight);

}

// This function evaluates the ratio of a protein.
//void QuanticProteinRatio::getProDataStrct(proDataStrct *data)
void QuanticProteinRatio::getProQuantStrct(proQuantStrct *data)
{

  vector<vector<double>*> quants;
  vector<vector<double>*> stdvs;
  vector<double> weights;

  vector<double> errors;
  int *outliers, *seqIndx;
  int seqNum, tmpNum;
  double cmnErr;

  int num[2] = {0, 0};

  //  vector<double> quants;
  //vector<double> stdvs;

  int i, j;

  if (data->indx_ == -1) {
    for (i=0; i<data->quantmeans_.size(); i++) {
      data->quantmeans_[i] = 1 / data->quantmeans_.size();  //flat priors
      data->quantstdvs_[i] = -999.;
    }

    return;
  }

  // collect all sequences
  seqIndx = (int *) calloc(data->peps_.size(), sizeof(int));
  seqNum = 0;
  for (i = 0; i < data->peps_.size(); ++i) {
    getSeqQuantStrct(data->peps_[i]);
    // double check dataCnts
    if(data->peps_[i]->dataInc_) {
      data->peps_[i]->dataInc_ = false;
      for (j = 0; j < data->peps_[i]->psms_.size(); ++j){
	if(data->peps_[i]->psms_[j]->dataInc_ == 1){
	  data->peps_[i]->dataInc_ = true;
	  break;
	}
      }
    }
    if(data->peps_[i]->dataInc_ == 1 && data->peps_[i]->indx_ != -1) {
      seqIndx[i] = 1;
      ++seqNum;
    }
    else
      seqIndx[i] = 0;
  }

  // only 0 or 1 valid sequences
  if(seqNum < 2) {
    if(seqNum < 1) {
      for (i=0; i<data->quantmeans_.size(); i++) {
	data->quantmeans_[i] = 1 / data->quantmeans_.size();  //flat priors
	data->quantstdvs_[i] = -999.;
      } 
    } // if(seqNum < 1) {
    else {
      for (i = 0; i < data->peps_.size(); ++i) {
	if(seqIndx[i] == 1) {
	  for (int ii = 0; ii < data->peps_[i]->quantmeans_.size(); ii++) {
	    if (ii >= data->quantmeans_.size()) {
	      data->quantmeans_.push_back(0.);
	      data->quantstdvs_.push_back(0.);
	    }
	    data->quantmeans_[ii] =  data->peps_[i]->quantmeans_[ii];
	    data->quantstdvs_[ii] =  data->peps_[i]->quantstdvs_[ii];
	  }

	  break;
	}
      }
    }

    if(data->indx_ == 0)
      data->indx_ = 1;
    
    free(seqIndx);

    return;
  }


  // allocate memory
  //  ratios = (double *) calloc(seqNum, sizeof(double));
  //  errors = (double *) calloc(seqNum, sizeof(double));
  //  inv_ratios = (double *) calloc(seqNum, sizeof(double));
  //  inv_errors = (double *) calloc(seqNum, sizeof(double));
  //  weights = (double *) calloc(seqNum, sizeof(double));
  //  outliers = (int *) calloc(seqNum, sizeof(int));

  // get data
  cmnErr = 0.; //total error
  tmpNum = 0;
  for(i = 0; i < data->peps_.size(); ++i) {
    quants.push_back(new vector<double>);
    stdvs.push_back(new vector<double>);
    weights.push_back(0.);
    errors.push_back(0.);
    if(seqIndx[i] == 1) {
      weights[i] = data->peps_[i]->weight_;
      for (int ii = 0; ii < data->peps_[i]->quantmeans_.size(); ii++) {
	quants[i]->push_back(data->peps_[i]->quantmeans_[ii]);
	stdvs[i]->push_back(data->peps_[i]->quantstdvs_[ii]);

	cmnErr += data->peps_[i]->quantstdvs_[ii];
	errors[i] += data->peps_[i]->quantstdvs_[ii];
	++tmpNum;
      }


    }
  }
  cmnErr /= tmpNum;

  // get weight
  // for(i = 0; i < seqNum; ++i) {
  //   if(errors[i] < cmnErr)
  //     weights[i] = 1./cmnErr/cmnErr;
  //   else
  //     weights[i] = 1./errors[i]/errors[i];
  //   if (isinf(weights[i])) {
  //     weights[i] = data->peps_[i]->weight_;
  //   }
  //   else {
  //     data->peps_[i]->weight_ = weights[i];
  //   }
  // }

  
 
  
  getDataQuants(&data->quantmeans_,
		&data->quantstdvs_,
		&quants, &stdvs, &weights);

  
  // reset indx
  if(data->indx_ == 0)
    data->indx_ = 1;

  // free memory
  
}

// This function calculates the ratio of a unique sequence.
//void QuanticProteinRatio::getSeqDataStrctRatio(seqDataStrct *data)
void QuanticProteinRatio::getSeqQuantStrct(pepQuantStrct *data)
{
  // double *ratios, *errors, *inv_ratios, *inv_errors, *weights;

  //int *outliers, *peakIndx;
  int peakNum, tmpNum;

  vector<vector<double>*> quants;
  vector<vector<double>*> stdvs;
  vector<double> weights;

  vector<double> errors;
  vector<int> peakIndx;

  peakIndx.reserve(data->psms_.size());

  double cmnErr = 0.;
  int i, j;

  int runseq_idx = 0;

  for (runseq_idx  = 0; runseq_idx  < runseqdata_->size(); runseq_idx++) {
    if (data->chg_ == 0 || (data->chg_ == (*runseqdata_)[runseq_idx]->chg_))
      if ( !strcmp((*runseqdata_)[runseq_idx]->modseq_, data->pep_.c_str()) ) break;
  } 

  data->psms_.clear();  data->weight_ = 0;
  
  for (i=0; i < (*(*(*runseqdata_)[runseq_idx]).data_).size(); i++) {
    if (data->chg_ == 0 || (data->chg_ == (*(*(*runseqdata_)[runseq_idx]).data_)[i]->chg_)) {
      
      if (i >= data->psms_.size()) {
	  data->psms_.push_back(new psmQuantStrct());
	  data->psms_[data->psms_.size()-1]->specName_ = (*(*(*runseqdata_)[runseq_idx]).data_)[i]->specName_;
	  data->psms_[data->psms_.size()-1]->chg_ = (*(*(*runseqdata_)[runseq_idx]).data_)[i]->chg_;
	  data->psms_[data->psms_.size()-1]->indx_ = (*(*(*runseqdata_)[runseq_idx]).data_)[i]->indx_;
	  data->psms_[data->psms_.size()-1]->dataInc_ = true ;
	  data->psms_[data->psms_.size()-1]->weight_ =
	    data->psms_[data->psms_.size()-1]->weight_ > (*((*runseqdata_)[runseq_idx]->data_))[i]->getTotalQuant() ? 
	    data->psms_[data->psms_.size()-1]->weight_ : (*((*runseqdata_)[runseq_idx]->data_))[i]->getTotalQuant();	  
	  data->weight_ += data->psms_[data->psms_.size()-1]->weight_;
	}



	for (int ii=0; ii< (*((*runseqdata_)[runseq_idx]->data_))[i]->quantmeans_->size(); ii++) {
	  data->psms_[data->psms_.size()-1]->quantmeans_.push_back((*(*(*(*runseqdata_)[runseq_idx]->data_)[i]).quantmeans_)[ii]);
	  data->psms_[data->psms_.size()-1]->quantstdvs_.push_back((*(*(*(*runseqdata_)[runseq_idx]->data_)[i]).quantstdvs_)[ii]);
	}
    }
  }
  
  if (data->indx_ == -1) {
    for (i=0; i<data->quantmeans_.size(); i++) {
      data->quantmeans_[i] = 1 / data->quantmeans_.size();  //flat priors
      data->quantstdvs_[i] = -999.;
    }
    return;
  }

  // // collect all peaks
  // // peakIndx = (int *) calloc(data->psms_.size(), sizeof(int));
  // peakNum = 0;
  // for (i = 0; i < data->psms_.size(); ++i) {
  //   //getDataStrctRatio(&(data->peaks[i]));

  //   //DDS: SHOULD BE PARSED NOT CALCULATED
  //   //getPSMQuantStrct(data->psms_[i]); 

   
  //   if(data->psms_[i]->dataInc_ == true && data->psms_[i]->indx_ != -1) {
  //     peakIndx[i] = 1;
  //     ++peakNum;
  //   }
  //   else
  //     peakIndx[i] = 0;
  // }

  // // only 0 or 1 valid peaks
  // if(peakNum < 2) {
  //   if(peakNum < 1) {

  //     for (i=0; i<data->quantmeans_.size(); i++) {
  // 	data->quantmeans_[i] = 1 / data->quantmeans_.size();  //flat priors
  // 	data->quantstdvs_[i] = -999.;
  //     }

  //   }
  //   else {
  //     for (i = 0; i < data->psms_.size(); ++i) {
  // 	if(peakIndx[i] == 1) {
  // 	  for(int ii = 0; ii < data->psms_[i]->quantmeans_.size(); ++ii) {
  // 	    data->quantmeans_[ii] = data->psms_[i]->quantmeans_[ii];
  // 	    data->quantstdvs_[ii] = data->psms_[i]->quantstdvs_[ii];

  // 	  }
  // 	  break;
  // 	}
  //     }
  //   }

  //   if(data->indx_ == 0)
  //     data->indx_ = 1;

  //   //free(peakIndx);
  //   return;
  // }


  
  //// ************************************ DDS ****************** WORKING HERE !!!

  // allocate memory
  //ratios = (double *) calloc(peakNum, sizeof(double));
  //errors = (double *) calloc(peakNum, sizeof(double));
  //inv_ratios = (double *) calloc(peakNum, sizeof(double));
  //inv_errors = (double *) calloc(peakNum, sizeof(double));
  //weights = (double *) calloc(peakNum, sizeof(double));
  //outliers = (int *) calloc(peakNum, sizeof(int));

  // get data
  tmpNum = 0;
  for(i = 0; i < data->psms_.size(); ++i) {
    quants.push_back(new vector<double>);
    stdvs.push_back(new vector<double>);
    weights.push_back(0.);
    errors.push_back(0.);

    //if(peakIndx[i] == 1) {
    weights[i] = data->psms_[i]->weight_;
    for (int ii = 0; ii < data->psms_[i]->quantmeans_.size(); ii++) {
	quants[i]->push_back(data->psms_[i]->quantmeans_[ii]);
	stdvs[i]->push_back(data->psms_[i]->quantstdvs_[ii]);

	errors[i] += data->psms_[i]->quantstdvs_[ii];
      	cmnErr += data->psms_[i]->quantstdvs_[ii];
	++tmpNum;
      }
      
      //}
  }
  
  cmnErr /= tmpNum;

  // get weight
  // for(i = 0; i < data->psms_.size(); ++i) {
  //   if(errors[i] < cmnErr)
  //     weights[i] = 1./cmnErr/cmnErr;
  //   else
  //     weights[i] = 1./errors[i]/errors[i];

  //   if (isinf(weights[i])) {
  //     weights[i] = data->psms_[i]->weight_;
  //   }
  //   else {
  //     data->psms_[i]->weight_ = weights[i];
  //   }
  // }

  // calculate ratio and error
  //  if(data->indx_ == 0) {
  //  getDataQuants(&(data->quantmeans_), &(data->quantstdvs_),
  //		 quants, stdvs, weights);
  //    tmpNum = 0;
  //    for(i = 0; i < data->dataNum; ++i) {
  //      if(peakIndx[i] == 1) {
  //	data->dataCnts[i] = 1 - outliers[tmpNum];
  //	++tmpNum;
  //      }
  //    }      
  //  }
  //  else
  
  //TODO: DDS DETERMINE OUTLIERS

  getDataQuants(&(data->quantmeans_), &(data->quantstdvs_),
		&quants, &stdvs, &weights);
  
  // reset indx
  if(data->indx_ == 0)
    data->indx_ = 1;


  for (int i=0; i<quants.size(); i++) {
    quants[i]->clear();    delete quants[i];
    stdvs[i]->clear();    delete stdvs[i];
  }

}


// This function calculates the ratio of a unique peak.
//void QuanticProteinRatio::getDataStrctRatio(dataStrct *data)
// void QuanticProteinRatio::getPSMQuantStrct(psmQuantStrct *data) 
// {
//   //pepDataStrct *peptides;
//   //pepQuantStrct *peptides;
//   double *quants, *stdev, *weights;
//   int *outliers, *pepIndx;
//   int pepNum, tmpNum;
//   double tol = 1.e-2;

//   int i;

//   vector<vector<double>*> quants;
//   vector<vector<double>*> stdvs;
//   vector<double> weights;

//   vector<double> errors;
  
//   if (data->indx == -1) {
//     for (i=0; i<data->quantmeans_.size(); i++) {
//       data->quantmeans_[i] = 1 / data->quantmeans_.size();  //flat priors
//       data->quantstdvs_[i] = -999.;
//     }

    
//     return;
//   }

//   // collect all valid peptides

//   // read peptides
//   //peptides = (pepDataStrct *) calloc(data->dataNum, sizeof(pepDataStrct));
//   // peptides = (pepQuantStrct *) calloc(data->dataIndx_.size(), sizeof(pepQuantStrct));
 
//   // calculate quantities
//   // only 0 or 1 valid peptides
//   if(pepNum < 2) {
//     if(pepNum < 1) {
//       for (i=0; i<data->quantmeans_.size(); i++) {
// 	data->quantmeans_[i] = 1 / data->quantmeans_.size();  //flat priors
// 	data->quantstdvs_[i] = -999.;
//       }
//       data->weight_ = 0.;
//     }
//     else {

//       for (i=0; i<data->quantmeans_.size(); i++) {
// 	data->quantmeans_[i] = peptides[0]->quantmeans_[i];  
// 	data->quantstdvs_[i] = peptides[0]->quantstdvs_[i];
	
//       }
//        data->weight = peptides[0]->weight;
//     }

//     if(data->indx_ == 0) {
//       for (i = 0; i < data->dataIndx_.size(); ++i) {
// 	if(data->dataInc_[i] == 1 && pepIndx[i] == 0)  
// 	  data->dataInc_[i] = 0;
//       }
//       data->indx_ = 1;
//     }

//     free(pepIndx);
//     //free(peptides);
//     return;
//   }


  
//   // weight
//   data->weight = 0.;
//   for (i = 0; i < pepNum; ++i) {
//     if(outliers[i] == 0 && weights[i] > data->weight)
//       data->weight = weights[i];
//   }

//   // if identical ratios
//   tmpNum = 0;
//   for (i = 0; i < pepNum; ++i) {
//     if(outliers[i] == 0 
//        && fabs(ratios[i]-data->ratio[0]) > tol*data->ratio[0]) {
//       ++tmpNum;
//     }
//   }
//   if(tmpNum == 0) {
//     data->ratio[1] = 0.;
//     data->inv_ratio[1] = 0.;
//     tmpNum = 0;
//     for (i = 0; i < pepNum; ++i) {
//       if(outliers[i] == 0) {
// 	data->ratio[1] += errors[i];
// 	data->inv_ratio[1] += inv_errors[i];
// 	++tmpNum;
//       }
//     }
//     if(tmpNum > 0) {
//       data->ratio[1] /= tmpNum;
//       data->inv_ratio[1] /= tmpNum;
//     }
//     else {
//       for (i = 0; i < pepNum; ++i) {
// 	if(outliers[i] == 0) {
// 	  data->ratio[1] = errors[i];
// 	  data->inv_ratio[1] = inv_errors[i];
// 	  break;
// 	}
//       }
//     }
//   }
  
//   // reset indx
//   if(data->indx_ == 0)
//     data->indx_ = 1;

//   // free memory
//   free(peptides);
//   free(pepIndx);
//   free(ratios);
//   free(errors);
//   free(inv_ratios);
//   free(inv_errors);
//   free(weights);
//   free(outliers);

//   return;
// }


// This function searches pepDataStrct from dreamStrct.
//void QuanticProteinRatio::searchPepDataStrct(pepDataStrct *pepData, int pepIndx, int xmlIndx)
// void QuanticProteinRatio::searchPepQuantStrct(pepQuantStrct *pepData, int pepIndx, int xmlIndx)
// {
//   int indx = 0;
//   int i, j;

//   for (i = 0; indx == 0 && i < runseqdata_->size(); ++i){
//     if((*runseqdata_)[i]->xml_index_ != xmlIndx)
//       continue;
//     for (j = 0; indx == 0 && j < (*runseqdata_)[i]->data_->size(); ++j){
//       if((*((*runseqdata_)[i]->indices_))[j] == pepIndx){
// 	//pepData = &(*((*runseqdata_)[i]->data_))[j];
// 	*pepData = (*((*runseqdata_)[i]->data_))[j];
// 	indx = 1;
// 	//cout << "found it with ratio: " << pepData->pepRatio[0] << " from " << (*((*runseqdata_)[i]->data_))[j].pepRatio[0] << endl;
//       }
//     }
//   } // for (i = 0; i < dreamNum; ++i){

//   if(! indx) {
//     cerr << "error, could not find pepdatastrct for pepindex " << pepIndx << " and xmlindex " << xmlIndx << endl;
//     exit(1);
//   }

// }



// ************************** WORKING HERE >>> DDS !!!
//void QuanticProteinRatio::getDataRatio(double *ratio, double *error, double *inv_ratio, double *inv_error, double confL, 
//				    double *data, double *dataErrs, double *inv_data, double *inv_dataErrs, 
//				    double *dataWghs, int *dataIndx, int dataSize, int testType)
  
void QuanticProteinRatio::getDataQuants(vector<double>* quants, vector<double> *stdvs,
					vector<vector<double>*>* dataQuant,
					vector<vector<double>*>* dataStdev,
					vector<double> *dataWghs)
{
  vector<double> wts;
  vector<double> vals;
  
  double wtsum = 0;
  double wtsumsq = 0;
  double wtsqsum = 0;
  
  vector<double> sums;
  vector<double> sumsqrs;

  int SZ = 0;
  
  for (int n=0; n < dataQuant->size(); n++) {
    if ((*dataQuant)[n]->size() > SZ) {
      SZ = (*dataQuant)[n]->size();
    }
    for (int i=0; i<(*dataQuant)[n]->size(); i++) {
      if (i>=sums.size()) {
	sums.push_back(0);
      }
     
      sums[i] += (*dataWghs)[n] * (*(*dataQuant)[n])[i];
      
      if (i==0) {
	wtsum += (*dataWghs)[n];
	wtsqsum += (*dataWghs)[n]*(*dataWghs)[n];
      }
    }
  }

  
  if (dataQuant->size())
    for (int i=0; i<SZ; i++) {
      if (i>=quants->size()) {
	quants->push_back(0);
      }
      if (wtsum > 0.) 
	(*quants)[i] = sums[i] / wtsum;
      else
	(*quants)[i] = 0;
    }
  
  
  for (int n=0; n < dataQuant->size(); n++) {
    for (int i=0; i<(*dataQuant)[n]->size(); i++) {
      if (n==0) sums[i]=0;
      
      sums[i] += (*dataWghs)[n] *
	( (*(*dataQuant)[n])[i] - (*quants)[i] ) *
	( (*(*dataQuant)[n])[i] - (*quants)[i] );
      
    }
  }
  
  for (int i=0; i<SZ; i++) {
    if (i>=stdvs->size()) {
      stdvs->push_back(0);
    }
    
    if (fabs (wtsum*wtsum - wtsqsum)  > 0.) 
      (*stdvs)[i] = sqrt( sums[i] * wtsum / (wtsum*wtsum - wtsqsum) );
    else
      (*stdvs)[i] = 0;
  }
  
  return;
}



proQuantStrct* QuanticProteinRatio::getProQuantStrct() {
  if(! set_)
    computeProteinQuants();
  //  cout << pro_->ratio[0] << " +- " << pro_->ratio[1] << endl;

  return pro_;
}

bool compareSeqs(const QuanticMSMSSeqData* a, const QuanticMSMSSeqData* b) {
   int result = strcmp(a->modseq_, b->modseq_);
   return result > 0;
}



void QuanticProteinRatio::sortMSMSSeqData() {
  std::sort(runseqdata_->begin(),runseqdata_->end(), compareSeqs);
}

// This function frees a matrix.
void QuanticProteinRatio::freeMtrx(void **mtrx, int size)
{
  int i;
  
  for (i = 0; i < size; ++i)
    free(mtrx[i]);
  free(mtrx);

  return;
}

// This function gets rid of any space at the beginning or end of a string.
void QuanticProteinRatio::getRidOfSpace(char *string) 
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

