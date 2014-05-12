#include "PepRanker.h"
#include "modifications.h"

PepRanker::PepRanker() :  
  seed(0),
  selectionfdr(0.01),
  num_hu(4),
  mu(0.01),
  weightDecay(0.000),
  max_net_gen(NULL),
  max_net_targ(NULL),
  net_clones(NULL),
  psm_count(0),
  in_dir(""), 
  out_dir(""), 
  skip_cleanup_flag(0),
  overwrite_flag(0),
  fileroot(""),
  feature_file_flag(0),
  file_format_("")
{
}

PepRanker::~PepRanker()
{
  delete [] max_net_gen;
  delete [] max_net_targ;
  delete [] net_clones;
}

int PepRanker :: getOverFDRPSM(PSMScores &set, NeuralNet &n, double fdr)
{
  double *r;
  double* featVec;

  for(int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = n.fprop(featVec);
      set[i].score = r[0];

    }
  return set.calcOverFDR(fdr);
}

/***********************************************************/

double PepRanker :: get_peptide_score_xcorr(int pepind)
{
  int num_psm = d.pepind2num_psm(pepind);
  int *psminds = d.pepind2psminds(pepind);
  double max_sc = -100000000.0;
  int max_ind = 0;
  for(int i = 0; i < num_psm; i++)
    {
      int psmind = psminds[i];
      double *featVec = d.psmind2features(psmind);
      double sc = featVec[3];
      if(max_sc < sc)
	{
	  max_sc = sc;
	  max_ind = i;
	}
    }
  return max_sc;
}

void PepRanker :: getMultiFDRXCorr(PepScores &set, vector<double> &qvalues)
{
  int pepind = 0;
  for(int i = 0; i < set.size(); i++)
    {
      pepind = set[i].pepind;
      double sc = get_peptide_score_xcorr(pepind);
      set[i].score = sc;
    }
  
  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti);
}


double PepRanker :: get_peptide_score(int pepind, NeuralNet &n)
{
  int num_psm = d.pepind2num_psm(pepind);
  int *psminds = d.pepind2psminds(pepind);
  double max_sc = -100000000.0;
  int max_ind = 0;
  for(int i = 0; i < num_psm; i++)
    {
      int psmind = psminds[i];
      double *feat = d.psmind2features(psmind);
      double *sc = n.fprop(feat);
      if(max_sc < sc[0])
	{
	  max_sc = sc[0];
	  max_ind = i;
	}
    }
  if((int)pepind_to_max_psmind.size() == d.get_num_peptides())
    pepind_to_max_psmind[pepind] = psminds[max_ind];
  return max_sc;
}


int PepRanker :: getOverFDR(PepScores &s, NeuralNet &n,double fdr)
{
  int pepind = 0;
  for(int i = 0; i < s.size(); i++)
    {
      pepind = s[i].pepind;
      double sc = get_peptide_score(pepind,n);
      s[i].score = sc;
    }

  int overFDR = s.calcOverFDR(fdr);

  return overFDR;
}

void PepRanker :: getMultiFDR(PepScores &set, NeuralNet &n, vector<double> &qvalues)
{
  int pepind = 0;
  for(int i = 0; i < set.size(); i++)
    {
      pepind = set[i].pepind;
      double sc = get_peptide_score(pepind,n);
      set[i].score = sc;
    }
  
  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti);
}

void PepRanker :: write_max_nets(string filename, NeuralNet* max_net)
{
  ofstream f1(filename.c_str());
  f1 << "FDR thresh" << "\t" << "PSMs trn" << "\t" << "PSMs tst" << endl;
  for(int count = 0; count < num_qvals; count++)
    {
      net = max_net[count];
      int r = getOverFDR(testset,net, qvals[count]);
      int r1 = getOverFDR(trainset,net, qvals[count]);
      f1 << qvals[count] << "\t" << r1 << "\t" << r << "\n";
      
      double qn;
      if(qvals[count] < 0.01)
        qn = 0.0012;
      else
	qn = 0.005;
      r = getOverFDR(testset,net, qvals[count]+qn);
      r1 = getOverFDR(trainset,net, qvals[count]+qn);
      f1 << qvals[count]+qn << "\t" << r1 << "\t" << r << "\n";
      
    }
  f1.close();
}

/***************************************************************************/

double PepRanker :: get_protein_score(int pepind)
{
  int num_psms = d.pepind2num_psm(pepind);
  int *psminds = d.pepind2psminds(pepind);
  double max_sc = -1000000.0;
  int max_ind = 0;
  for (int j = 0; j < num_psms; j++)
    {
      double *feat = d.psmind2features(psminds[j]);
      double *sc = net_clones[psm_count].fprop(feat);
      if(sc[0] > max_sc)
	{
	  max_sc = sc[0];
	  max_ind = j;
	}
      psm_count++;
    }
  max_psm_inds.push_back(max_ind);
      
  return max_sc;
}

void PepRanker :: calc_gradients(int pepind, int label, int i)
{
  double *gc = new double[1];
  gc[0] = -label;
  int num_psms = d.pepind2num_psm(pepind);
  int clone_ind = psm_count+max_psm_inds[i];
  net_clones[clone_ind].bprop(gc);
  psm_count += num_psms;
  delete[] gc;
}



void PepRanker :: train_net_ranking(PepScores &set, int interval)
{
  double r1;
  double r2;
  double diff = 0;
  int label = 1;
  double *gc = new double[1];


  for(int i = 0; i < set.size(); i++)
    { 
      int ind1, ind2;
      int label_flag = 1;
      //get the first index
      if(interval == 0)
	ind1 = 0;
      else
	ind1 = myrandom_limit(interval);
      if(ind1>set.size()-1) continue;
      if(set[ind1].label == 1)
	label_flag = -1;
      else
	label_flag = 1;
      
      int cn = 0;
      while(1)
	{
	  ind2 = myrandom_limit(interval);
	  if(ind2>set.size()-1) continue;
	  if(set[ind2].label == label_flag) break;
	  if(cn > 1000)
	    {
	      ind2 = myrandom_limit(set.size());
	      break;
	    }
	  cn++;
	}
      
      //get the scores
      max_psm_inds.erase(max_psm_inds.begin(),max_psm_inds.end());
      psm_count = 0;
      //pass both through the net
      r1 = get_protein_score(set[ind1].pepind);
      assert(psm_count == d.pepind2num_psm(set[ind1].pepind));
      r2 = get_protein_score(set[ind2].pepind);;
      assert((int)max_psm_inds.size() == 2);
      
      
      diff = r1-r2;
      

      label=0;
      if(  set[ind1].label==1 && set[ind2].label==-1)
	label=1;
      if( set[ind1].label==-1 && set[ind2].label==1)
	    label=-1;
      
      if(label != 0)
	{
	  if(label*diff<1)
	    {
	      net.clear_gradients();
	      psm_count = 0;
	      calc_gradients(set[ind1].pepind, label, 0);
	      calc_gradients(set[ind2].pepind, -1.0*label, 1);
	      net.update(mu,weightDecay);
	    }
	  
	}
    }
  delete[] gc;
  
}

void PepRanker :: train_many_target_nets()
{

  int  thr_count = num_qvals-1;
  while (thr_count > 0)
    {
      net.copy(max_net_gen[thr_count]);
        
      carp(CARP_INFO, "training threshold %d", thr_count);
      interval = max_overFDR[thr_count];
      for(int i=switch_iter;i<niter;i++) {
	//sorts the examples in the training set according to the current net scores
	getMultiFDR(trainset,net,qvals);
	train_net_ranking(trainset, interval);
			
	for(int count = 0; count < num_qvals;count++)
	  {
	    if(overFDRmulti[count] > max_overFDR[count])
	      {
		max_overFDR[count] = overFDRmulti[count];
		max_net_targ[count] = net;
	      }
	  }

	if((i % 3) == 0)
	  {
	    carp(CARP_INFO, "Iteration %d :", i);
	    getMultiFDR(trainset,net,qvals);
	    carp(CARP_INFO, "trainset %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d ", 
		 qvals[0], overFDRmulti[0], qvals[1], overFDRmulti[1], qvals[2], overFDRmulti[2],
		 qvals[3], overFDRmulti[3], qvals[4], overFDRmulti[4], qvals[5], overFDRmulti[5],
		 qvals[6], overFDRmulti[6], qvals[7], overFDRmulti[7], qvals[8], overFDRmulti[8],
		 qvals[9], overFDRmulti[9], qvals[10], overFDRmulti[10], qvals[11], overFDRmulti[11],
		 qvals[12], overFDRmulti[12], qvals[13], overFDRmulti[13]);
	    getMultiFDR(testset,net,qvals);
	    carp(CARP_INFO, "testset %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d\n ", 
		 qvals[0], overFDRmulti[0], qvals[1], overFDRmulti[1], qvals[2], overFDRmulti[2],
		 qvals[3], overFDRmulti[3], qvals[4], overFDRmulti[4], qvals[5], overFDRmulti[5],
		 qvals[6], overFDRmulti[6], qvals[7], overFDRmulti[7], qvals[8], overFDRmulti[8],
		 qvals[9], overFDRmulti[9], qvals[10], overFDRmulti[10], qvals[11], overFDRmulti[11],
		 qvals[12], overFDRmulti[12], qvals[13], overFDRmulti[13]);

	  }

      }
      thr_count -= 3;
    }
}


void PepRanker :: train_many_general_nets()
{
  interval = trainset.size();
  for(int i=0;i<switch_iter;i++) {
    train_net_ranking(trainset, interval);
           
    //record the best result
    getMultiFDR(thresholdset,net,qvals);
    for(int count = 0; count < num_qvals;count++)
      {
	if(overFDRmulti[count] > max_overFDR[count])
	  {
	    max_overFDR[count] = overFDRmulti[count];
	    max_net_gen[count] = net;
	  }
      }
    if((i % 10) == 0)
      {
	carp(CARP_INFO, "Iteration %d :", i);
	getMultiFDR(trainset,net,qvals);
	carp(CARP_INFO, "trainset %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d ", 
	     qvals[0], overFDRmulti[0], qvals[1], overFDRmulti[1], qvals[2], overFDRmulti[2],
	     qvals[3], overFDRmulti[3], qvals[4], overFDRmulti[4], qvals[5], overFDRmulti[5],
	     qvals[6], overFDRmulti[6], qvals[7], overFDRmulti[7], qvals[8], overFDRmulti[8],
	     qvals[9], overFDRmulti[9], qvals[10], overFDRmulti[10], qvals[11], overFDRmulti[11],
	     qvals[12], overFDRmulti[12], qvals[13], overFDRmulti[13]);
	getMultiFDR(testset,net,qvals);
	carp(CARP_INFO, "testset %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d\n ", 
	     qvals[0], overFDRmulti[0], qvals[1], overFDRmulti[1], qvals[2], overFDRmulti[2],
	     qvals[3], overFDRmulti[3], qvals[4], overFDRmulti[4], qvals[5], overFDRmulti[5],
	     qvals[6], overFDRmulti[6], qvals[7], overFDRmulti[7], qvals[8], overFDRmulti[8],
	     qvals[9], overFDRmulti[9], qvals[10], overFDRmulti[10], qvals[11], overFDRmulti[11],
	     qvals[12], overFDRmulti[12], qvals[13], overFDRmulti[13]);
      }
  }

}

void PepRanker::train_many_nets()
{
  carp(CARP_INFO, "Before Iterating");
  getMultiFDRXCorr(trainset,qvals);
  carp(CARP_INFO, "trainset %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d ", 
       qvals[0], overFDRmulti[0], qvals[1], overFDRmulti[1], qvals[2], overFDRmulti[2],
       qvals[3], overFDRmulti[3], qvals[4], overFDRmulti[4], qvals[5], overFDRmulti[5],
       qvals[6], overFDRmulti[6], qvals[7], overFDRmulti[7], qvals[8], overFDRmulti[8],
       qvals[9], overFDRmulti[9], qvals[10], overFDRmulti[10], qvals[11], overFDRmulti[11],
       qvals[12], overFDRmulti[12], qvals[13], overFDRmulti[13]);
  getMultiFDRXCorr(testset,qvals);
  carp(CARP_INFO, "trainset %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d  %.2f:%d %.2f:%d %.2f:%d %.2f:%d %.2f:%d \n", 
       qvals[0], overFDRmulti[0], qvals[1], overFDRmulti[1], qvals[2], overFDRmulti[2],
       qvals[3], overFDRmulti[3], qvals[4], overFDRmulti[4], qvals[5], overFDRmulti[5],
       qvals[6], overFDRmulti[6], qvals[7], overFDRmulti[7], qvals[8], overFDRmulti[8],
       qvals[9], overFDRmulti[9], qvals[10], overFDRmulti[10], qvals[11], overFDRmulti[11],
       qvals[12], overFDRmulti[12], qvals[13], overFDRmulti[13]);

  train_many_general_nets();
  
  //copy the general net into target nets;
  for(int count = 0; count < num_qvals; count++)
    max_net_targ[count] = max_net_gen[count];
  
  train_many_target_nets();  
  
  ostringstream fname;
  fname << out_dir << "/" << fileroot << "pep-ranker.peps.at.fdr.thresholds.txt";;
  write_max_nets(fname.str(), max_net_targ);

}

void PepRanker :: setup_for_training()
{
  //mysrandom(seed); This is set by CruxApplication::initialize()
  carp(CARP_INFO, "reading data");
  
  ostringstream res;
  res << out_dir << "/qranker_output";
  res_prefix = res.str();
    
  d.load_data_pep_training();
  if(feature_file_flag)
    {
      string str = feature_file_name.str();
      if(!d.print_features(str))
	carp(CARP_INFO, "could not open file %s for writing features. Feature file will not be written", feature_file_name.str().c_str());
    }
  d.normalize_psms();
  d.load_labels_pep_training();
  PepScores::fillFeaturesSplit(trainset, testset, d, 0.75);
  thresholdset = trainset;
  d.clear_labels_pep_training();

  int max_psm = 0;
  for(int i = 0; i < trainset.size(); i++)
    {
      int pepind = trainset[i].pepind;
      int num_psm = d.pepind2num_psm(pepind);
      if(num_psm > max_psm)
	max_psm=num_psm;
    }

  switch_iter =30;
  niter = 40;
   
  num_qvals = 14;
  qvals.resize(num_qvals,0.0);
  overFDRmulti.resize(num_qvals,0);
  max_overFDR.resize(num_qvals,0); 

  //initialize q-value arrays
  double q = 0.0;
  for(int count = 0; count < num_qvals; count++)
    {
      qvals[count] = q;
      if(q < 0.01)
	q+=0.0025;
      else
	q+=0.01;
    }
  
  max_net_gen = new NeuralNet[num_qvals];
  max_net_targ = new NeuralNet[num_qvals];
  //set the linear flag: 1 if linear, 0 otherwise
  int lf = 0; num_hu = 3;
  if(num_hu == 1)
    lf = 1;
  //set whether there is bias in the linear units: 1 if yes, 0 otherwise
  int bs = 0;
  
  net.initialize(d.get_num_features(),num_hu,lf,bs);
  for(int count = 0; count < num_qvals; count++){
    max_net_gen[count] = net;
  }

  //create the net clones
  net_clones = new NeuralNet[2*max_psm];
  for (int i = 0; i < 2*max_psm;i++)
    net_clones[i].clone(net);
  
}


int PepRanker::run( ) {
 
  setup_for_training();
  train_many_nets();
  report_results_xml_tab();
  return 0;
}

/**************************************************************************/

void PepRanker :: get_pep_seq(string &pep, string &seq, string &n, string &c)
{
  string tmp;
  int pos;
  pos = pep.find('.');
  n = pep.at(pos-1); 
  tmp = pep.substr(pos+1, pep.size());

  pos = tmp.find('.');
  c = tmp.at(pos+1);
  seq = tmp.substr(0, pos);
}

void PepRanker :: get_tab_delim_proteins(string protein_str, vector<string> &proteins)
{
  proteins.clear(); 
  string str = protein_str;
  size_t pos = str.find("\t");
  while(pos != string::npos){
    if(pos == 0){
      str = str.substr(1,str.size()-1);
   } else{
      string prot = str.substr(0,pos);
      str = str.substr(pos+1,str.size()-1); 
      proteins.push_back(prot);
    }
    pos = str.find("\t");
  }
  proteins.push_back(str);
}


void PepRanker :: setup_for_reporting_results()
{
  d.load_labels_pep_training();
  trainset.clear();
  testset.clear();
  thresholdset.clear();
  PepScores::fillFeaturesFull(trainset, d);
  psmtrainset.clear();
  psmtestset.clear();
  PSMScores::fillFeaturesFull(psmtrainset, d);
  d.clear_labels_pep_training();
  
  //choose the best net for the selectionfdr
  int max_fdr = 0;
  int fdr = 0;
  int ind = 0;
  for(unsigned int count = 0; count < qvals.size();count++)
    {
      fdr = getOverFDR(trainset, max_net_targ[count], selectionfdr);
      if(fdr > max_fdr) {
        max_fdr = fdr;
        ind = count;
      }
    }
  net = max_net_targ[ind];
  pepind_to_max_psmind.clear();
  pepind_to_max_psmind.resize(d.get_num_peptides(),-1);
  int fdr_trn = getOverFDR(trainset,net, selectionfdr);
  carp(CARP_INFO, "total peptides at q<%.2f: %d", selectionfdr, fdr_trn);
  getOverFDRPSM(psmtrainset, net, selectionfdr);
  d.clear_data_pep_training();
  d.load_data_pep_results();

  computePepNSAF();
  computePEP();
  
}

/*****************************************************/

void PepRanker :: write_results_peptides_xml(ofstream &os)
{
  string protein_str;
  vector<string> tab_delim_proteins;
  os << "<peptides>" <<endl;
  int cn = 0;
  for(int i = 0; i < trainset.size(); i++)
    {
      int pepind = trainset[i].pepind;
      if(trainset[i].label == 1)
	{
	  cn++;
	  //write out proteins
	  string pep = d.ind2pep(pepind);
	  string seq, n, c;
	  get_pep_seq(pep, seq, n, c);
	  os << " <peptide peptide_id=\"" << seq << "\">" << endl;
	  os << "  <q_value>" << trainset[i].q << "</q_value>" << endl;
	  os << "  <score>" << trainset[i].score << "</score>" << endl;
	  os << "  <nsaf>" << trainset[i].nsaf << "</nsaf>" << endl;
          os << " <barista PEP>" <<trainset[i].PEP<<"</barista PEP>"<<endl; 
	  //write out peptides
	  int psmind = pepind_to_max_psmind[pepind];
	  if(psmind > -1)
	    os << "  <main_psm_id>"<< psmind << "</main_psm_id>" << endl;
	  else
	    cout << "warning: did not assign peptide" << pep  << " ind " << pepind << " max psmind\n";
	    	  
	  //print out all the psms in which this peptide is present
	  int num_psm = d.pepind2num_psm(pepind);
	  int *psminds = d.pepind2psminds(pepind);
	  os << "  <psm_ids>" << endl;
	  for(int j = 0; j < num_psm; j++)
	    os << "   <psm_id>" << psminds[j] << "</psm_id>" << endl;
	  os << "  </psm_ids>" << endl;
	  

	  //print out all the proteins
	  int num_prot = d.pepind2num_prot(pepind);
	  int *protinds = d.pepind2protinds(pepind);
	  os << "  <protein_ids>" << endl;
	  for(int j = 0; j < num_prot; j++)
	    {
	      protein_str = d.ind2prot(protinds[j]);
	      get_tab_delim_proteins(protein_str, tab_delim_proteins);
	      for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
		os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
	      //os << "   <protein_id>" << d.ind2prot(protinds[j]) << "</protein_id>" << endl;
	    }
	  
	  os << "  </protein_ids>" << endl;
	  os << " </peptide>" << endl;  
	}
      
    }
  os << "</peptides>" << endl;
}

void PepRanker :: write_results_psm_xml(ofstream &os)
{
  os << "<psms>" <<endl;
  int cn = 0;
  for(int i = 0; i < psmtrainset.size(); i++)
    {
      int psmind = psmtrainset[i].psmind;
      if(psmtrainset[i].label == 1)
	{
	  cn++;
	  //write out proteins
	  os << " <psm psm_id=" << "\"" << psmind << "\"" << ">" << endl;
	  os << "  <q_value>" << psmtrainset[i].q << "</q_value>" << endl;
	  os << "  <score>" << psmtrainset[i].score << "</score>" << endl;
	  os << "  <scan>" << d.psmind2scan(psmind) << "</scan>" << endl;
	  os << "  <charge>" << d.psmind2charge(psmind) << "</charge>" << endl;
	  os << "  <precursor_mass>" << d.psmind2precursor_mass(psmind) << "</precursor_mass>" << endl;
	  int pepind = d.psmind2pepind(psmind);
	  string pep = d.ind2pep(pepind);
	  string seq, n,c;
	  get_pep_seq(pep,seq,n,c);
	  os << "  <peptide_seq n =\"" << n << "\" c=\"" << c << "\" seq=\"" << seq << "\"/>" << endl;
	  os << "  <file_name>" << d.psmind2fname(psmind) << "</file_name>" << endl;
	  os << " </psm>" << endl;  
	}
    }
  os << "</psms>" << endl;
}


void PepRanker :: report_results_xml()
{

  ostringstream fname;
  fname << out_dir << "/" << fileroot << "barista.xml";
  ofstream of(fname.str().c_str());
  of << "<barista_output>" << endl;
  of << endl;
 
  write_results_peptides_xml(of);

  write_results_psm_xml(of);

  ostringstream xml_file_name;
  xml_file_name << out_dir << "/" << fileroot << "barista.target.pep.xml";
  PepXMLWriter xmlfile;
  xmlfile.openFile(xml_file_name.str().c_str(), overwrite_flag);

  //...
  xmlfile.closeFile();
  write_results_pep_xml(xmlfile);

  
  of << endl;
  of << "</barista_output>" << endl;
  of.close();
    
}

/***********************************************************************************/

void PepRanker :: write_results_peptides_tab(ofstream &os)
{
  int ps = pepind_to_max_psmind[trainset[0].pepind]; 
  os<< "q-value" << "\t" << "barista score" << "\t";
  os << "NSAF score" << "\t";
  os<<"barista PEP\t";
  os << "scan" << "\t" ;
  os<< "charge" << "\t";
  os<<"spectrum precursor m/z"<<"\t";
  os<<"spectrum neutral mass"<<"\t"; 
  os<<"peptide mass"<<"\t";
  os<< "delta_cn"<< "\t";
  if(d.psmind2spscore(ps)!=-1){
    os<<"sp score"<<"\t";
    os<<"sp rank\t";
  }
  os<<"xcorr score"<<"\t";
  os<<"xcorr rank"<<"\t";
  if(d.psmind2spscore(ps)!=-1){
    os<<"b/y ions matched"<<"\t";
    os<<"b/y ions total"<<"\t";
  }
  os<<"matches/spectrum"<<"\t";
  os<<"sequence\t"; 
  os<<"cleavage type"<<"\t"; 
  os<<"protein id"<<"\t";
  os<<"flanking_aa"<<endl;
  int cn = 0;
  for(int i = 0; i < trainset.size(); i++){
    if(trainset[i].label == 1){
      cn++;
      //write out proteins
      int pepind = trainset[i].pepind;
      string pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      //q-value
      os << trainset[i].q << "\t";
      //score 
      os << trainset[i].score << "\t";
      ///nsaf
      os<<trainset[i].nsaf<<"\t";
     //PEP
      os<<trainset[i].PEP<<"\t"; 
   
      //write out peptides
      int psmind = pepind_to_max_psmind[pepind];
      if(psmind > -1){
       //scan  
	os << d.psmind2scan(psmind) << "\t" ;
       //charge
        os<< d.psmind2charge(psmind)<<"\t";
        //mass-to-charge ratio 
        os<<(d.psmind2precursor_mass(psmind)+
        d.psmind2charge(psmind)*MASS_PROTON)/
        d.psmind2charge(psmind)<<"\t";
        //Spectrum Neutral Mass 
        os<<d.psmind2precursor_mass(psmind)<<"\t";
        //Peptide Mass
        os<<d.psmind2peptide_mass(psmind)<<"\t";
        //DELTA CN
        os <<d.psmind2deltaCn(psmind)<< "\t";
        if(d.psmind2spscore(ps)!=-1){
          //Sp Score
          os<<d.psmind2spscore(psmind)<<"\t";
          //Sp Rank 
          os<<d.psmind2sp_rank(psmind)<<"\t";
        }
        //xcorr Score
      	os<<d.psmind2xcorr(psmind)<<"\t";
      	//xcorr rank
      	os<<d.psmind2xcorr_rank(psmind)<<"\t";
        if(d.psmind2spscore(ps)!=-1){
          //by ions match 
          os<<d.psmind2by_ions_matched(psmind)<<"\t"; 
          //by ions total 
          os<<d.psmind2by_ions_total(psmind)<<"\t";
        }
      	//Matches/Spectrum 
      	os<<d.psmind2matches_spectrum(psmind)<<"\t";
        //sequence 
        os<<seq<<"\t"; 
        //cleavage type
        os<<cleavage_type<<"\t"; 
      	//protein id
      	vector<string> prots;  
        get_protein_id(pepind,prots);
        print_protein_ids(prots,os,psmind);
        //Flanking_aa 
      	os<<n<<c<<endl;     
      }else{
	  cout<< "waning: did not assign peptide max psmind\n";
          os<<"not assignd\t";
       }

    }
  }
}

void PepRanker :: write_results_psm_tab(ofstream &os)
{
  int ps= psmtrainset[0].psmind; 
  os << "scan" << "\t" << "charge" << "\t";
  os << "q-value" << "\t" << "barista score" << "\t";
  os << "barista PEP\t";
  os<<"spectrum precursor m/z"<<"\t";
  os<<"spectrum neutral mass"<<"\t"; 
  os<<"peptide mass"<<"\t";
  os<< "delta_cn"<< "\t";
  if(d.psmind2spscore(ps)!=-1){
    os<<"sp score"<<"\t";
    os<<"sp rank\t";
  }
  os<<"xcorr score"<<"\t";
  os<<"xcorr rank"<<"\t";
  if(d.psmind2spscore(ps)!=-1){
    os<<"b/y ions matched"<<"\t";
    os<<"b/y ions total"<<"\t";
  }
  os<<"matches/spectrum"<<"\t";
  os<<"sequence"<<"\t";
  os<<"cleavage type"<<"\t"; 
  os<<"protein id"<<"\t";
  os<<"flanking aa"<<"\t";
  os << "filename" << endl;
 int cn = 0;
  for(int i = 0; i < psmtrainset.size(); i++)
    {
      if(psmtrainset[i].label == 1)
	{
	  cn++;
	  //write out proteins
	  int psmind = psmtrainset[i].psmind; 
	  int pepind = d.psmind2pepind(psmind);
	  string pep = d.ind2pep(pepind);
	  string seq, n,c;
	  get_pep_seq(pep,seq,n,c);
	  //os << psmind << "\t";
	  os << d.psmind2scan(psmind) << "\t";
	  os << d.psmind2charge(psmind) << "\t";
	  os << psmtrainset[i].q << "\t";
	  os << psmtrainset[i].score << "\t";
          os << psmtrainset[i].PEP << "\t";
          //mass-to-charge ratio 
          os<<(d.psmind2precursor_mass(psmind)+
          d.psmind2charge(psmind)*MASS_PROTON)/
          d.psmind2charge(psmind)<<"\t";
          //Spectrum Neutral Mass 
          os<<d.psmind2precursor_mass(psmind)<<"\t";
          //Peptide Mass
          os<<d.psmind2peptide_mass(psmind)<<"\t";
          //DELTA CN
          os <<d.psmind2deltaCn(psmind)<< "\t";
          if(d.psmind2spscore(ps)!=-1){
            //Sp Score
            os<<d.psmind2spscore(psmind)<<"\t";
            //Sp Rank 
            os<<d.psmind2sp_rank(psmind)<<"\t";
          }
          //xcorr Score
      	  os<<d.psmind2xcorr(psmind)<<"\t";
      	  //xcorr rank
      	  os<<d.psmind2xcorr_rank(psmind)<<"\t";
          if(d.psmind2spscore(ps)!=-1){
            //by ions match 
            os<<d.psmind2by_ions_matched(psmind)<<"\t"; 
            //by ions total 
            os<<d.psmind2by_ions_total(psmind)<<"\t";
          }
      	  //Matches/Spectrum 
      	  os<<d.psmind2matches_spectrum(psmind)<<"\t";
      	  get_pep_seq(pep,seq,n,c);
      	  //Sequence 
      	  os<<seq<<"\t";
      	  //cleavage type
      	  os<<cleavage_type<<"\t";
      	  //protein id
      	  vector<string> prots;  
          get_protein_id(pepind,prots);
          print_protein_ids(prots,os,psmind);
      	  //Flanking_aa 
      	  os<<n<<c<<"\t";   
	  os << d.psmind2fname(psmind) << endl;
	}
    }
}


void PepRanker :: report_results_tab()
{
  ofstream of;
  ostringstream fname;
  
  fname << out_dir << "/" << fileroot << "barista.target.peptides.txt";
  of.open(fname.str().c_str());
  write_results_peptides_tab(of);
  of.close();
  fname.str("");
  
  fname << out_dir << "/" << fileroot << "barista.target.psms.txt";
  of.open(fname.str().c_str());
  write_results_psm_tab(of);
  of.close();
  fname.str("");
  
}
/******************************************************************/


void PepRanker :: report_results_xml_tab()
{
  setup_for_reporting_results();
  report_results_xml();
  report_results_tab();
  d.clear_data_pep_results(); 
  
}


void PepRanker :: print_description()
{
  cout << endl;
  cout << "\t crux q-ranker [options] <spectra> <search results>" << endl <<endl;
  cout << "REQUIRED ARGUMENTS:" << endl << endl;
  cout << "\t <spectra> Directory with ms2 files, list of ms2 files or a single ms2 file used for database search." << endl;
  cout << "\t <search results> Directory with sqt files, list of sqt files or a single sqt file with psms generated during search." << endl;
  cout << endl;
  
  cout << "OPTIONAL ARGUMENTS:" << endl << endl;
  cout << "\t [--enzyme <string>] \n \t     The enzyme used to digest the proteins in the experiment. Default trypsin." << endl;
  cout << "\t [--decoy-prefix <string>] \n \t     Specifies the prefix of the protein names that indicates a decoy. Default decoy_" << endl;
  cout << "\t [--separate-searches <string>] \n \t     If the target and decoy searches were run separately, the option then allows the user to specify the location of the decoy search results, the target database search should be provided as required argument." << endl;
  cout << "\t [--fileroot <string>] \n \t     The fileroot string will be added as a prefix to all output file names. Default = none." <<endl;
  cout << "\t [--output-dir <directory>] \n \t     The name of the directory where output files will be created. Default = crux-output." << endl;
  cout << "\t [--overwrite <T/F>] \n \t     Replace existing files (T) or exit if attempting to overwrite (F). Default=F." << endl;
  cout << "\t [--skip-cleanup <T/F>] \n \t     When set to T, prevents the deletion of lookup tables created during the preprocessing step. Default = F." << endl; 
  cout << "\t [--re-run <directory>] \n \t      Re-run PepRanker analysis using a previously computed set of lookup tables." <<endl;  
  cout << "\t [--use-spec-features <T/F>] \n \t      When set to F, use minimal feature set. Default T." <<endl;  
  cout << endl; 

}
string PepRanker:: file_extension (string str){
  string file_exten= str.substr(str.rfind('.')+1);
  return file_exten; 
}


FILE_FORMAT_T PepRanker::check_file_format(string& source) {
  string ext = file_extension(source);

    if (ext=="sqt") 
      return SQT_FORMAT;
    else if (ext=="txt") 
      return DELIMITED_FORMAT;
   return INVALID_FORMAT;
}

int PepRanker :: crux_set_command_line_options(int argc, char *argv[])
{
  const char* option_list[] = {
    "enzyme",
    "decoy-prefix",
    "separate-searches",
    "fileroot",
    "output-dir",
    "overwrite",
    "skip-cleanup",
    "re-run",
    "use-spec-features",
    "parameter-file",
    "verbosity",
     "list-of-files",
    "feature-file"
  };
  int num_options = sizeof(option_list)/sizeof(char*);

  const char* argument_list[] = {
    "spectra",
    "search results"
  };
  int num_arguments = sizeof(argument_list)/sizeof(char*);

  initialize(argument_list, num_arguments, 
	     option_list, num_options,
	     argc, argv);

  
  string sqt_source;
  string ms2_source;
  string sqt_decoy_source;
  int separate_search_flag;
  string output_directory;
  string enzyme;
  string decoy_prefix;

  string dir_with_tables;
  int found_dir_with_tables;

  bool spec_features_flag;
  bool list_of_files_flag; 
  overwrite_flag = get_boolean_parameter("overwrite");
  
  fileroot = get_string_parameter_pointer("fileroot");
  if(fileroot != "__NULL_STR")
    fileroot.append(".");
  else
    fileroot = "";

  decoy_prefix = get_string_parameter_pointer("decoy-prefix");


  enzyme = get_string_parameter_pointer("enzyme");

  spec_features_flag = get_boolean_parameter("use-spec-features");

  skip_cleanup_flag = get_boolean_parameter("skip-cleanup");
  

  dir_with_tables = get_string_parameter_pointer("re-run"); 
  if(dir_with_tables != "__NULL_STR")
    found_dir_with_tables = 1;
  else
    found_dir_with_tables = 0;

  output_directory = get_string_parameter_pointer("output-dir");

  feature_file_flag = get_boolean_parameter("feature-file");
  feature_file_name << output_directory << "/" << fileroot << "q-ranker.features.txt";

  if(found_dir_with_tables)
    {
      //set input and output dirs
      parser = new SQTParser();
      parser -> set_output_dir(dir_with_tables);
      set_input_dir(dir_with_tables);
      set_output_dir(output_directory);

      carp(CARP_INFO, "directory with tables: %s", dir_with_tables.c_str());
      carp(CARP_INFO, "output_directory: %s", output_directory.c_str());
    }else{
      ms2_source = get_string_parameter_pointer("spectra");
      sqt_source = get_string_parameter_pointer("search results");
      list_of_files_flag=get_boolean_parameter("list-of-files");
      //check file format 
     
      vector<string> files; 
      string files_list; 
      if(list_of_files_flag==1){
        ifstream f(sqt_source.c_str());
        f>>files_list; 
        while(!f.eof()){
          files.push_back(files_list);
          f>> files_list; 
        } 
      }else{
         files.push_back(sqt_source);
       }
      
      for(unsigned i=0; i<files.size();i++){
        FILE_FORMAT_T format = check_file_format(files[i]);
        switch (format) {
        case SQT_FORMAT:
          file_format_="sqt"; 
          parser = new SQTParser();
	  break;
        case DELIMITED_FORMAT:
          file_format_="txt";
	  parser = new CruxParser();
	  break;
        case INVALID_FORMAT:
          file_format_="NULL";
        default:
	  carp(CARP_FATAL, "Please enter .sqt or .txt search results"); 
        }

        parser->set_decoy_prefix(decoy_prefix);
        parser->set_enzyme(enzyme);
        if(enzyme.find("elastase") != string::npos){
           cleavage_type="elastase-full-digest";
        }else if (enzyme.find("chymotrypsin") != string::npos){
           cleavage_type="chymotrypsin-full-digest";
        }else if (enzyme.find("trypsin") != string::npos){
           cleavage_type="trypsin-full-digest";
        }else{
           cleavage_type="Null";
        }

        //num of spec features
        if(spec_features_flag)
	  parser->set_num_spec_features(3);
        else 
	  parser->set_num_spec_features(0);
        parser->write_features_header();
        
        if(parser->get_use_quadratic_features())
          parser->add_quadratic_features_header(); 
        d.get_features_header(parser->get_final_features_header());
        sqt_decoy_source = get_string_parameter_pointer("separate-searches"); 
        if(sqt_decoy_source != "__NULL_STR")
	  separate_search_flag = 1;
        else
	  separate_search_flag = 0;
      
        //set the output directory for the parser
        if(!parser->set_output_dir(output_directory))
	  return 0;
        //set input and output for the leaning algo (in and out are the same as the out for the parser)
        set_input_dir(output_directory);
        set_output_dir(output_directory);
      
        if(separate_search_flag){
	  if(!parser->set_input_sources(ms2_source, sqt_source, sqt_decoy_source)){
	    carp(CARP_FATAL, "could not extract features for training");
          }
	  parser->set_num_hits_per_spectrum(1);
	}else{
	  if(!parser->set_input_sources(ms2_source, sqt_source)){
	    carp(CARP_FATAL, "could not extract features for training");
            }
	 }
      
        //print some info
        if(format==SQT_FORMAT)
	  carp(CARP_INFO, "sqt source: %s", sqt_source.c_str());
        if(format==DELIMITED_FORMAT)
	  carp(CARP_INFO, "delmited source: %s", sqt_source.c_str());
        carp(CARP_INFO, "ms2 source: %s", ms2_source.c_str());
        carp(CARP_INFO, "output_directory: %s", output_directory.c_str());
        carp(CARP_INFO, "enzyme: %s", enzyme.c_str());
        carp(CARP_INFO, "decoy prefix: %s", decoy_prefix.c_str());
      
        if(!parser->run())
	  carp(CARP_FATAL, "Could not proceed with training.");
        parser->clear();
      }
   }
 
   if(!parser->check_input_dir(in_dir))
    carp(CARP_FATAL, "Please re-run with ms2 input and txt input.");

  return 1;

  
}



int PepRanker :: computePepNSAF()
{

  //create auxiliary psmind_to_ind index
  vector<int>psmind_to_ind;
  psmind_to_ind.resize(psmtrainset.size(),0);
  for(int i = 0; i < psmtrainset.size(); i++)
    psmind_to_ind[psmtrainset[i].psmind] = i;

  //compute NSAF for each peptide
  double sum = 0.0;
  for(int k = 0; k < trainset.size(); k++)
    {
      if(trainset[k].label == 1)
	{
	  int pepind = trainset[k].pepind;
	  string pep = d.ind2pep(pepind);
	  string seq, n,c;
	  get_pep_seq(pep,seq,n,c);
	  int len = seq.size();
	  
	  int cnt = 0;
	  int num_psms = d.pepind2num_psm(pepind);
	  int *psminds = d.pepind2psminds(pepind);
	  for (int j = 0; j < num_psms; j++)
	    {
	      int psmind = psminds[j];
	      int ind  = psmind_to_ind[psmind];
	      assert(psmtrainset[ind].psmind == psmind);
	      if(psmtrainset[ind].q <= 0.01)
		cnt++;
	    }
	
	  //compute NSAF for this protein
	  double nsaf = 0.0;
	  if(len != 0)
	    nsaf = (double)cnt/(double)len;
	  trainset[k].nsaf = nsaf;
	  sum += nsaf;
 	}
    }
  
  if(sum > 0)
    {
      for(int i = 0; i < trainset.size(); i++)
	{
	  if(trainset[i]. label == 1)
	    trainset[i].nsaf /= sum;
	}
    }
  else
    return 0;

  /*
  double check = 0.0;
  for(int i = 0; i < trainset.size(); i++)
    {
      if(trainset[i]. label == 1)
	check += trainset[i].nsaf;
    }
  cout << check << endl;
  */
  psmind_to_ind.clear();
  return 1;
}


void PepRanker :: computePEP(){
  carp(CARP_DEBUG, "Computing PEPs");
  vector<double> target_scores_vect;
  vector<double> decoy_scores_vect;

  // pull out the target and decoy scores
  for(int i = 0; i < psmtrainset.size(); i++){
    if( psmtrainset[i].label == 1 ){
      target_scores_vect.push_back(psmtrainset[i].score);
    } else { // == -1
      decoy_scores_vect.push_back(psmtrainset[i].score);
    }
  }

  int num_targets = target_scores_vect.size();
  int num_decoys = decoy_scores_vect.size();
  carp(CARP_DEBUG, "Found %d targets and %d decoys", num_targets, num_decoys); 

  // copy them to an array as required by the compute_PEP method
  double* target_scores = new double[num_targets];
  copy(target_scores_vect.begin(), target_scores_vect.end(), target_scores);
  double* decoy_scores = new double[num_decoys];
  copy(decoy_scores_vect.begin(), decoy_scores_vect.end(), decoy_scores);

  double* PEPs = compute_PEP(target_scores, num_targets, 
                             decoy_scores, num_decoys);

  // fill in the data set with the new scores for the targets
  int target_idx = 0;
  for(int full_idx = 0; full_idx < psmtrainset.size(); full_idx++){
    if( psmtrainset[full_idx].label == 1 ){
      psmtrainset[full_idx].PEP = PEPs[target_idx];
      target_idx++; 
    } // else, skip decoys
  }

  delete target_scores;
  delete decoy_scores;
  delete PEPs;

  /** 
   * Calculate Peptide PEPs
   */
  target_scores_vect.clear();
  decoy_scores_vect.clear();
  
  //pull out the target and decoy scores for peptide PEPs
  for(int i=0;i<trainset.size();i++){
    if(trainset[i].label==1)
      target_scores_vect.push_back(trainset[i].score);
    else 
      decoy_scores_vect.push_back(trainset[i].score);
  }
  num_targets = target_scores_vect.size();
  num_decoys = decoy_scores_vect.size();
  carp(CARP_DEBUG, "Found %d targets and %d decoys", num_targets, num_decoys);  

  // copy them to an array as required by the compute_PEP method
  target_scores = new double[num_targets];
  copy(target_scores_vect.begin(), target_scores_vect.end(), target_scores);
  decoy_scores = new double[num_decoys];
  copy(decoy_scores_vect.begin(), decoy_scores_vect.end(), decoy_scores);
  PEPs = compute_PEP(target_scores, num_targets, 
                             decoy_scores, num_decoys);

  // fill in the data set with the new scores for the targets
  target_idx = 0;
  for(int full_idx = 0; full_idx < trainset.size(); full_idx++){
    if( trainset[full_idx].label == 1 ){
      trainset[full_idx].PEP = PEPs[target_idx];
      target_idx++; 
    } // else, skip decoys
  }

  delete target_scores;
  delete decoy_scores;
  delete PEPs;

}


void PepRanker :: write_results_pep_xml(PepXMLWriter& xmlfile)
{
  xmlfile.writeHeader();

  bool* scores_to_print = new bool[NUMBER_SCORER_TYPES];
  for(int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
    scores_to_print[score_idx] = false;
  }
  scores_to_print[SP] = true; 
  scores_to_print[XCORR] = true;
  scores_to_print[BARISTA_SCORE] = true;
  scores_to_print[BARISTA_QVALUE] = true;
  scores_to_print[BARISTA_PEP] = true;

  xmlfile.SetScoresComputed(scores_to_print);

  double* scores = new double[NUMBER_SCORER_TYPES];

  for(int i = 0; i < psmtrainset.size(); i++) {
    // only print target psms
    if( psmtrainset[i].label == -1 ){
      continue;
    }

      int psmind = psmtrainset[i].psmind;

      // spectrum info
      int scan = d.psmind2scan(psmind);
      const char* filename = d.psmind2fname(psmind).c_str();
      char** path_name = parse_filename_path_extension(filename, NULL);
      filename = path_name[0];
      double spectrum_mass = d.psmind2precursor_mass(psmind); 
      int charge = d.psmind2charge(psmind);

      // peptide info
      int pepind = d.psmind2pepind(psmind);
      string pep = d.ind2pep(pepind);
      string modified_sequence, n,c;
      get_pep_seq(pep, modified_sequence, n, c);
      char* sequence = unmodify_sequence(modified_sequence.c_str());
      string flanking_aas = n + c;
      double peptide_mass = d.psmind2peptide_mass(psmind);

      // protein info
      int num_proteins = d.pepind2num_prot(pepind); 
      int* protein_indexes = d.pepind2protinds(pepind);
      vector<string> protein_names;
      vector<string> protein_descriptions;
      for(int prot_idx = 0; prot_idx < num_proteins; prot_idx++){
        string protein_name = d.ind2prot(protein_indexes[prot_idx]);
	protein_names.push_back(protein_name);
        protein_descriptions.push_back("");
        flanking_aas += "," + n + c; // todo
      }

      // psm info
      int* psm_rank=new int[NUMBER_SCORER_TYPES] ;
      psm_rank[XCORR]=1;
      double delta_cn = d.psmind2deltaCn(psmind);
      scores[XCORR] = d.psmind2xcorr(psmind);
      scores[SP] = d.psmind2spscore(psmind);
      scores[QRANKER_SCORE] = psmtrainset[i].score;
      scores[QRANKER_QVALUE] = psmtrainset[i].q;
      scores[QRANKER_PEP] = psmtrainset[i].PEP;
        
      xmlfile.writePSM(scan, filename, spectrum_mass, charge, psm_rank,
                       sequence, modified_sequence.c_str(),
                       peptide_mass, num_proteins,
                       flanking_aas.c_str(), protein_names, 
                       protein_descriptions, delta_cn, scores_to_print, scores,d.psmind2by_ions_matched(psmind),
                       d.psmind2by_ions_total(psmind), d.psmind2matches_spectrum(psmind));

      free(sequence);
      if( path_name[0] ){
        free(path_name[0]);
      }
      if( path_name[1] ){
        free(path_name[1]);
      }
      free(path_name);
  }

  xmlfile.writeFooter();
}



int PepRanker::main(int argc, char **argv) {

  
  if(!crux_set_command_line_options(argc, argv))
    return 1;
  
  run();
  if(skip_cleanup_flag != 1)
    parser->clean_up(out_dir);
    
  return 0;
}   

bool PepRanker :: needsOutputDirectory()
{
  return true;
}


string PepRanker::getName() {
  return "q-ranker";
}
string PepRanker::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the new q-ranker algorithm.";

}

COMMAND_T PepRanker::getCommand(){
  return QRANKER_COMMAND;
}


void PepRanker :: get_protein_id(int pepind, vector<string> &prot){
  int * protinds= d.pepind2protinds(pepind) ;
  unsigned prot_length=d.pepind2num_prot(pepind);  
  for (unsigned k=0;k<prot_length;k++){
    string protein_str= d.ind2prot( protinds[k]);
    prot.push_back(protein_str); 
  }
}


void PepRanker :: print_protein_ids(vector<string> &prots, ofstream &os, int psmind){
    
 // TODO: find peptide position in SQT files 
 //sqt files do not return peptide position in the protein.
 //if the search result file is txt we can find file peptide_pos and print it 
 // in front of protein, else do not print anything

  if(file_format_=="txt"){
    for (unsigned int j=0;j<prots.size();j++){
      os << prots[j];
      if (d.psmind2peptide_position(psmind) > -1) {
	os << "("<<d.psmind2peptide_position(psmind)<<")";
      }
      if (j == prots.size() - 1) {
	os << "\t";
      } else {
	os << ",";
      }
    }
  }else if(file_format_=="sqt"){
    for (unsigned int j=0;j<prots.size();j++){
      if(j==prots.size()-1)
      	os<<prots[j]<<"\t";
      else
        os<<prots[j]<<",";
    }
  }
  prots.clear();
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
