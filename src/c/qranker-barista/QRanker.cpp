#include "QRanker.h"
#include "modifications.h"

QRanker::QRanker() :  
  seed(0),
  selectionfdr(0.01),
  num_hu(4),
  mu(0.01),
  weightDecay(0.000),
  in_dir(""), 
  out_dir(""), 
  skip_cleanup_flag(0),
  overwrite_flag(0),
  fileroot(""),
  feature_file_flag(0)
{
}

QRanker::~QRanker()
{
  delete [] max_net_gen;
  delete [] max_net_targ;
  delete [] nets;
}

int QRanker :: getOverFDR(PSMScores &set, NeuralNet &n, double fdr)
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


void QRanker :: getMultiFDR(PSMScores &set, NeuralNet &n, vector<double> &qvalues)
{
  double *r;
  double* featVec;
   
  for(int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = n.fprop(featVec);
      set[i].score = r[0];
    }
 
  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti);
}

void QRanker :: getMultiFDRXCorr(PSMScores &set, vector<double> &qvalues)
{
  double r;
  double* featVec;
   
  for(int i = 0; i < set.size(); i++)
    {
      featVec = d.psmind2features(set[i].psmind);
      r = featVec[3];
      set[i].score = r;
    }
 
  for(unsigned int ct = 0; ct < qvalues.size(); ct++)
    overFDRmulti[ct] = 0;
  set.calcMultiOverFDR(qvalues, overFDRmulti);
}


void QRanker :: printNetResults(vector<int> &scores)
{
  double qv; int fdr;
  cerr << "QVALS SCORES:: ";
  for(unsigned int count = 0; count < qvals.size();count++)
    {  
      qv=qvals[count]; 
      fdr = scores[count];
      cerr << qv << ":" << fdr << " ";
    }
  cerr << endl;
}




void QRanker :: write_results()
{
  //write out the results of the general net
  trainset.clear();
  testset.clear();
  thresholdset.clear();
  d.load_labels_psm_training();
  PSMScores::fillFeaturesFull(fullset, d);
  d.clear_labels_psm_training();

  //choose the best net for the selectionfdr
  int max_fdr = 0;
  int fdr = 0;
  int ind = 0;
  for(unsigned int count = 0; count < qvals.size();count++)
    {
      fdr = getOverFDR(fullset, max_net_targ[count], selectionfdr);
      if(fdr > max_fdr)
	{
	  max_fdr = fdr;
	  ind = count;
	}
    }
  net = max_net_targ[ind];
  getOverFDR(fullset,net, selectionfdr);
  d.clear_data_psm_training();
  
  ostringstream fname;
  fname << out_dir << "/" << fileroot << "q-ranker.target.psms.txt";

  d.load_data_psm_results();
  computePEP();

  ofstream f1(fname.str().c_str()); 
  write_results_psm_tab(f1);
  f1.close();
  fname.str("");

  fname << out_dir << "/" << fileroot << "q-ranker_output.xml";
  PepXMLWriter xmlfile;
  xmlfile.openFile(fname.str().c_str(), overwrite_flag);
  write_results_psm_xml(xmlfile);
  xmlfile.closeFile();
  d.clear_data_psm_results();
}

void QRanker :: write_results_psm_tab(ofstream &os)
{
  os << "scan" << "\t" << "charge" << "\t" << "q-ranker q-value" << "\t" 
     << "q-ranker score" << "\t" << "PEP\t" << "sequence" << "\t" 
     << "filename" << endl;

  for(int i = 0; i < fullset.size(); i++)
    {
      if( fullset[i].label == 1 ){ // only print target psms
        int psmind = fullset[i].psmind;
        int pepind = d.psmind2pepind(psmind);
        os << d.psmind2scan(psmind) << "\t" << d.psmind2charge(psmind) << "\t" 
           << fullset[i].q << "\t" << fullset[i].score << "\t"  
           << fullset[i].PEP << "\t"
           << d.ind2pep(pepind) << "\t" << d.psmind2fname(psmind) << endl;
      }
    }
}


void QRanker :: get_pep_seq(string &pep, string &seq, string &n, string &c)
{
  string tmp;
  int pos;
  pos = pep.find(".");
  n = pep.at(pos-1); 
  tmp = pep.substr(pos+1, pep.size());

  pos = tmp.find(".");
  c = tmp.at(pos+1);
  seq = tmp.substr(0, pos);
}


void QRanker ::write_results_psm_xml(PepXMLWriter& xmlfile)
{
  xmlfile.writeHeader();

  bool* scores_to_print = new bool[NUMBER_SCORER_TYPES];
  for(int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
    scores_to_print[score_idx] = false;
  }
  scores_to_print[SP] = true; 
  scores_to_print[XCORR] = true;
  scores_to_print[QRANKER_SCORE] = true;
  scores_to_print[QRANKER_QVALUE] = true;
  scores_to_print[QRANKER_PEP] = true;

  xmlfile.SetScoresComputed(scores_to_print);

  double* scores = new double[NUMBER_SCORER_TYPES];

  for(int i = 0; i < fullset.size(); i++)
    {
      // only print target psms
      if( fullset[i].label == -1 ){
        continue;
      }

      int psmind = fullset[i].psmind;

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
      int psm_rank = 1;
      double delta_cn = d.psmind2deltaCn(psmind);
      scores[XCORR] = d.psmind2xcorr(psmind);
      scores[SP] = d.psmind2spscore(psmind);
      scores[QRANKER_SCORE] = fullset[i].score;
      scores[QRANKER_QVALUE] = fullset[i].q;
      scores[QRANKER_PEP] = fullset[i].PEP;

      xmlfile.writePSM(scan, filename, spectrum_mass, charge, psm_rank,
                       sequence, modified_sequence.c_str(),
                       peptide_mass, num_proteins,
                       flanking_aas.c_str(), protein_names, 
                       protein_descriptions, delta_cn, scores_to_print, scores);

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



void QRanker :: write_max_nets(string filename, NeuralNet* max_net)
{
  //write out the results of the general net

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




void QRanker :: write_unique_peptides(string filename, NeuralNet* max_net)
{
  //write out the results of the general net
  ostringstream s1;
  s1 << filename << ".txt";
  ofstream f1(s1.str().c_str());
  
  for(int count = 0; count < num_qvals; count++)
    {
      net = max_net[count];
      int r = getOverFDR(testset,net, qvals[count]);
      set<int> peps;
      int cn = 0;
      for(int i = 0; i < testset.size();i++)
	{
	  if(testset[i].label == 1)
	    {
	      cn++;
	      int pepind = d.psmind2pepind(testset[i].psmind);
	      peps.insert(pepind);
	    }
	  if(cn > r) break;
	}
      int num_tst = 0;
      for(set<int>::iterator it = peps.begin(); it != peps.end(); it++)
	num_tst++;
      peps.clear();

      int r1 = getOverFDR(trainset,net, qvals[count]);
      cn = 0;
      for(int i = 0; i < trainset.size();i++)
	{
	  if(trainset[i].label == 1)
	    {
	      cn++;
	      int pepind = d.psmind2pepind(trainset[i].psmind);
	      peps.insert(pepind);
	    }
	  if(cn > r1) break;
	}
      int num_trn = 0;
      for(set<int>::iterator it = peps.begin(); it != peps.end(); it++)
	num_trn++;
      peps.clear();
      f1 << qvals[count] << " " << num_trn << " " << num_tst << "\n";
    }

  f1.close();
}


void QRanker :: write_num_psm_per_spectrum(NeuralNet* max_net)
{
  //write out the results of the general net
  //ostringstream s1;
  //s1 << filename << ".txt";
  //ofstream f1(s1.str().c_str());
  
  int count = 5; 
  net = max_net[count];
  int r = getOverFDR(trainset,net, qvals[count]);
      
  map<int,set<int> > scan_to_pepinds;
  int cn = 0;
  for(int i = 0; i < trainset.size();i++)
    {
      if(trainset[i].label == 1)
	{
	  cn++;
	  int pepind = d.psmind2pepind(trainset[i].psmind);
	  int scan = d.psmind2scan(trainset[i].psmind);
	  if(scan_to_pepinds.find(scan) == scan_to_pepinds.end())
	    {
	      set<int> peps;
	      scan_to_pepinds[scan] = peps;
	    }
	  (scan_to_pepinds[scan]).insert(pepind);
	  
	}
      if(cn > r) break;
    }
    
  vector<int> counts;
  counts.resize(11,0);

  for(map<int,set<int> >::iterator it = scan_to_pepinds.begin();
      it != scan_to_pepinds.end(); it++)
    {
      int cnt =  (it->second).size();
      counts[cnt]++;
    }
  for(unsigned int i = 0; i < counts.size();i++)
    cout << i << " " << counts[i] << endl;

  //f1.close();
}



/*********************** training net functions*******************************************/

void QRanker :: train_net_sigmoid(PSMScores &set, int interval)
{
  double *r;
  int label;
  double *gc = new double[1];
  for(int i = 0; i < set.size(); i++)
    { 
      unsigned int ind;
      ind = rand()%(interval);
      //pass both through the net
      r = net.fprop(d.psmind2features(set[ind].psmind));
      label = set[ind].label;
      double a = exp(label*r[0]);
      net.clear_gradients();
      gc[0] = -a/((1+a)*(1+a))*label;
      net.bprop(gc);
      net.update(mu,weightDecay);

    }
  delete[] gc;
}

void QRanker :: train_net_hinge(PSMScores &set, int interval)
{
  double *r;
  int label;
  double *gc = new double[1];
  for(int i = 0; i < set.size(); i++)
    { 
      unsigned int ind;
      ind = rand()%(interval);
      //pass both through the net
      r = net.fprop(d.psmind2features(set[ind].psmind));
      label = set[ind].label;
      if(label*r[0]<1)
	{
	  net.clear_gradients();
	  gc[0] = -1.0*label;
	  net.bprop(gc);
	  net.update(mu,weightDecay);
	}
    }
  delete[] gc;
}


void QRanker :: train_net_ranking(PSMScores &set, int interval)
{
  double *r1;
  double *r2;
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
	ind1 = rand()%(interval);
      if(ind1>set.size()-1) continue;
      if(set[ind1].label == 1)
	label_flag = -1;
      else
	label_flag = 1;
      
      int cn = 0;
      while(1)
	{
	  ind2 = rand()%(interval);
	  if(ind2>set.size()-1) continue;
	  if(set[ind2].label == label_flag) break;
	  if(cn > 1000)
	    {
	      ind2 = rand()%set.size();
	      break;
	    }
	  cn++;
	}
      
      //pass both through the net
      r1 = nets[0].fprop(d.psmind2features(set[ind1].psmind));
      r2 = nets[1].fprop(d.psmind2features(set[ind2].psmind));
      diff = r1[0]-r2[0];
      

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
	      gc[0] = -1.0*label;
	      nets[0].bprop(gc);
	      gc[0] = 1.0*label;
	      nets[1].bprop(gc);
	      net.update(mu,weightDecay);
	    }
	  
	}
    }
  delete[] gc;
  
}


void QRanker :: count_pairs(PSMScores &set, int interval)
{
  double *r1;
  double *r2;
  double diff = 0;
  int label = 1;
  double err = 0.0;

  for(int i = 0; i < set.size()*100; i++)
    { 
      int ind1, ind2;
      int label_flag = 1;
      //get the first index
      if(interval == 0)
	ind1 = 0;
      else
	ind1 = rand()%(interval);
      if(ind1>set.size()-1) continue;
      if(set[ind1].label == 1)
	label_flag = -1;
      else
	label_flag = 1;
      
      int cn = 0;
      while(1)
	{
	  ind2 = rand()%(interval);
	  if(ind2>set.size()-1) continue;
	  if(set[ind2].label == label_flag) break;
	  if(cn > 1000)
	    {
	      ind2 = rand()%set.size();
	      break;
	    }
	  cn++;
	}
      
      //pass both through the net
      r1 = nets[0].fprop(d.psmind2features(set[ind1].psmind));
      r2 = nets[1].fprop(d.psmind2features(set[ind2].psmind));
      diff = r1[0]-r2[0];
      

      label=0;
      if(  set[ind1].label==1 && set[ind2].label==-1)
	label=1;
      if( set[ind1].label==-1 && set[ind2].label==1)
	    label=-1;
      
      if(label != 0)
	{
	  if(label*diff<1)
	    {
	      err += 1-label*diff;
	      
	    }
	  
	}
    }

  cout << "er" << err/trainset.size() << endl;
}



void QRanker :: train_many_general_nets()
{
  interval = trainset.size();
  for(int i=0;i<switch_iter;i++) {
    //train_net_hinge(trainset, interval);
    train_net_ranking(trainset, interval);
    //train_net_sigmoid(trainset, interval);

       
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
	
	/*
	cout << "Iteration " << i << " : \n";
	cout << "trainset: ";
	getMultiFDR(trainset,net,qvals);
	printNetResults(overFDRmulti);
	cout << "\n";
	cout << "testset: ";
	getMultiFDR(testset,net,qvals);
	printNetResults(overFDRmulti);
	cout << "\n";
	*/
      }
  }

}

void QRanker :: train_many_target_nets()
{

  int  thr_count = num_qvals-1;
  while (thr_count > 0)
    {
      net.copy(max_net_gen[thr_count]);
        
      carp(CARP_INFO, "training threshold %d", thr_count);
      //cout << "training thresh " << thr_count  << "\n";
      //interval = getOverFDR(trainset, net, qvals[thr_count]);
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
	    
	    /*
	    cout << "Iteration " << i << " : \n";
	    getMultiFDR(testset,net,qvals);
	    cout << "testset: ";
	    printNetResults(overFDRmulti);
	    cout << "\n";
	    */
	  }

      }
      thr_count -= 3;
    }
}



void QRanker::train_many_nets()
{
  switch_iter =30;
  niter = 40;
   
  num_qvals = 14;
  qvals.resize(num_qvals,0.0);
  qvals1.resize(num_qvals,0.0);
  qvals2.resize(num_qvals,0.0);

  overFDRmulti.resize(num_qvals,0);
  ave_overFDR.resize(num_qvals,0);
  max_overFDR.resize(num_qvals,0); 
  max_net_gen = new NeuralNet[num_qvals];
  max_net_targ = new NeuralNet[num_qvals];
  
  double q = 0.0;
  for(int count = 0; count < num_qvals; count++)
    {
      qvals[count] = q;
      if (count < 2)
	qvals1[count] = q;
      else
	qvals1[count] = q-0.005;
      qvals2[count] = q+0.005;

      if(q < 0.01)
	q+=0.0025;
      else
	q+=0.01;
    }
  
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

  nets = new NeuralNet[2];
  nets[0].clone(net);
  nets[1].clone(net);

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
  
  /*
  cout << "Before iterating\n";
  cout << "trainset: ";
  getMultiFDRXCorr(trainset,qvals);
  printNetResults(overFDRmulti);
  cout << "\n";
  for(int count = 0; count < num_qvals;count++)
    if(overFDRmulti[count] > max_overFDR[count])
      max_overFDR[count] = overFDRmulti[count];
  cout << "testset: ";
  getMultiFDR(testset,net,qvals);
  printNetResults(overFDRmulti);
  cout << "\n";
  */

  train_many_general_nets();
  
  //copy the general net into target nets;
  for(int count = 0; count < num_qvals; count++)
    max_net_targ[count] = max_net_gen[count];

  train_many_target_nets();  
 
  

  ostringstream fname;
  fname << out_dir << "/" << fileroot << "q-ranker.psms.at.fdr.thresholds.txt";;
  write_max_nets(fname.str(), max_net_targ);

  
}

int QRanker::run( ) {
  carp(CARP_INFO, "reading data");
  
  ostringstream res;
  res << out_dir << "/qranker_output";
  res_prefix = res.str();
    
  d.load_data_psm_training();
  if(feature_file_flag)
    {
      string str = feature_file_name.str();
      if(!d.print_features(str))
	carp(CARP_INFO, "could not open file %s for writing features. Feature file will not be written", feature_file_name.str().c_str());
    }
  d.normalize_psms();
  d.load_labels_psm_training();
  PSMScores::fillFeaturesSplit(trainset, testset, d, 0.75);
  thresholdset = trainset;
  d.clear_labels_psm_training();
  train_many_nets();
  
  write_results();

  return 0;
}


/*
int QRanker::set_command_line_options(int argc, char **argv)
{
  vector<string> fnames;
  int arg = 1;
  while(arg < argc)
    {
      fnames.push_back(argv[arg]);
      arg++;
    }
  string out_dir = "crux-output";
  pars.set_output_dir(out_dir);
  if(!pars.run_on_xlink(fnames))
    return 0;
  pars.clear();
  return 1;
}
*/

void QRanker :: print_description()
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
  cout << "\t [--re-run <directory>] \n \t      Re-run QRanker analysis using a previously computed set of lookup tables." <<endl;  
  cout << "\t [--use-spec-features <T/F>] \n \t      When set to F, use minimal feature set. Default T." <<endl;  
  cout << endl; 

}

int QRanker :: crux_set_command_line_options(int argc, char *argv[])
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

  overwrite_flag = get_boolean_parameter("overwrite");

  fileroot = get_string_parameter_pointer("fileroot");
  if(fileroot != "__NULL_STR")
    fileroot.append(".");
  else
    fileroot = "";

  decoy_prefix = get_string_parameter_pointer("decoy-prefix");
  sqtp.set_decoy_prefix(decoy_prefix);

  enzyme = get_string_parameter_pointer("enzyme");
  sqtp.set_enzyme(enzyme);

  spec_features_flag = get_boolean_parameter("use-spec-features");
  //num of spec features
  if(spec_features_flag)
    sqtp.set_num_spec_features(3);
  else
    sqtp.set_num_spec_features(0);

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
      sqtp.set_output_dir(dir_with_tables);
      set_input_dir(dir_with_tables);
      set_output_dir(output_directory);

      carp(CARP_INFO, "directory with tables: %s", dir_with_tables.c_str());
      carp(CARP_INFO, "output_directory: %s", output_directory.c_str());
    }
  else
    {
      ms2_source = get_string_parameter_pointer("spectra");
      sqt_source = get_string_parameter_pointer("search results");
      sqt_decoy_source = get_string_parameter_pointer("separate-searches"); 
      if(sqt_decoy_source != "__NULL_STR")
	separate_search_flag = 1;
      else
	separate_search_flag = 0;
      
      //set the output directory for the parser
      if(!sqtp.set_output_dir(output_directory))
	return 0;
      //set input and output for the leaning algo (in and out are the same as the out for the parser)
      set_input_dir(output_directory);
      set_output_dir(output_directory);
      
      if(separate_search_flag)
	{
	  if(!sqtp.set_input_sources(ms2_source, sqt_source, sqt_decoy_source))
	    carp(CARP_FATAL, "could not extract features for training");
	  sqtp.set_num_hits_per_spectrum(1);
	}
      else
	{
	  if(!sqtp.set_input_sources(ms2_source, sqt_source))
	    carp(CARP_FATAL, "could not extract features for training");
	}
      
      //print some info
      carp(CARP_INFO, "sqt source: %s", sqt_source.c_str()); 
      carp(CARP_INFO, "ms2 source: %s", ms2_source.c_str());
      carp(CARP_INFO, "output_directory: %s", output_directory.c_str());
      carp(CARP_INFO, "enzyme: %s", enzyme.c_str());
      carp(CARP_INFO, "decoy prefix: %s", decoy_prefix.c_str());
      
      if(!sqtp.run())
	carp(CARP_FATAL, "Could not proceed with training.");
      sqtp.clear();
    }
 
   if(!sqtp.check_input_dir(in_dir))
    carp(CARP_FATAL, "Please re-run with ms2 input and sqt input.");

  return 1;

  
}

/**
 * Uses the target and decoy scores to compute posterior error
 * probabilities. 
 */
void QRanker::computePEP(){
  carp(CARP_DEBUG, "Computing PEPs");
  vector<double> target_scores_vect;
  vector<double> decoy_scores_vect;

  // pull out the target and decoy scores
  for(int i = 0; i < fullset.size(); i++){
    if( fullset[i].label == 1 ){
      target_scores_vect.push_back(fullset[i].score);
    } else { // == -1
      decoy_scores_vect.push_back(fullset[i].score);
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
  for(int full_idx = 0; full_idx < fullset.size(); full_idx++){
    if( fullset[full_idx].label == 1 ){
      fullset[full_idx].PEP = PEPs[target_idx];
      target_idx++; 
    } // else, skip decoys
  }

  delete target_scores;
  delete decoy_scores;
  delete PEPs;
}



/*
int QRanker :: set_command_line_options(int argc, char *argv[])
{
  string db_source;
  string sqt_source;
  string sqt_decoy_source;
  int separate_search_flag = 0;
  string ms2_source;
  string output_directory = "crux-output";
  string enzyme = "trypsin";
  string decoy_prefix = "reverse_";
  string dir_with_tables = "";
  int found_dir_with_tables = 0;
  int spec_features_flag = 1;
  int arg = 1;
  int database_exists = 0;
  
  while (arg < argc)
    {
      string str = argv[arg];
      //parse the options
      if(str.find("--") != string::npos)
	{
	  //found enzyme
	  if(str.find("enzyme") != string::npos)
	    {
	      arg++;
	      enzyme = argv[arg];
	      sqtp.set_enzyme(enzyme);
	    }
	  //found decoy prefix
	  else if(str.find("decoy-prefix") != string::npos)
	    {
	      arg++;
	      decoy_prefix = argv[arg];
	      sqtp.set_decoy_prefix(decoy_prefix);
	    }
	  //found output directory
	  else if(str.find("output-dir") != string::npos)
	    {
	      arg++;
	      output_directory = argv[arg];
	      set_output_dir(output_directory);
	    }
	  //found overwrite directory
	  else if(str.find("overwrite") != string::npos)
	    {
	      arg++;
	      string opt = argv[arg];
	      if (opt.compare("T") == 0)
		overwrite_flag = 1;
	      else
		overwrite_flag = 0;
	    }
	  //found fileroot
	  else if(str.find("fileroot") != string::npos)
	    {
	      arg++;
	      fileroot = argv[arg];
	      fileroot.append(".");
	    }
	  //no cleanup
	  else if(str.find("skip-cleanup") != string::npos)
	    {
	      arg++;
	      string opt = argv[arg];
	      if (opt.compare("T") == 0)
		{
		  skip_cleanup_flag = 1;
		  cout << "INFO: will not do cleanup" << endl;
		}
	    }
	  //
	  else if(str.find("re-run") != string::npos)
	    {
	      arg++;
	      dir_with_tables = argv[arg];
	      found_dir_with_tables = 1;
	      cout << "INFO: directory with preprocessed data: " << dir_with_tables << endl;
	    }
	  //found overwrite directory
	  else if(str.find("spec-features") != string::npos)
	    {
	      arg++;
	      string opt = argv[arg];
	      if (opt.compare("T") == 0)
		spec_features_flag = 1;
	      else
		spec_features_flag = 0;
	    }
	  //found separate search
	  else if(str.find("separate-search") != string::npos)
	    {
	      arg++;
	      string opt = argv[arg];
	      sqt_decoy_source = opt;
	      separate_search_flag = 1;
	      cout << sqt_decoy_source << endl;
	    }
	  else
	    {
	      cout << "FATAL: option " << str << " does not exist" << endl;
	      return 0;
	    }
	}
      else
	break;
      arg++;
    }
  
  ostringstream cmd;
  for(int k = 0; k < argc; k++)
    cmd << argv[k] << " ";

  if(found_dir_with_tables)
    {
      sqtp.set_output_dir(dir_with_tables, overwrite_flag);
      set_input_dir(dir_with_tables);
      set_output_dir(output_directory);

      set_verbosity_level(CARP_INFO);
      initialize_parameters();
      //open log file
      set_boolean_parameter("overwrite", overwrite_flag);
      set_string_parameter("output-dir", output_directory.c_str());
      ostringstream logfname;
      logfname << fileroot << "qranker.log.txt";
      string str = logfname.str();
      char *log_file = my_copy_string(str.c_str());
      open_log_file(&log_file);
      free(log_file);

      carp(CARP_INFO, "COMMAND: %s", cmd.str().c_str());
      carp(CARP_INFO, "directory with tables: %s", dir_with_tables.c_str());
      carp(CARP_INFO, "output_directory: %s", output_directory.c_str());
      carp(CARP_INFO, "enzyme: %s", enzyme.c_str());
      carp(CARP_INFO, "decoy prefix: %s", decoy_prefix.c_str());
      if(fileroot.compare("") != 0)
	carp(CARP_INFO, "fileroot: %s", fileroot.c_str());

      //write out a parameter file
      stringstream fname;
      fname << output_directory << "/" << fileroot << "qranker.params.txt";
      ofstream fparam(fname.str().c_str());
      fparam << "enzyme=" << enzyme << endl;
      fparam << "decoy prefix=" << decoy_prefix << endl;
      if(separate_search_flag)
	fparam << "separate search=" << sqt_decoy_source << endl;
      fparam << "fileroot=" << fileroot << endl;
      fparam << "output directory=" << output_directory << endl;
      if(skip_cleanup_flag)
	fparam << "skip-cleanup=T" << endl;
      else
	fparam << "skip-cleanup=F" << endl;
      fparam << "re-run=" << dir_with_tables << endl;
      if(spec_features_flag)
	fparam << "use spec features=T" << endl;
      else
	fparam << "use spec features=F" << endl;
      fparam.close();
    }
  else
    {
      //if(argc-arg < 3)
      if(argc-arg < 2)
	{
	  print_description();
	  return 0;
	} 
      if(argc-arg == 2)
	{
	  //db_source = argv[arg]; arg++; 
	  ms2_source = argv[arg]; arg++;
	  sqt_source = argv[arg];
	}
      if(argc-arg == 3)
	{
	  db_source = argv[arg]; arg++; 
	  ms2_source = argv[arg]; arg++;
	  sqt_source = argv[arg];
	  database_exists = 1;
	}
      //set the output directory
      if(!sqtp.set_output_dir(output_directory, overwrite_flag))
      	return 0;
      set_input_dir(output_directory);
      set_output_dir(output_directory);

      set_verbosity_level(CARP_INFO);
      initialize_parameters();
      //open log file
      set_boolean_parameter("overwrite", overwrite_flag);
      set_string_parameter("output-dir", output_directory.c_str());
      ostringstream logfname;
      logfname << fileroot << "qranker.log.txt";
      string str = logfname.str();
      char *log_file = my_copy_string(str.c_str());
      open_log_file(&log_file);
      free(log_file);
      
      //write out a parameter file
      stringstream fname;
      fname << output_directory << "/" << fileroot << "qranker.params.txt";
      ofstream fparam(fname.str().c_str());
      fparam << "enzyme=" << enzyme << endl;
      fparam << "decoy prefix=" << decoy_prefix << endl;
      if(separate_search_flag)
	fparam << "separate search=" << sqt_decoy_source << endl;
      fparam << "fileroot=" << fileroot << endl;
      fparam << "output directory=" << output_directory << endl;
      if(skip_cleanup_flag)
	fparam << "skip-cleanup=T" << endl;
      else
	fparam << "skip-cleanup=F" << endl; 
      if(spec_features_flag)
	fparam << "use spec features=T" << endl;
      else
	fparam << "use spec features=F" << endl;
      fparam.close();

      if(database_exists)
	{
	  if(!sqtp.set_database_source(db_source))
	    carp(CARP_FATAL, "could not extract features for training");
	}
      if(separate_search_flag)
	{
	  if(!sqtp.set_input_sources(ms2_source, sqt_source, sqt_decoy_source))
	    carp(CARP_FATAL, "could not extract features for training");
	  sqtp.set_num_hits_per_spectrum(1);
	}
      else
	{
	  if(!sqtp.set_input_sources(ms2_source, sqt_source))
	    carp(CARP_FATAL, "could not extract features for training");
	}

       carp(CARP_INFO, "COMMAND: %s", cmd.str().c_str());
       carp(CARP_INFO, "database source: %s", db_source.c_str());
       carp(CARP_INFO, "sqt source: %s", sqt_source.c_str()); 
       carp(CARP_INFO, "ms2 source: %s", ms2_source.c_str());
       carp(CARP_INFO, "output_directory: %s", output_directory.c_str());
       carp(CARP_INFO, "enzyme: %s", enzyme.c_str());
       carp(CARP_INFO, "decoy prefix: %s", decoy_prefix.c_str());
       if(fileroot.compare("") != 0)
	 carp(CARP_INFO, "fileroot: %s", fileroot.c_str());
       
       //num of spec features
       if(spec_features_flag)
	 sqtp.set_num_spec_features(7);
       else
	 sqtp.set_num_spec_features(0);
       if(!sqtp.run())
	 carp(CARP_FATAL, "Could not proceed with training.");
       sqtp.clear();
    }

  if(!sqtp.check_input_dir(in_dir))
    carp(CARP_FATAL, "Please re-run with database, ms2 input and sqt input.");

  return 1;
  

}
*/



int QRanker::main(int argc, char **argv) {

  if(!crux_set_command_line_options(argc, argv))
    return 1;

  run();
  if(skip_cleanup_flag != 1)
    sqtp.clean_up(out_dir);
    
  return 0;
}   

bool QRanker :: needsOutputDirectory()
{
  return true;
}


string QRanker::getName() {
  return "q-ranker";
}

string QRanker::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the new q-ranker algorithm.";

}

COMMAND_T QRanker::getCommand(){
  return QRANKER_COMMAND;
}
