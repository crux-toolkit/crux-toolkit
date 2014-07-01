#include "QRanker.h"
#include "modifications.h"

QRanker::QRanker() :  
  seed(0),
  selectionfdr(0.01),
  num_hu(4),
  mu(0.01),
  weightDecay(0.000),
  max_net_gen(NULL),
  max_net_targ(NULL),
  nets(NULL),
  in_dir(""), 
  out_dir(""), 
  skip_cleanup_flag(0),
  overwrite_flag(0),
  fileroot(""),
  feature_file_flag(0),
  file_format_("")
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
      if(fdr > max_fdr) {
        max_fdr = fdr;
        ind = count;
      }
    }
  net = max_net_targ[ind];
  getOverFDR(fullset,net, selectionfdr);
  d.clear_data_psm_training();
  
  d.load_data_psm_results();
  computePEP();

  ostringstream fname;
  if (get_boolean_parameter("txt-output")) {
    fname << out_dir << "/" << fileroot << "q-ranker.target.psms.txt";
    ofstream f1(fname.str().c_str()); 
    fname.str("");
    fname << out_dir << "/" << fileroot << "q-ranker.decoy.psms.txt";
    ofstream f2(fname.str().c_str());
    write_results_psm_tab(f1, f2);
    f1.close();
    f2.close();
    fname.str("");
  }

  if (get_boolean_parameter("pepxml-output")) {
    fname << out_dir << "/" << fileroot << "q-ranker.target.pep.xml";
    PepXMLWriter xmlfile1;
    xmlfile1.openFile(fname.str().c_str(), overwrite_flag);
    fname.str("");
    fname << out_dir << "/" << fileroot << "q-ranker.decoy.pep.xml";
    PepXMLWriter xmlfile2;
    xmlfile2.openFile(fname.str().c_str(), overwrite_flag);
    write_results_psm_xml(xmlfile1, xmlfile2);
    xmlfile1.closeFile();
    xmlfile2.closeFile();
  }
  d.clear_data_psm_results();
}

void QRanker :: write_results_psm_tab(ofstream &osTarget, ofstream &osDecoy)
{
  int ps =fullset[0].psmind; 
  ofstream* streams[2] = {&osTarget, &osDecoy};
  for (size_t i = 0; i < 2; ++i) {
    ofstream& os = *streams[i];
    os << "scan" << "\t" << "charge" << "\t" << "q-ranker q-value" << "\t" ;
    os<< "q-ranker score" << "\t" << "q-ranker PEP\t"; 
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
    os<< "filename" << endl;
  }
  for(int i = 0; i < fullset.size(); i++){
    ofstream& os = (fullset[i].label == 1) ? osTarget : osDecoy;
    int psmind = fullset[i].psmind;
    int pepind = d.psmind2pepind(psmind);
    string pep = d.ind2pep(pepind);
    string seq, n,c;
    os << d.psmind2scan(psmind) << "\t" ;
    os<< d.psmind2charge(psmind) << "\t"; 
    os<< fullset[i].q << "\t" ;
    os<< fullset[i].score << "\t"; 
    os<< fullset[i].PEP << "\t"; 
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
    //Flanking_aa 
    os<<n<<c;
    //TODO temp fix to write a set of flanking AAs per protein
    for(unsigned int j=1;j<prots.size();j++)
      os<<','<<n<<c;
    os<<'\t';
    os<< d.psmind2fname(psmind) << endl;
  }
}


void QRanker :: get_pep_seq(string &pep, string &seq, string &n, string &c)
{
  string tmp;
  int pos;
  pos = pep.find('.');
  n = pep.at(pos-1); 
  tmp = pep.substr(pos+1, pep.size());

  pos = tmp.rfind('.');
  c = tmp.at(pos+1);
  seq = tmp.substr(0, pos);
}


void QRanker ::write_results_psm_xml(
  PepXMLWriter& xmlfileTarget,
  PepXMLWriter& xmlfileDecoy)
{
  xmlfileTarget.writeHeader();
  xmlfileDecoy.writeHeader();

  bool* scores_to_print = new bool[NUMBER_SCORER_TYPES];
  for(int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
    scores_to_print[score_idx] = false;
  }
  scores_to_print[SP] = true; 
  scores_to_print[XCORR] = true;
  scores_to_print[QRANKER_SCORE] = true;
  scores_to_print[QRANKER_QVALUE] = true;
  scores_to_print[QRANKER_PEP] = true;

  xmlfileTarget.SetScoresComputed(scores_to_print);
  xmlfileDecoy.SetScoresComputed(scores_to_print);

  double* scores = new double[NUMBER_SCORER_TYPES];

  for(int i = 0; i < fullset.size(); i++)
    {
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
      int* psm_ranks = new int[NUMBER_SCORER_TYPES];
      psm_ranks[XCORR]=1; 
      double delta_cn = d.psmind2deltaCn(psmind);
      
      scores[XCORR] = d.psmind2xcorr(psmind);
      
      scores[SP] = d.psmind2spscore(psmind);
      scores[QRANKER_SCORE] = fullset[i].score;
      scores[QRANKER_QVALUE] = fullset[i].q;
      scores[QRANKER_PEP] = fullset[i].PEP;
      psm_ranks[SP]=-1; 
      
      psm_ranks[SP]=d.psmind2sp_rank(psmind);
      double by_ions_matched=d.psmind2by_ions_matched(psmind);
      
      psm_ranks[XCORR]=d.psmind2xcorr_rank(psmind);

      PepXMLWriter& xmlfile = (fullset[i].label == 1) ? xmlfileTarget : xmlfileDecoy;
      xmlfile.writePSM(scan, filename, spectrum_mass, charge, psm_ranks,
                       sequence, modified_sequence.c_str(),
                       peptide_mass, num_proteins,
                       flanking_aas.c_str(), protein_names, 
                       protein_descriptions, delta_cn, scores_to_print, scores,by_ions_matched,
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

  xmlfileTarget.writeFooter();
  xmlfileDecoy.writeFooter();
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
      ind = myrandom_limit(interval);
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
      ind = myrandom_limit(interval);
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
  //mysrandom(seed); This is set by CruxApplication::initialize()
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
string QRanker:: file_extension (string str){
  string file_exten= str.substr(str.rfind('.')+1);
  return file_exten; 
}


FILE_FORMAT_T QRanker::check_file_format(string& source) {
  string ext = file_extension(source);

    if (ext=="sqt") 
      return SQT_FORMAT;
    else if (ext=="txt") 
      return DELIMITED_FORMAT;
   return INVALID_FORMAT;
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
    "pepxml-output",
    "txt-output",
    "skip-cleanup",
    "re-run",
    "use-spec-features",
    "parameter-file",
    "verbosity",
     "list-of-files",
    "feature-file",
    "spectrum-parser"
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

/**
 * Uses the target and decoy scores to compute posterior error
 * probabilities. 
 */
void QRanker::computePEP(){
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


int QRanker::main(int argc, char **argv) {

  
  if(!crux_set_command_line_options(argc, argv))
    return 1;
  
  run();
  if(skip_cleanup_flag != 1)
    parser->clean_up(out_dir);
    
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
  return "Re-rank a collection of PSMs using the Q-ranker algorithm.";
}

COMMAND_T QRanker::getCommand(){
  return QRANKER_COMMAND;
}


void QRanker :: get_protein_id(int pepind, vector<string> &prot){
  int * protinds= d.pepind2protinds(pepind) ;
  unsigned prot_length=d.pepind2num_prot(pepind);  
  for (unsigned k=0;k<prot_length;k++){
    string protein_str= d.ind2prot( protinds[k]);
    prot.push_back(protein_str); 
  }
}


void QRanker :: print_protein_ids(vector<string> &prots, ofstream &os, int psmind){
    
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
