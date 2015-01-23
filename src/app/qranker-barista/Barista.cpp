#include "Barista.h"
#include "objects.h"
using namespace std; 
double Barista :: check_gradients_hinge_one_net(int protind, int label){
  int num_pep = d.protind2num_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  max_psm_inds.erase(max_psm_inds.begin(),max_psm_inds.end());
  max_psm_scores.erase(max_psm_scores.begin(),max_psm_scores.end());
  for (int i = 0; i < num_pep; i++){
    int pepind = pepinds[i];
    int num_psms = d.pepind2num_psm(pepind);
    int *psminds = d.pepind2psminds(pepind);
    double max_sc = -1000000.0;
    int max_ind = 0;
    for (int j = 0; j < num_psms; j++){
      double *feat = d.psmind2features(psminds[j]);
      double *sc = net.fprop(feat);
      if(sc[0] > max_sc){
        max_sc = sc[0];
	max_ind = j;
      }
    }
      max_psm_inds.push_back(max_ind);
      max_psm_scores.push_back(max_sc);
    }
  assert((int)max_psm_inds.size() == num_pep);
  assert((int)max_psm_scores.size() == num_pep);
  
  double sm = 0.0;
  double n = pow(num_pep,alpha);
  for(unsigned int i = 0; i < max_psm_inds.size() ; i++)
    sm+= max_psm_scores[i];
  sm /= n;

  //if(sm < 1)
  // {
  net.clear_gradients();
  double *gc = new double[1];
  gc[0] = -label/n;
  for(int i = 0; i < num_pep; i++){
    int pepind = pepinds[i];
    int *psminds = d.pepind2psminds(pepind);
    int max_psm_ind = psminds[max_psm_inds[i]];
    double *feat = d.psmind2features(max_psm_ind);
    net.fprop(feat);
    net.bprop(gc);
  }

  double h = 0.000001;
  double diff = -(1-sm*label);
  double err = 0.0;
  double *w = net.get_weights(1);
  double *dw = net.get_dweights(1);
  for (int k = 0; k < d.get_num_features(); k++){
      w[k] += h;
      max_psm_inds.erase(max_psm_inds.begin(),max_psm_inds.end());
      max_psm_scores.erase(max_psm_scores.begin(),max_psm_scores.end());
      for (int i = 0; i < num_pep; i++){
	  int pepind = pepinds[i];
	  int num_psms = d.pepind2num_psm(pepind);
	  int *psminds = d.pepind2psminds(pepind);
	  double max_sc = -1000000.0;
	  int max_ind = 0;
	  for (int j = 0; j < num_psms; j++){
	      double *feat = d.psmind2features(psminds[j]);
	      double *sc = net.fprop(feat);
	      if(sc[0] > max_sc){
		  max_sc = sc[0];
		  max_ind = j;
		}
	    }
	  max_psm_inds.push_back(max_ind);
	  max_psm_scores.push_back(max_sc);
	}
      double sm1 = 0.0;
      double n = pow(num_pep,alpha);
      for(unsigned int i = 0; i < max_psm_inds.size() ; i++)
	sm1 += max_psm_scores[i];
      sm1 /= n;
      
      diff += (1-sm1*label);
      //cout << dw[k] << " " << diff/h << " " << dw[k]-diff/h << endl;
      err += dw[k]-diff/h;
      w[k]-=h;
      diff -= (1-sm1*label);
    }

  
  double *bias = net.get_bias(1);
  double *dbias = net.get_dbias(1);
  for (int k = 0; k < num_hu; k++)
    {
      bias[k] += h;
      max_psm_inds.erase(max_psm_inds.begin(),max_psm_inds.end());
      max_psm_scores.erase(max_psm_scores.begin(),max_psm_scores.end());
      for (int i = 0; i < num_pep; i++)
	{
	  int pepind = pepinds[i];
	  int num_psms = d.pepind2num_psm(pepind);
	  int *psminds = d.pepind2psminds(pepind);
	  double max_sc = -1000000.0;
	  int max_ind = 0;
	  for (int j = 0; j < num_psms; j++)
	    {
	      double *feat = d.psmind2features(psminds[j]);
	      double *sc = net.fprop(feat);
	      if(sc[0] > max_sc)
		{
		  max_sc = sc[0];
		  max_ind = j;
		}
	    }
	  max_psm_inds.push_back(max_ind);
	  max_psm_scores.push_back(max_sc);
	}
      double sm1 = 0.0;
      double n = pow(num_pep,alpha);
      for(unsigned int i = 0; i < max_psm_inds.size() ; i++)
	sm1 += max_psm_scores[i];
      sm1 /= n;

      diff += 1-sm1*label;
      //cout << dbias[k] << " " << diff/h << " " << dbias[k]-diff/h << endl;
      err += dbias[k]-diff/h;
      bias[k]-=h;
      diff -= 1-sm1*label;
    }
  
  w = net.get_weights(2);
  dw = net.get_dweights(2);
  for (int k = 0; k < num_hu; k++)
    {
      w[k] += h;
      max_psm_inds.erase(max_psm_inds.begin(),max_psm_inds.end());
      max_psm_scores.erase(max_psm_scores.begin(),max_psm_scores.end());
      for (int i = 0; i < num_pep; i++)
	{
	  int pepind = pepinds[i];
	  int num_psms = d.pepind2num_psm(pepind);
	  int *psminds = d.pepind2psminds(pepind);
	  double max_sc = -1000000.0;
	  int max_ind = 0;
	  for (int j = 0; j < num_psms; j++)
	    {
	      double *feat = d.psmind2features(psminds[j]);
	      double *sc = net.fprop(feat);
	      if(sc[0] > max_sc)
		{
		  max_sc = sc[0];
		  max_ind = j;
		}
	    }
	  max_psm_inds.push_back(max_ind);
	  max_psm_scores.push_back(max_sc);
	}
      double sm1 = 0.0;
      double n = pow(num_pep,alpha);
      for(unsigned int i = 0; i < max_psm_inds.size() ; i++)
	sm1 += max_psm_scores[i];
      sm1 /= n;

      diff += 1-sm1*label;
      //cout << dw[k] << " " << diff/h << " " << dw[k]-diff/h << endl;
      err += dw[k]-diff/h;
      w[k]-=h;
      diff -= 1-sm1*label;
    }

  bias = net.get_bias(2);
  dbias = net.get_dbias(2);
  for (int k = 0; k < 1; k++)
    {
      bias[k] += h;
      max_psm_inds.erase(max_psm_inds.begin(),max_psm_inds.end());
      max_psm_scores.erase(max_psm_scores.begin(),max_psm_scores.end());
      for (int i = 0; i < num_pep; i++)
	{
	  int pepind = pepinds[i];
	  int num_psms = d.pepind2num_psm(pepind);
	  int *psminds = d.pepind2psminds(pepind);
	  double max_sc = -1000000.0;
	  int max_ind = 0;
	  for (int j = 0; j < num_psms; j++)
	    {
	      double *feat = d.psmind2features(psminds[j]);
	      double *sc = net.fprop(feat);
	      if(sc[0] > max_sc)
		{
		  max_sc = sc[0];
		  max_ind = j;
		}
	    }
	  max_psm_inds.push_back(max_ind);
	  max_psm_scores.push_back(max_sc);
	}
      double sm1 = 0.0;
      double n = pow(num_pep,alpha);
      for(unsigned int i = 0; i < max_psm_inds.size() ; i++)
	sm1 += max_psm_scores[i];
      sm1 /= n;

      diff += 1-sm1*label;
      //cout << dbias[k] << " " << diff/h << " " << dbias[k]-diff/h << endl;
      err += dbias[k]-diff/h;
      bias[k]-=h;
      diff -= 1-sm1*label;
    }
  //}
  return err;
}

double Barista :: check_gradients_hinge_clones(int protind, int label)
{
  double sm = get_protein_score(protind);
  net.clear_gradients();
  calc_gradients(protind,label);

  double h = 0.0000001;
  double diff = -(1-sm*label);
  double err = 0.0;
  double *w = net.get_weights(1);
  double *dw = net.get_dweights(1);
  for (int k = 0; k < d.get_num_features(); k++)
    {
      w[k] += h;
      double sm1 = get_protein_score(protind);
      diff += (1-sm1*label);
      //cout << dw[k] << " " << diff/h << " " << dw[k]-diff/h << endl;
      err += dw[k]-diff/h;
      w[k]-=h;
      diff -= (1-sm1*label);
    }
  
  double *bias = net.get_bias(1);
  double *dbias = net.get_dbias(1);
  for (int k = 0; k < num_hu; k++)
    {
      bias[k] += h;
      double sm1 = get_protein_score(protind);
      diff += 1-sm1*label;
      //cout << dbias[k] << " " << diff/h << " " << dbias[k]-diff/h << endl;
      err += dbias[k]-diff/h;
      bias[k]-=h;
      diff -= 1-sm1*label;
    }
  
  w = net.get_weights(2);
  dw = net.get_dweights(2);
  for (int k = 0; k < num_hu; k++)
    {
      w[k] += h;
      double sm1 = get_protein_score(protind);
      diff += 1-sm1*label;
      //cout << dw[k] << " " << diff/h << " " << dw[k]-diff/h << endl;
      err += dw[k]-diff/h;
      w[k]-=h;
      diff -= 1-sm1*label;
    }

  bias = net.get_bias(2);
  dbias = net.get_dbias(2);
  for (int k = 0; k < 1; k++)
    {
      bias[k] += h;
      double sm1 = get_protein_score(protind);
      diff += 1-sm1*label;
      //cout << dbias[k] << " " << diff/h << " " << dbias[k]-diff/h << endl;
      err += dbias[k]-diff/h;
      bias[k]-=h;
      diff -= 1-sm1*label;
    }
  
  net.update(mu);

  return err;
}

/**********************************************************/
int Barista :: getOverFDRPSM(PSMScores &s, NeuralNet &n,double fdr)
{
  double* featVec;
  for(int i = 0; i < s.size(); i++)
    {
      featVec = d.psmind2features(s[i].psmind);
      double *r = n.fprop(featVec);
      s[i].score = r[0];
    }

  int overFDR = s.calcOverFDR(fdr);
 
  if(verbose > 1)
    {
      cout << "psm over fdr: num psms " << overFDR << endl;
      int cn = 0;
      set<int> peptides;
      for(int i = 0; i < s.size(); i++)
	{
	  if(s[i].label == 1)
	    cn++;
	  peptides.insert(d.psmind2pepind(s[i].psmind));
	  if(cn >= overFDR)
	    break;
	}
      cout << "psm over fdr: num peptides " << peptides.size() << endl;
    }
  return overFDR;
}

double Barista :: get_peptide_score(int pepind, NeuralNet &n)
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


int Barista :: getOverFDRPep(PepScores &s, NeuralNet &n,double fdr)
{
  int pepind = 0;
  int label = 0;
  for(int i = 0; i < s.size(); i++)
    {
      pepind = s[i].pepind;
      double sc = get_peptide_score(pepind,n);
      s[i].score = sc;
    }

  int overFDR = s.calcOverFDR(fdr);

  if(verbose > 1)
    {
      int cn = 0;
      set<int> proteins_pos;
      set<int> proteins_neg;
      int pep_pos = 0;
      int pep_neg = 0;
      for(int i = 0; i < s.size(); i++)
	{
	  label = s[i].label;
	  if(label == 1)
	    cn++;
      
	  if(label == 1)
	    pep_pos++;
	  else
	    pep_neg++;
	  
	  int pepind = s[i].pepind;
	  int num_prot = d.pepind2num_prot(pepind);
	  int *protinds = d.pepind2protinds(pepind);
	  for(int j=0; j < num_prot;j++)
	    {
	      if(label == 1)
		proteins_pos.insert(protinds[j]);
	  else
	    proteins_neg.insert(protinds[j]);
	    }
	  if(cn >= overFDR)
	    break;
	}
      cout << "over FDR peptides:\n";
      cout << "num pos proteins  " << proteins_pos.size() << " num neg proteins " << proteins_neg.size() << " num pos peptides " << pep_pos << " num neg peptides " << pep_neg  << endl;

    }
  return overFDR;
}

/*****************************************************/
void Barista :: clear()
{
  trainset.clear();
  testset.clear();
  peptrainset.clear();
  peptestset.clear();
  psmtrainset.clear();
  psmtestset.clear();
  delete[] net_clones; net_clones = 0;
  max_psm_inds.clear();
  max_psm_scores.clear();
  used_peptides.clear();
  pepind_to_max_psmind.clear();
}


double Barista :: get_protein_score_parsimonious(int protind, NeuralNet &n)
{
  int num_pep = d.protind2num_pep(protind);
  int num_all_pep = d.protind2num_all_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  double sm = 0.0;
  double div = pow(num_all_pep,alpha);

  for (int i = 0; i < num_pep; i++)
    {
      int pepind = pepinds[i];

      if(used_peptides[pepind] == 0)
	{
	  used_peptides[pepind] = 1;
	  int num_psms = d.pepind2num_psm(pepind);
	  int *psminds = d.pepind2psminds(pepind);
	  double max_sc = -1000000.0;

	  for (int j = 0; j < num_psms; j++)
	    {
	      double *feat = d.psmind2features(psminds[j]);
	      double *sc = n.fprop(feat);
	      if(sc[0] > max_sc)
		{
		  max_sc = sc[0];
		}
	    }
	  sm += max_sc;
	}
    }
  
  sm /= div;
  return sm;
}

int Barista :: getOverFDRProtParsimonious(ProtScores &set, NeuralNet &n, double fdr)
{
  
  int total_num_pep = d.get_num_peptides();
  used_peptides.clear();
  used_peptides.resize(total_num_pep,0);
  double r = 0.0;
  for(int i = 0; i < set.size(); i++)
    {
      int protind = set[i].protind;
      r = get_protein_score_parsimonious(protind,n);
      set[i].score = r;
    }
  return set.calcOverFDR(fdr);
}

void Barista :: write_results_prot(string &out_dir, int fdr)
{
  ostringstream fn;
  fn << out_dir <<"/barista.target.proteins.txt";
  //cout << "writing results to " << fn.str() << endl;
  carp(CARP_INFO, "write_results to %s", fn.str().c_str());
  ofstream f_res(fn.str().c_str());

  int cn = 0;
  for(int i = 0; i < trainset.size(); i++)
    {
      if(trainset[i].label == 1)
	{
	  //write out proteins
	  f_res << "protein group " << cn+1  << " q=" << trainset[i].q << endl;
	  int protind = trainset[i].protind;
	  f_res << d.ind2prot(protind) << " ";
	  for(unsigned int j = 0; j < trainset[i].subset_protinds.size(); j++)
	    {
	      int ind = trainset[i].subset_protinds[j];
	      if(d.protind2label(ind) == 1)
		f_res << d.ind2prot(ind) << " ";
	    }
	  f_res << endl;
	
	  //write out peptides
	  f_res << "peptides\n";
	  int num_pep = d.protind2num_pep(protind);
	  int *pepinds = d.protind2pepinds(protind);
	  for( int j = 0; j < num_pep; j++)
	    {
	      int pepind = pepinds[j];
	      f_res << d.ind2pep(pepind); 
	      int psmind = pepind_to_max_psmind[pepind];
	      if(psmind > -1)
		f_res << "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
	      else
		cout << "warning: did not assign peptide max psmind\n";
	      f_res << " ";
	    }
	  f_res << endl << endl;
	  cn++;
	}
      if(cn >= fdr)
	break;
    }
  f_res.close();
}


void Barista :: report_all_results()
{
  trainset.clear();
  testset.clear();
  ProtScores::fillProteinsFull(trainset, d);
#ifdef CRUX
  carp(CARP_INFO, "finished training, making parsimonious protein set");
#else
  cout << "finished training, making parsimonious protein set\n";
#endif
  trainset.make_meta_set(d);
  int fdr_trn = getOverFDRProtParsimonious(trainset,max_net_prot,selectionfdr);
#ifdef CRUX
  carp(CARP_INFO, "total proteins parsimonious at q<%.2f: %d", selectionfdr, fdr_trn);
#else
  cout << "total proteins parsimonious at q<" << selectionfdr << ": " << fdr_trn << endl;
#endif

  int fdr_trn_psm = getOverFDRPSM(psmtrainset, max_net_psm, selectionfdr);
#ifdef CRUX
  carp(CARP_INFO, "peptides at q<%.2f: %d", selectionfdr, getOverFDRPep(peptrainset, max_net_pep, selectionfdr)); 
  carp(CARP_INFO, "psms at q<%.2f: %d", selectionfdr, fdr_trn_psm);
#else  
  cout << "peptides at q< " << selectionfdr << ": " <<  << endl;
  cout << "psms at q< " << selectionfdr << ": " << fdr_trn_psm << endl;
#endif
  
  write_results_prot(out_dir, fdr_trn);
}

void Barista :: get_pep_seq(string &pep, string &seq, string &n, string &c)
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

void Barista :: get_tab_delim_proteins(string protein_str, vector<string> &proteins)
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


void Barista :: write_protein_special_case_xml(ofstream &os, int i)
{

  string protein_str;
  vector<string> tab_delim_proteins;
  
  int protind = trainset[i].protind;
  int group = trainset[i].group_number;
  os << " <protein_group group_id=" << "\"" << group << "\"" << ">" << endl;
  os << "  <q_value>" << trainset[i].q << "</q_value>" << endl;
  os << "  <score>" << trainset[i].score << "</score>" << endl;
  os << "  <protein_ids>" << endl;
  protein_str = d.ind2prot(protind);
  get_tab_delim_proteins(protein_str, tab_delim_proteins);
  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
  
  vector<int> complement = trainset[i].protind2complement[protind];
  if(complement.size() != 0)
    {
      for(unsigned int k = 0; k < complement.size(); k++)
	{
	  int pepind = complement[k];
	  string pep = d.ind2pep(pepind);
	  string seq, n, c;
	  get_pep_seq(pep, seq, n, c);
	  os << "   <alternative_peptide_id>" << seq << "</alternative_peptide_id>"<< endl;
	}
    }
  
  for(unsigned int j = 0; j < trainset[i].indistinguishable_protinds.size(); j++)
    {
      int ind = trainset[i].indistinguishable_protinds[j];
      if(d.protind2label(ind) == 1)
	{
	  //os << "   <protein_id>" << d.ind2prot(ind) << "</protein_id> " << endl;;
	  protein_str = d.ind2prot(ind);
	  get_tab_delim_proteins(protein_str, tab_delim_proteins);
	  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
	    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
	  
	  vector<int> complement = trainset[i].protind2complement[ind];
	  if(complement.size() != 0)
	    {
	      for(unsigned int k = 0; k < complement.size(); k++)
		{
		  int pepind = complement[k];
		  string pep = d.ind2pep(pepind);
		  string seq, n, c;
		  get_pep_seq(pep, seq, n, c);
		  os << "   <alternative_peptide_id>" << seq << "</alternative_peptide_id>"<< endl;
		}
	    }
	}
    }
  os << "  </protein_ids>" << endl;
  
  //write out peptides
  os << "  <peptide_ids>" << endl;
  vector<int> intersection = trainset[i].intersection;
  for( unsigned int j = 0; j < intersection.size(); j++)
    {
      int pepind = intersection[j];
      string pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      os << "   <peptide_id>" << seq << "</peptide_id>"<< endl; 
    }
  os << "  </peptide_ids>" << endl;
  os << " </protein_group>" << endl;
}


void Barista :: write_results_prot_xml(ofstream &os)
{
  string protein_str;
  vector<string> tab_delim_proteins;
  os << "<proteins>" <<endl;
  int cn = 0;
  for(int i = 0; i < trainset.size(); i++)
    {
      int protind = trainset[i].protind;
      int group = trainset[i].group_number;
      if(d.protind2label(protind) == 1)
	{
	  cn++;
	  if(trainset[i].has_complement == 1)
	    {
	      write_protein_special_case_xml(os, i);
	      continue;
	    }
	  //write out proteins
	  os << " <protein_group group_id=" << "\"" << group << "\"" << ">" << endl;
	  os << "  <q_value>" << trainset[i].q << "</q_value>" << endl;
	  os << "  <score>" << trainset[i].score << "</score>" << endl;
	  os << "  <nsaf>" << trainset[i].nsaf << "</nsaf>" << endl;
          os<<"<barista PEP>"<<trainset[i].PEP<<"</barista PEP>"<<endl; 
	  os << "  <protein_ids>" << endl;
	  protein_str = d.ind2prot(protind);
	  get_tab_delim_proteins(protein_str, tab_delim_proteins);
	  //os << "   <protein_id>" << d.ind2prot(protind) << "</protein_id>" <<endl;
	  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
	    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
	  
	  for(unsigned int j = 0; j < trainset[i].indistinguishable_protinds.size(); j++)
	    {
	      int ind = trainset[i].indistinguishable_protinds[j];
	      if(d.protind2label(ind) == 1)
		{
		  protein_str = d.ind2prot(ind);
		  get_tab_delim_proteins(protein_str, tab_delim_proteins);
		  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
		    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
		}
	    }
	  os << "  </protein_ids>" << endl;
	
	  //write out peptides

	  os << "  <peptide_ids>" << endl;
	  int num_pep = d.protind2num_pep(protind);
	  int *pepinds = d.protind2pepinds(protind);
	  for( int j = 0; j < num_pep; j++)
	    {
	      int pepind = pepinds[j];
	      string pep = d.ind2pep(pepind);
	      string seq, n, c;
	      get_pep_seq(pep, seq, n, c);
	      os << "   <peptide_id>" << seq << "</peptide_id>"<< endl; 
	    }
	  os << "  </peptide_ids>" << endl;
	  os << " </protein_group>" << endl;
	}

    }
  os << "</proteins>" << endl;
}

void Barista :: write_subset_protein_special_case_xml(ofstream &os, int i)
{
  string protein_str;
  vector<string> tab_delim_proteins;

  //write out proteins
  int protind = trainset.get_subset_prot(i).protind;
  int group = trainset.get_subset_prot(i).group_number;
  os << " <protein_group group_id=" << "\"" << group << "\"";
  vector<int> parent_groups = trainset.get_subset_prot(i).parent_groups;
  if(parent_groups.size() == 0)
    cout << "problem: no parent groups" << endl;
  os << " parent_group_ids=" << "\""; 
  for(unsigned int k = 0; k < parent_groups.size()-1; k++)
    os << parent_groups[k] << ","; 
  os << parent_groups[parent_groups.size()-1] << "\"";
  os << ">" << endl;

  os << "  <protein_ids>" << endl;
  protein_str = d.ind2prot(protind);
  get_tab_delim_proteins(protein_str, tab_delim_proteins);
  //os << "   <protein_id>" << d.ind2prot(protind) << "</protein_id>" <<endl;
  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;

  vector<int> complement = trainset.get_subset_prot(i).protind2complement[protind];
  if(complement.size() != 0)
    {
      for(unsigned int k = 0; k < complement.size(); k++)
	{
	  int pepind = complement[k];
	  string pep = d.ind2pep(pepind);
	  string seq, n, c;
	  get_pep_seq(pep, seq, n, c);
	  os << "   <alternative_peptide_id>" << seq << "</alternative_peptide_id>"<< endl;
	}
    }
  
  for(unsigned int j = 0; j < trainset.get_subset_prot(i).indistinguishable_protinds.size(); j++)
    {
      int ind = trainset.get_subset_prot(i).indistinguishable_protinds[j];
      if(d.protind2label(ind) == 1)
	{
	  protein_str = d.ind2prot(ind);
	  get_tab_delim_proteins(protein_str, tab_delim_proteins);
	  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
	    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
	  
	  vector<int> complement = trainset.get_subset_prot(i).protind2complement[ind];
	  if(complement.size() != 0)
	    {
	      for(unsigned int k = 0; k < complement.size(); k++)
		{
		  int pepind = complement[k];
		  string pep = d.ind2pep(pepind);
		  string seq, n, c;
		  get_pep_seq(pep, seq, n, c);
		  os << "   <alternative_peptide_id>" << seq << "</alternative_peptide_id>"<< endl;
		}
	    }
	  
	}
    }
  os << "  </protein_ids>" << endl;
  
  //write out peptides
  os << "  <peptide_ids>" << endl;
  vector<int> intersection = trainset.get_subset_prot(i).intersection;
  for( unsigned int j = 0; j < intersection.size(); j++)
    {
      int pepind = intersection[j];
      string pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      os << "   <peptide_id>" << seq << "</peptide_id>"<< endl; 
    }
  os << "  </peptide_ids>" << endl;
  os << " </protein_group>" << endl;
  
}



void Barista :: write_results_subset_prot_xml(ofstream &os)
{
  string protein_str;
  vector<string> tab_delim_proteins;
  os << "<subset_proteins>" <<endl;

  int num_subset_prot = trainset.get_num_subsets();

  int cn = 0;
  for(int i = 0; i < num_subset_prot; i++)
    {
      int protind = trainset.get_subset_prot(i).protind;
      int group = trainset.get_subset_prot(i).group_number;
      if(d.protind2label(protind) == 1)
	{
	  cn++;

	  if(trainset.get_subset_prot(i).has_complement == 1)
	    {
	      write_subset_protein_special_case_xml(os, i);
	      continue;
	    }

	  //write out proteins
	  os << " <protein_group group_id=" << "\"" << group << "\"";
	  vector<int> parent_groups = trainset.get_subset_prot(i).parent_groups;
	  if(parent_groups.size() == 0)
	    cout << "problem: no parent groups" << endl;
	  os << " parent_group_ids=" << "\""; 
	  for(unsigned int k = 0; k < parent_groups.size()-1; k++)
	    os << parent_groups[k] << ","; 
	  os << parent_groups[parent_groups.size()-1] << "\"";
	  os << ">" << endl;
	  os << "  <protein_ids>" << endl;
	  protein_str = d.ind2prot(protind);
	  get_tab_delim_proteins(protein_str, tab_delim_proteins);
	  //os << "   <protein_id>" << d.ind2prot(protind) << "</protein_id>" <<endl;
	  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
	    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
	  
	  for(unsigned int j = 0; j < trainset.get_subset_prot(i).indistinguishable_protinds.size(); j++)
	    {
	      int ind = trainset.get_subset_prot(i).indistinguishable_protinds[j];
	      if(d.protind2label(ind) == 1)
		{
		  protein_str = d.ind2prot(ind);
		  get_tab_delim_proteins(protein_str, tab_delim_proteins);
		  for(unsigned int k = 0; k < tab_delim_proteins.size(); k++)
		    os << "   <protein_id>" << tab_delim_proteins[k] << "</protein_id>" <<endl;
		}
	    }
	  os << "  </protein_ids>" << endl;
	
	  //write out peptides
	  os << "  <peptide_ids>" << endl;
	  int num_pep = d.protind2num_pep(protind);
	  int *pepinds = d.protind2pepinds(protind);
	  for( int j = 0; j < num_pep; j++)
	    {
	      int pepind = pepinds[j];
	      string pep = d.ind2pep(pepind);
	      string seq, n, c;
	      get_pep_seq(pep, seq, n, c);
	      os << "   <peptide_id>" << seq << "</peptide_id>"<< endl; 
	    }
	  os << "  </peptide_ids>" << endl;
	  os << " </protein_group>" << endl;
	}
    }
  os << "</subset_proteins>" << endl;
}


void Barista :: write_results_peptides_xml(ofstream &os)
{
  string protein_str;
  vector<string> tab_delim_proteins;
  os << "<peptides>" <<endl;
  int cn = 0;
  for(int i = 0; i < peptrainset.size(); i++)
    {
      int pepind = peptrainset[i].pepind;
      if(peptrainset[i].label == 1)
	{
	  cn++;
	  //write out proteins
	  string pep = d.ind2pep(pepind);
	  string seq, n, c;
	  get_pep_seq(pep, seq, n, c);
	  os << " <peptide peptide_id=\"" << seq << "\">" << endl;
	  os << "  <q_value>" << peptrainset[i].q << "</q_value>" << endl;
	  os << "  <score>" << peptrainset[i].score << "</score>" << endl;
	  os << "  <nsaf>" << peptrainset[i].nsaf << "</nsaf>" << endl;
          os << " <barista PEP>" <<peptrainset[i].PEP<<"</barista PEP>"<<endl; 
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

void Barista :: write_results_psm_xml(ofstream &os)
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


void Barista :: write_results_pep_xml(PepXMLWriter& xmlfile)
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
      int* psm_ranks = new int[NUMBER_SCORER_TYPES];
      double delta_cn = d.psmind2deltaCn(psmind);
     
      scores[XCORR] = d.psmind2xcorr(psmind);
      scores[SP] = d.psmind2spscore(psmind);
      scores[QRANKER_SCORE] = psmtrainset[i].score;
      scores[QRANKER_QVALUE] = psmtrainset[i].q;
      scores[QRANKER_PEP] = psmtrainset[i].PEP;
      psm_ranks[SP]=d.psmind2sp_rank(psmind); 
      psm_ranks[XCORR]=d.psmind2xcorr_rank(psmind);
      double by_ion_matched=d.psmind2by_ions_matched(psmind);
      xmlfile.writePSM(scan, filename, spectrum_mass, charge, psm_ranks,
                       sequence, modified_sequence.c_str(),
                       peptide_mass, num_proteins,
                       flanking_aas.c_str(), protein_names, 
                       protein_descriptions, delta_cn, scores_to_print, scores,by_ion_matched, 
                       d.psmind2by_ions_total(psmind), 
                       d.psmind2matches_spectrum(psmind));

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



/*************************************************************************/

void Barista :: write_results_prot_special_case_tab(ofstream &os, int i)
{
  //write out protein
  int protind = trainset[i].protind;
  int group = trainset[i].group_number;
  os << group << "\t";
  os << trainset[i].q << "\t";
  os << trainset[i].score << "\t";
  os << d.ind2prot(protind); 
  vector<int> complement = trainset[i].protind2complement[protind];
  if(complement.size() != 0)
    {
      os << "(";
      int pepind = complement[0];
      string pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      os << pep;
      for(unsigned int k = 1; k < complement.size(); k++)
	{
	  pepind = complement[k];
	  pep = d.ind2pep(pepind);
	  get_pep_seq(pep, seq, n, c);
	  os << "," << pep;
	}
      os << ")";
    }
  
  for(unsigned int j = 0; j < trainset[i].indistinguishable_protinds.size(); j++)
    {
      int ind = trainset[i].indistinguishable_protinds[j];
      if(d.protind2label(ind) == 1)
	{
	  os << "," << d.ind2prot(ind);
	  vector<int> complement = trainset[i].protind2complement[protind];
	  if(complement.size() != 0)
	    {
	      os << "(";
	      int pepind = complement[0];
	      string pep = d.ind2pep(pepind);
	      string seq, n, c;
	      get_pep_seq(pep, seq, n, c);
	      int psmind = pepind_to_max_psmind[pepind];
	      if(psmind > -1)
		os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
	      for(unsigned int k = 1; k < complement.size(); k++)
		{
		  pepind = complement[k];
		  pep = d.ind2pep(pepind);
		  get_pep_seq(pep, seq, n, c);
		  psmind = pepind_to_max_psmind[pepind];
		  if(psmind > -1)
		    os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
      		}
	      os << ")";
	    }
	}
    }
  os << "\t";
  
  //write out peptides
  int num_pep = d.protind2num_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  //write the first one
  int pepind = pepinds[0];
  string pep = d.ind2pep(pepind);
  string seq, n, c;
  get_pep_seq(pep, seq, n, c);
  int psmind = pepind_to_max_psmind[pepind];
  if(psmind > -1)
    os << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
  else
    {
#ifndef CRUX
      cout << "warning: did not assign peptide max psmind\n";
#endif
    }
  for( int j = 1; j < num_pep; j++)
    {
      pepind = pepinds[j];
      pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      psmind = pepind_to_max_psmind[pepind];
      if(psmind > -1)
	os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
      else
	{
#ifndef CRUX
	  cout << "warning: did not assign peptide" << pep  << " ind " << pepind << " max psmind\n";
#endif
	}
    }
  os << endl;


}

void Barista :: write_results_prot_tab(ofstream &os)
{
  os << "group number" << "\t" << "q-value" << "\t" << "barista score" << "\t" << "NSAF score" << "\t";
  os<<"barista PEP\t";
  os << "proteins" << "\t" << "peptides-scan.charge"<<endl;

  int cn = 0;
  for(int i = 0; i < trainset.size(); i++)
    {
      if(trainset[i].label == 1)
	{
	  cn++;
	  if(trainset[i].has_complement == 1)
	    {
	      write_results_prot_special_case_tab(os, i);
	      continue;
	    }
	  //write out proteins
	  int protind = trainset[i].protind;
	  int group = trainset[i].group_number;
	  os << group << "\t";
	  os << trainset[i].q << "\t";
	  os << trainset[i].score << "\t";
	  os << trainset[i].nsaf << "\t";
          os<<trainset[i].PEP<<"\t"; 
	  os << d.ind2prot(protind); 
	  for(unsigned int j = 0; j < trainset[i].indistinguishable_protinds.size(); j++)
	    {
	      int ind = trainset[i].indistinguishable_protinds[j];
	      if(d.protind2label(ind) == 1)
		{
		  os << "," << d.ind2prot(ind);
		}
	    }
	  os << "\t";
	  
	  //write out peptides
	  int num_pep = d.protind2num_pep(protind);
	  int *pepinds = d.protind2pepinds(protind);
	  //write the first one
	  int pepind = pepinds[0];
	  string pep = d.ind2pep(pepind);
	  string seq, n, c;
	  get_pep_seq(pep, seq, n, c);
	  int psmind = pepind_to_max_psmind[pepind];
	  if(psmind > -1)
	    os << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
	  else
	    {
#ifndef CRUX
	      cout << "warning: did not assign peptide max psmind\n";
#endif
	    }
	  for( int j = 1; j < num_pep; j++)
	    {
	      pepind = pepinds[j];
	      pep = d.ind2pep(pepind);
	      string seq, n, c;
	      get_pep_seq(pep, seq, n, c);
	      psmind = pepind_to_max_psmind[pepind];
	      if(psmind > -1)
		os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
	      else
		{
#ifndef CRUX
		  cout << "warning: did not assign peptide" << pep  << " ind " << pepind << " max psmind\n";
#endif
		}
	    }
	  os << endl;
	}
    }
}


void Barista :: write_subset_protein_special_case_tab(ofstream &os, int i)
{
  //write out protein
  int protind = trainset[i].protind;
  int group = trainset[i].group_number;
  os << group << "\t";
  
  vector<int> parent_groups = trainset.get_subset_prot(i).parent_groups;
  if(parent_groups.size() > 0)
    {
      os << parent_groups[0]; 
      for(unsigned int k = 1; k < parent_groups.size(); k++)
	os << "," << parent_groups[k];
    } 
  os << "\t";

  os << d.ind2prot(protind); 
  vector<int> complement = trainset[i].protind2complement[protind];
  if(complement.size() != 0)
    {
      os << "(";
      int pepind = complement[0];
      string pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      os << pep;
      for(unsigned int k = 1; k < complement.size(); k++)
	{
	  pepind = complement[k];
	  pep = d.ind2pep(pepind);
	  get_pep_seq(pep, seq, n, c);
	  os << "," << pep;
	}
      os << ")";
    }
  
  for(unsigned int j = 0; j < trainset[i].indistinguishable_protinds.size(); j++)
    {
      int ind = trainset[i].indistinguishable_protinds[j];
      if(d.protind2label(ind) == 1)
	{
	  os << "," << d.ind2prot(ind);
	  vector<int> complement = trainset[i].protind2complement[protind];
	  if(complement.size() != 0)
	    {
	      os << "(";
	      int pepind = complement[0];
	      string pep = d.ind2pep(pepind);
	      string seq, n, c;
	      get_pep_seq(pep, seq, n, c);
	      int psmind = pepind_to_max_psmind[pepind];
	      if(psmind > -1)
		os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
	      for(unsigned int k = 1; k < complement.size(); k++)
		{
		  pepind = complement[k];
		  pep = d.ind2pep(pepind);
		  get_pep_seq(pep, seq, n, c);
		  psmind = pepind_to_max_psmind[pepind];
		  if(psmind > -1)
		    os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
      		}
	      os << ")";
	    }
	}
    }
  os << "\t";
  
  //write out peptides
  int num_pep = d.protind2num_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  //write the first one
  int pepind = pepinds[0];
  string pep = d.ind2pep(pepind);
  string seq, n, c;
  get_pep_seq(pep, seq, n, c);
  int psmind = pepind_to_max_psmind[pepind];
  if(psmind > -1)
    os << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
  else
    {
#ifndef CRUX
      cout << "warning: did not assign peptide max psmind\n";
#endif
    }
  for( int j = 1; j < num_pep; j++)
    {
      pepind = pepinds[j];
      pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      psmind = pepind_to_max_psmind[pepind];
      if(psmind > -1)
	os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
      else
	{
#ifndef CRUX
	  cout << "warning: did not assign peptide" << pep  << " ind " << pepind << " max psmind\n";
#endif
	}
    }
  os << endl;


}


void Barista :: write_results_subset_prot_tab(ofstream &os)
{
  string protein_str;
  vector<string> tab_delim_proteins;
  
  os << "group number" << "\t" << "parent group numbers" << "\t";
  os << "proteins" << "\t" << "peptides-scan.charge" << endl;

  int num_subset_prot = trainset.get_num_subsets();
  
  int cn = 0;
  for(int i = 0; i < num_subset_prot; i++)
    {
      int protind = trainset.get_subset_prot(i).protind;
      int group = trainset.get_subset_prot(i).group_number;
      if(d.protind2label(protind) == 1)
	{
	  cn++;
	  
	  if(trainset.get_subset_prot(i).has_complement == 1)
	    {
	      write_subset_protein_special_case_tab(os, i);
	      continue;
	    }
	  
	  os << group << "\t";
	  vector<int> parent_groups = trainset.get_subset_prot(i).parent_groups;
	  if(parent_groups.size() > 0)
	    {
	      os << parent_groups[0]; 
	      for(unsigned int k = 1; k < parent_groups.size(); k++)
		os << "," << parent_groups[k];
	    } 
	  os << "\t";
	  
	  //write out proteins
	  os << d.ind2prot(protind); 
	  for(unsigned int j = 0; j < trainset[i].indistinguishable_protinds.size(); j++)
	    {
	      int ind = trainset[i].indistinguishable_protinds[j];
	      if(d.protind2label(ind) == 1)
		{
		  os << "," << d.ind2prot(ind);
		}
	    }
	  os << "\t";
	  
	  //write out peptides
	  int num_pep = d.protind2num_pep(protind);
	  int *pepinds = d.protind2pepinds(protind);
	  //write the first one
	  int pepind = pepinds[0];
	  string pep = d.ind2pep(pepind);
	  string seq, n, c;
	  get_pep_seq(pep, seq, n, c);
	  int psmind = pepind_to_max_psmind[pepind];
	  if(psmind > -1)
	    os << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
	  else
	    {
#ifndef CRUX
	      cout << "warning: did not assign peptide max psmind\n";
#endif
	    }
	  for( int j = 1; j < num_pep; j++)
	    {
	      pepind = pepinds[j];
	      pep = d.ind2pep(pepind);
	      string seq, n, c;
	      get_pep_seq(pep, seq, n, c);
	      psmind = pepind_to_max_psmind[pepind];
	      if(psmind > -1)
		os << "," << pep <<  "-" << d.psmind2scan(psmind) << "." << d.psmind2charge(psmind);
	      else
		{
#ifndef CRUX
		  cout << "warning: did not assign peptide" << pep  << " ind " << pepind << " max psmind\n";
#endif
		}
	    }
	  
	  os << endl;
	}
    }

}



void Barista :: write_results_peptides_tab(ofstream &os)
{
  int ps = pepind_to_max_psmind[peptrainset[0].pepind]; 
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
  for(int i = 0; i < peptrainset.size(); i++){
    if(peptrainset[i].label == 1){
      cn++;
      //write out proteins
      int pepind = peptrainset[i].pepind;
      string pep = d.ind2pep(pepind);
      string seq, n, c;
      get_pep_seq(pep, seq, n, c);
      //q-value
      os << peptrainset[i].q << "\t";
      //score 
      os << peptrainset[i].score << "\t";
      ///nsaf
      os<<peptrainset[i].nsaf<<"\t";
     //PEP
      os<<peptrainset[i].PEP<<"\t"; 
   
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

void Barista :: write_results_psm_tab(ofstream &os)
{
  int ps= psmtrainset[0].psmind; 
  os << "scan" << "\t" << "charge" << "\t";
  os << "barista q-value" << "\t" << "barista score" << "\t";
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



void Barista :: report_all_results_xml()
{
  ostringstream fname;
  fname << out_dir << "/" << fileroot << "barista.xml";
  ofstream of(fname.str().c_str());
  of << "<barista_output>" << endl;
  of << endl;
  
  write_results_prot_xml(of);
  write_results_subset_prot_xml(of);
  
  write_results_peptides_xml(of);

  write_results_psm_xml(of);

  if (get_boolean_parameter("pepxml-output")) {
    ostringstream xml_file_name;
    xml_file_name << out_dir << "/" << fileroot << "barista.target.pep.xml";
    PepXMLWriter xmlfile;
    xmlfile.openFile(xml_file_name.str().c_str(), overwrite_flag);
    write_results_pep_xml(xmlfile);

    xmlfile.closeFile();
  }

  
  of << endl;
  of << "</barista_output>" << endl;
  of.close();
    
}


void Barista :: report_all_results_tab()
{
  ofstream of;
  ostringstream fname;
    
  fname << out_dir << "/" << fileroot << "barista.target.proteins.txt";
  of.open(fname.str().c_str());
  write_results_prot_tab(of);
  of.close();
  fname.str("");

  fname << out_dir << "/" << fileroot << "barista.target.subset-proteins.txt";
  of.open(fname.str().c_str());
  write_results_subset_prot_tab(of);
  of.close();
  fname.str("");
  
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


void Barista :: report_all_results_xml_tab()
{
  setup_for_reporting_results();
  stringstream fname;
  report_all_results_xml();
  if (get_boolean_parameter("txt-output")) {
    report_all_results_tab();
  }
  d.clear_data_all_results(); 
  
}

void Barista :: setup_for_reporting_results()
{
  d.load_labels_prot_training();
  trainset.clear();
  testset.clear();

  ProtScores::fillProteinsFull(trainset, d);
#ifdef CRUX
  carp(CARP_INFO, "finished training, making parsimonious protein set");
#else
  cout << "finished training, making parsimonious protein set\n";
#endif
  trainset.make_meta_set(d);
  
  psmtrainset.clear();
  psmtestset.clear();
  PSMScores::fillFeaturesFull(psmtrainset, d);
  peptrainset.clear();
  peptestset.clear();
  PepScores::fillFeaturesFull(peptrainset, d);
  d.clear_labels_prot_training();
  
  int fdr_trn = getOverFDRProt(trainset,max_net_prot,selectionfdr);
  fdr_trn = getOverFDRProtParsimonious(trainset,max_net_prot,selectionfdr);
#ifdef CRUX
  carp(CARP_INFO, "total proteins parsimonious at q<%.2f: %d", selectionfdr, fdr_trn);
#else
  cout << "total proteins parsimonious at q<" << selectionfdr << ": " << fdr_trn << endl;
#endif
  int fdr_trn_psm = getOverFDRPSM(psmtrainset, max_net_psm, selectionfdr);

  pepind_to_max_psmind.clear();
  pepind_to_max_psmind.resize(d.get_num_peptides(),-1);
  int fdr_trn_pep = getOverFDRPep(peptrainset, max_net_pep, selectionfdr);
#ifdef CRUX
  carp(CARP_INFO, "peptides at q<%.2f: %d", selectionfdr, fdr_trn_pep);
  carp(CARP_INFO, "psms at q<%.2f: %d", selectionfdr, fdr_trn_psm);
#else 
  cout << "peptides at q< " << selectionfdr << ": " << getOverFDRPep(peptrainset, max_net_pep, selectionfdr) << endl;
  cout << "psms at q< " << selectionfdr << ": " << fdr_trn_psm << endl;
#endif
  d.clear_labels_psm_training();
  
  d.clear_data_prot_training();
  d.load_data_all_results();
    
  computeNSAF();
  computePepNSAF();
  computePEP();
  
}

int Barista :: computeNSAF()
{
  //check that the lengths of the proteins are available
  int cnt = 0;
  for(int i = 0; i < trainset.size(); i++)
    {
      if(trainset[i].label == 1)
	{
	  int protind = trainset[i].protind;
	  int len = d.protind2length(protind);
	  if(len < 0)
	    cnt++;
	}
    }
  if(cnt > 0)
    return 0;

  //create auxiliary psmind_to_ind index
  vector<int>psmind_to_ind;
  psmind_to_ind.resize(psmtrainset.size(),0);
  for(int i = 0; i < psmtrainset.size(); i++)
    psmind_to_ind[psmtrainset[i].psmind] = i;

  //compute NSAF for each protein
  double sum = 0.0;
  for(int k = 0; k < trainset.size(); k++)
    {
      if(trainset[k].label == 1)
	{
	  int protind = trainset[k].protind;
	  int len = d.protind2length(protind);
	  int num_pep = d.protind2num_pep(protind);
	  int *pepinds = d.protind2pepinds(protind);
	  int cnt = 0;
	  for (int i = 0; i < num_pep; i++)
	    {
	      int pepind = pepinds[i];
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

int Barista :: computePepNSAF()
{

  //create auxiliary psmind_to_ind index
  vector<int>psmind_to_ind;
  psmind_to_ind.resize(psmtrainset.size(),0);
  for(int i = 0; i < psmtrainset.size(); i++)
    psmind_to_ind[psmtrainset[i].psmind] = i;

  //compute NSAF for each peptide
  double sum = 0.0;
  for(int k = 0; k < peptrainset.size(); k++)
    {
      if(peptrainset[k].label == 1)
	{
	  int pepind = peptrainset[k].pepind;
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
	  peptrainset[k].nsaf = nsaf;
	  sum += nsaf;
 	}
    }
  
  if(sum > 0)
    {
      for(int i = 0; i < peptrainset.size(); i++)
	{
	  if(peptrainset[i]. label == 1)
	    peptrainset[i].nsaf /= sum;
	}
    }
  else
    return 0;

  /*
  double check = 0.0;
  for(int i = 0; i < peptrainset.size(); i++)
    {
      if(peptrainset[i]. label == 1)
	check += peptrainset[i].nsaf;
    }
  cout << check << endl;
  */
  psmind_to_ind.clear();
  return 1;
}


void Barista :: computePEP(){
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

  delete[] target_scores;
  delete[] decoy_scores;
  delete[] PEPs;

  /** 
   * Calculate Peptide PEPs
   */
  target_scores_vect.clear();
  decoy_scores_vect.clear();
  
  //pull out the target and decoy scores for peptide PEPs
  for(int i=0;i<peptrainset.size();i++){
    if(peptrainset[i].label==1)
      target_scores_vect.push_back(peptrainset[i].score);
    else 
      decoy_scores_vect.push_back(peptrainset[i].score);
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
  for(int full_idx = 0; full_idx < psmtrainset.size(); full_idx++){
    if( peptrainset[full_idx].label == 1 ){
      peptrainset[full_idx].PEP = PEPs[target_idx];
      target_idx++; 
    } // else, skip decoys
  }

  delete[] target_scores;
  delete[] decoy_scores;
  delete[] PEPs;


  /** 
   * Calculate Protein PEPs
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
  carp(CARP_DEBUG, "Found %d targets and %d decoys of Protein", num_targets, num_decoys);  

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

  delete[] target_scores;
  delete[] decoy_scores;
  delete[] PEPs;


}

/*************************************************************************/
void Barista :: report_prot_fdr_counts(vector<double> &qvals, ofstream &of)
{
  for(unsigned int count = 0; count < qvals.size(); count++)
    {
      double q = qvals[count];
      if(trainset.size() > 0)
	{
	  int fdr_trn = getOverFDRProtParsimonious(trainset,max_net_prot,q);
	  of << q << " " << fdr_trn;
	  cout << q << " " << fdr_trn;
	}
      
      if(testset.size() > 0)
	{
	  int fdr_tst = getOverFDRProtParsimonious(testset,max_net_prot,q);
	  of << " " << fdr_tst;
	  cout << " " << fdr_tst;
	}
      of << endl;
      cout << endl;
    }
} 

void Barista :: report_psm_fdr_counts(vector<double> &qvals, ofstream &of)
{
  for(unsigned int count = 0; count < qvals.size(); count++)
    {
      double q = qvals[count];
      if(psmtrainset.size() > 0)
	{
	  int fdr_trn = getOverFDRPSM(psmtrainset,max_net_prot,q);
	  of << q << " " << fdr_trn;
	  cout << q << " " << fdr_trn;
	}
      
      if(psmtestset.size() > 0)
	{
	  int fdr_tst = getOverFDRPSM(psmtestset,max_net_prot,q);
	  of << " " << fdr_tst;
	  cout << " " << fdr_tst;
	}
      of << endl;
      cout << endl;
    }
} 


void Barista :: report_pep_fdr_counts(vector<double> &qvals, ofstream &of)
{
  for(unsigned int count = 0; count < qvals.size(); count++)
    {
      double q = qvals[count];
      if(peptrainset.size() > 0)
	{
	  int fdr_trn = getOverFDRPep(peptrainset,max_net_prot,q);
	  of << q << " " << fdr_trn;
	  cout << q << " " << fdr_trn;
	}
      
      if(psmtestset.size() > 0)
	{
	  int fdr_tst = getOverFDRPep(peptestset,max_net_prot,q);
	  of << " " << fdr_tst;
	  cout << " " << fdr_tst;
	}
      of << endl;
      cout << endl;
    }
} 



void Barista :: report_all_fdr_counts()
{
  int num_qvals = 14;
  vector<double> qvals;
  qvals.resize(num_qvals,0.0);
  double q = 0.0;
  for(int count = 0; count < num_qvals; count++)
    {
      qvals[count] = q;
      if(q < 0.01)
	q+=0.0025;
      else
	q+=0.01;
    }
    
  ostringstream fname;
  fname << out_dir << "/barista_prot_m7" << ".txt";
  ofstream of(fname.str().c_str());
  if(trainset.size() > 0)
    trainset.make_meta_set(d);
  if(testset.size() > 0)
    testset.make_meta_set(d);
  report_prot_fdr_counts(qvals,of);
  of.close();

  fname.str("");
  fname << out_dir << "/barista_psm_m7" << ".txt";
  ofstream ofpsm(fname.str().c_str());
  report_psm_fdr_counts(qvals,ofpsm);
  ofpsm.close();

  fname.str("");
  fname << out_dir << "/barista_pep_m7" << ".txt";
  ofstream ofpep(fname.str().c_str());
  report_pep_fdr_counts(qvals,ofpep);
  ofpep.close();


}



/*******************************************************************************/
double Barista :: get_protein_score(int protind, NeuralNet &n)
{
  int num_pep = d.protind2num_pep(protind);
  int num_all_pep = d.protind2num_all_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  double sm = 0.0;
  double div = pow(num_all_pep,alpha);

  for (int i = 0; i < num_pep; i++)
    {
      int pepind = pepinds[i];
      int num_psms = d.pepind2num_psm(pepind);
      int *psminds = d.pepind2psminds(pepind);
      double max_sc = -1000000.0;

      for (int j = 0; j < num_psms; j++)
	{
	  double *feat = d.psmind2features(psminds[j]);
	  double *sc = n.fprop(feat);
	  if(sc[0] > max_sc)
	    {
	      max_sc = sc[0];

	    }
	}
      sm += max_sc;
    }
  sm /= div;
  return sm;
}

int Barista :: getOverFDRProt(ProtScores &set, NeuralNet &n, double fdr)
{
  double r = 0.0;
  for(int i = 0; i < set.size(); i++)
    {
      int protind = set[i].protind;
      r = get_protein_score(protind,n);
      set[i].score = r;
    }
  return set.calcOverFDR(fdr);
  
}

int Barista :: getOverFDRProt(ProtScores &set, double fdr)
{
  double r = 0.0;
  for(int i = 0; i < set.size(); i++)
    {
      int protind = set[i].protind;
      r = get_protein_score(protind);
      set[i].score = r;
    }
  return set.calcOverFDR(fdr);
  
}

double Barista :: get_protein_score(int protind)
{
  int num_pep = d.protind2num_pep(protind);
  int num_all_pep = d.protind2num_all_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  max_psm_inds.erase(max_psm_inds.begin(),max_psm_inds.end());
  max_psm_scores.erase(max_psm_scores.begin(),max_psm_scores.end());
  int psm_count = 0;
  for (int i = 0; i < num_pep; i++)
    {
      int pepind = pepinds[i];
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
      max_psm_scores.push_back(max_sc);
    }
  assert((int)max_psm_inds.size() == num_pep);
  assert((int)max_psm_scores.size() == num_pep);
  
  double sm = 0.0;
  double n = pow(num_all_pep,alpha);
  for(unsigned int i = 0; i < max_psm_inds.size() ; i++)
    sm+= max_psm_scores[i];
  sm /= n;
  return sm;
}

void Barista :: calc_gradients(int protind, int label)
{
  int num_pep = d.protind2num_pep(protind);
  int num_all_pep = d.protind2num_all_pep(protind);
  int *pepinds = d.protind2pepinds(protind);
  double n = pow(num_all_pep,alpha);

  double *gc = new double[1];
  gc[0] = -label/n;
  int psm_count = 0;
  for(int i = 0; i < num_pep; i++)
    {
      int pepind = pepinds[i];
      int num_psms = d.pepind2num_psm(pepind);
      int clone_ind = psm_count+max_psm_inds[i];
      net_clones[clone_ind].bprop(gc);
      psm_count += num_psms;
    }
  delete[] gc;
}

double Barista :: train_hinge(int protind, int label)
{
  double sm = get_protein_score(protind);
    
  double err = max(0.0,1.0-sm*label);

  if(sm*label < 1)
    {
      net.clear_gradients();
      calc_gradients(protind,label);
      net.update(mu);
    }
  return err;
}


double Barista :: train_hinge_psm(int psmind, int label)
{
  double *x = d.psmind2features(psmind);
  double *c = net.fprop(x);
  double err = 1.0-c[0]*label;

  if(c[0]*label < 1)
    {
      double *gc = new double[1];
      gc[0] = -1*label;
      net.clear_gradients();
      net.bprop(gc);
      net.update(mu);
      delete[] gc;
    }
  return err;
}

void Barista :: train_net(double selectionfdr)
{
  for (int k = 0; k < nepochs; k++)
    {
      //fdebug << "epoch " << k << endl;
    if(verbose > 0)
	cout << "epoch " << k << endl;
      double err_sum = 0.0;
      for(int i = 0; i < trainset.size(); i++)
	{
	  int ind = myrandom_limit(trainset.size());
	  int protind = trainset[ind].protind;
	  int label = trainset[ind].label;
	  err_sum += train_hinge(protind,label);
	}      
      //int fdr_trn = getOverFDRProt(trainset,net,selectionfdr);
      int fdr_trn = getOverFDRProt(trainset,selectionfdr);
      
      //fdebug << fixed;
      //for(int j = 0; j < trainset.size(); j++)
	//fdebug << j << " " << trainset[j].protind << " " << setprecision(16) << trainset[j].score << endl;
	//fdebug << j << " " << trainset[j].protind  << endl;
      
      if(verbose > 0)
	{
	  cout << "err " << err_sum << "  ";
	  cout << selectionfdr << " " << fdr_trn;
	  if(testset.size() > 0)
	    cout << " " << getOverFDRProt(testset,net,selectionfdr);
	  cout << endl;
	}
	if(fdr_trn > max_fdr)
	{
	  max_net_prot = net;
	  max_fdr = fdr_trn;
	  if(verbose == 0)
	    {
#ifdef CRUX
	      if(testset.size() > 0)
		carp(CARP_INFO, "q<%.2f: max non-parsimonious so far %d %d", selectionfdr, max_fdr, getOverFDRProt(testset,max_net_prot,selectionfdr));
	      else
		carp(CARP_INFO, "q<%.2f: max non-parsimonious so far %d", selectionfdr, max_fdr);
#else
	      cout << "q< " << selectionfdr << ": max non-parsimonious so far " << max_fdr;
	      if(testset.size() > 0)
		cout << " " << getOverFDRProt(testset,max_net_prot,selectionfdr);
	      cout << endl;
#endif
	    }
	}
      if(verbose > 0)
	{
	  cout << "q< " << selectionfdr << ": max non-parsimonious so far " << max_fdr;
	  if(testset.size() > 0)
	    cout << " " << getOverFDRProt(testset,max_net_prot,selectionfdr);
	  cout << endl;
	}
      if(1)
	{
	  int fdr_trn_psm = getOverFDRPSM(psmtrainset, net, selectionfdr); 
	  if(fdr_trn_psm > max_fdr_psm)
	    {
	      max_net_psm = net;
	      max_fdr_psm = fdr_trn_psm;
	    }
	}
      if(1)
	{
	  int fdr_trn_pep = getOverFDRPep(peptrainset, net, selectionfdr); 
	  if(fdr_trn_pep > max_fdr_pep)
	    {
	      max_net_pep = net;
	      max_fdr_pep = fdr_trn_pep;
	    }
	}
    }
  if(verbose > 0)
    {
      cout << "max peptides so far at q< " << selectionfdr << ": " << getOverFDRPep(peptrainset, max_net_psm, selectionfdr) << endl;
      cout << "max psms so far at q<" << selectionfdr << ": " << max_fdr_psm << endl;
    }
}



void Barista :: train_net_multi_task(double selectionfdr, int interval)
{
  for (int k = 0; k < nepochs; k++)
    {
      //fdebug << "epoch " << k << endl;
    if(verbose > 0)
	cout << "epoch " << k << endl;
      double err_sum = 0.0;
      for(int i = 0; i < trainset.size(); i++)
	{
	  int ind = myrandom_limit(trainset.size());
	  int protind = trainset[ind].protind;
	  int label = trainset[ind].label;
	  err_sum += train_hinge(protind,label);

	  ind = myrandom_limit(interval);
	  int psmind = psmtrainset[ind].psmind;
	  label = psmtrainset[ind].label;
	  train_hinge_psm(psmind,label);	
	}      
      //int fdr_trn = getOverFDRProt(trainset,net,selectionfdr);
      int fdr_trn = getOverFDRProt(trainset,selectionfdr);
      
      //fdebug << fixed;
      //for(int j = 0; j < trainset.size(); j++)
	//fdebug << j << " " << trainset[j].protind << " " << setprecision(16) << trainset[j].score << endl;
	//fdebug << j << " " << trainset[j].protind  << endl;
      
      if(verbose > 0)
	{
	  if(interval == trainset.size())
	    cout << "err " << err_sum << "  ";
	  cout << selectionfdr << " " << fdr_trn;
	  if(testset.size() > 0)
	    cout << " " << getOverFDRProt(testset,net,selectionfdr);
	  cout << endl;
	}
	if(fdr_trn > max_fdr)
	{
	  max_net_prot = net;
	  max_fdr = fdr_trn;
	  if(verbose == 0)
	    {
#ifdef CRUX
	      if(testset.size() > 0)
		carp(CARP_INFO, "q<%.2f: max non-parsimonious so far %d %d", selectionfdr, max_fdr, getOverFDRProt(testset,max_net_prot,selectionfdr));
	      else
		carp(CARP_INFO, "q<%.2f: max non-parsimonious so far %d", selectionfdr, max_fdr);
#else
	      cout << "q< " << selectionfdr << ": max non-parsimonious so far " << max_fdr;
	      if(testset.size() > 0)
		cout << " " << getOverFDRProt(testset,max_net_prot,selectionfdr);
	      cout << endl;
#endif
	    }
	}
      if(verbose > 0)
	{
	  cout << "q< " << selectionfdr << ": max non-parsimonious so far " << max_fdr;
	  if(testset.size() > 0)
	    cout << " " << getOverFDRProt(testset,max_net_prot,selectionfdr);
	  cout << endl;
	}
      //cout << "max prot " << getOverFDRProtMax(trainset,max_net_prot,selectionfdr) << endl;
      if(1)
	{
	  int fdr_trn_psm = getOverFDRPSM(psmtrainset, net, selectionfdr); 
	  if(fdr_trn_psm > max_fdr_psm)
	    {
	      max_net_psm = net;
	      max_fdr_psm = fdr_trn_psm;
	    }
	}
      if(1)
	{
	  int fdr_trn_pep = getOverFDRPep(peptrainset, net, selectionfdr); 
	  if(fdr_trn_pep > max_fdr_pep)
	    {
	      max_net_pep = net;
	      max_fdr_pep = fdr_trn_pep;
	    }
	}
    }
  if(verbose > 0)
    {
      cout << "max peptides so far at q< " << selectionfdr << ": " << getOverFDRPep(peptrainset, max_net_psm, selectionfdr) << endl;
      cout << "max psms so far at q<" << selectionfdr << ": " << max_fdr_psm << endl;
    }
}



void Barista :: setup_for_training(int trn_to_tst)
{
#ifdef CRUX
  carp(CARP_INFO, "loading and normalizing data");
#else
  cout << "loading data" << endl;
#endif
   
  d.load_data_prot_training();
  d.load_labels_prot_training();
  
  if(trn_to_tst)
    {
#ifdef CRUX
      carp(CARP_INFO, "splitting data into training and testing sets");
#else
      cout << "splitting data into training and testing sets" << endl;
#endif
      ProtScores::fillProteinsSplit(trainset, testset, d, trn_to_tst);
#ifdef CRUX
      carp(CARP_INFO, "trainset size %d testset size %d", trainset.size(), testset.size());
#else
      cout << "trainset size " << trainset.size() << " testset size " << testset.size() << endl;
#endif
    }
  else
    {
      ProtScores::fillProteinsFull(trainset, d);
#ifdef CRUX
      carp(CARP_INFO, "trainset size %d", trainset.size());
#else
      cout << "trainset size " << trainset.size() << endl;
#endif
    }
  thresholdset = trainset;
  

  if(trn_to_tst)
    PSMScores::fillFeaturesSplit(psmtrainset, psmtestset, d, 0.5);
  else
    PSMScores::fillFeaturesFull(psmtrainset, d);
  if(trn_to_tst)
    PepScores::fillFeaturesSplit(peptrainset, peptestset, d, 0.5);
  else
    PepScores::fillFeaturesFull(peptrainset, d);
  d.clear_labels_prot_training();
  
  if(feature_file_flag)
    {
      string str = feature_file_name.str();
      if(!d.print_features(str))
	carp(CARP_INFO, "could not open file %s for writing features. Feature file will not be written", feature_file_name.str().c_str());
    }

  d.normalize_psms();

  num_features = d.get_num_features();
  int has_bias = 1;
  int is_lin = 1;
  if(is_lin)
    num_hu = 1;
  net.initialize(num_features, num_hu, is_lin, has_bias);
  max_net_prot = net;
  max_net_psm = net;
  max_net_pep = net;
  
  //get the max num peptides
  max_peptides = 0;
  max_psms_in_prot = 0;
  for (int i = 0; i < d.get_num_proteins();i++)
    {
      int num_pep = d.protind2num_pep(i);
      if(num_pep > max_peptides)
	max_peptides = num_pep;
      int *pepinds = d.protind2pepinds(i);
      int sum_num_psms = 0;
      for(int j = 0; j < num_pep; j++)
	{
	  int pepind = pepinds[j];
	  sum_num_psms += d.pepind2num_psm(pepind);
	}
      if(max_psms_in_prot < sum_num_psms)
	max_psms_in_prot = sum_num_psms;
    }
  max_psm_inds.reserve(max_peptides);
  max_psm_scores.reserve(max_peptides);
 
  //create the net clones
  net_clones = new NeuralNet[max_psms_in_prot];
  for (int i = 0; i < max_psms_in_prot;i++)
    net_clones[i].clone(net);
}


int Barista :: run()
{
  //mysrandom(seed); This is set by CruxApplication::initialize()
  setup_for_training(0);
  
  train_net(selectionfdr);
  
  return 0;

}


int Barista :: run_tries()
{
  setup_for_training(0);
    
  int tries = 3;
  vector<double> mu_choices;
  mu_choices.resize(3,0.0);
  mu_choices[0] = 0.005;
  mu_choices[1] = 0.01;
  mu_choices[2] = 0.05;
  for(int k = 0; k < tries; k++)
    {
      mu = mu_choices[k];
      net.make_random();
      train_net(selectionfdr);
    }
  
  report_all_results_xml_tab();

  return 0;

}

int Barista :: run_tries_multi_task()
{
  setup_for_training(0);
    
  //fdebug.open("debug.txt");

#ifdef CRUX
  carp(CARP_INFO, "training the model");
#else
  cout << "training the model" << endl;
#endif

  int tries = 3;
  vector<double> mu_choices;
  mu_choices.resize(3,0.0);
  mu_choices[0] = 0.005;
  mu_choices[1] = 0.01;
  mu_choices[2] = 0.05;
  for(int k = 1; k < tries; k++)
    {
      mu = mu_choices[k];
      net.make_random();
      train_net_multi_task(selectionfdr, psmtrainset.size());
    }
    
  //net.copy(max_net_prot);

  int interval= getOverFDRPSM(psmtrainset,max_net_psm,0.01)*2;
  if(interval > psmtrainset.size())
    interval = psmtrainset.size()/4;
  if(interval < 50)
    interval = psmtrainset.size()/4;
  train_net_multi_task(selectionfdr, interval);

  //report_all_fdr_counts();
  report_all_results_xml_tab();
  
  //fdebug.close();

   return 0;

}



void Barista :: print_description()
{
  cout << endl;
  cout << "\t crux barista [options] <protein database> <spectra> <search results>" << endl <<endl;
  cout << "REQUIRED ARGUMENTS:" << endl << endl;
  cout << "\t <protein database> Directory with FASTA files , list of FASTA files or a single FASTA file with the protein database used for the search." << endl;
  cout << "\t <spectra> Directory with ms2 files, list of ms2 files or a single ms2 file used for database search." << endl;
  cout << "\t <search results> Directory with sqt files, list of sqt files or a single sqt file with psms generated during search." << endl;
  cout << endl;
	  
  cout << "OPTIONAL ARGUMENTS:" << endl << endl;
  cout << "\t [--enzyme <string>] \n \t     The enzyme used to digest the proteins in the experiment. Default trypsin." << endl;
  cout << "\t [--decoy-prefix <string>] \n \t     Specifies the prefix of the protein names that indicates a decoy. Default decoy_" << endl;
  cout << "\t [--optimization <string>] \n \t     Specifies whether to do optimization at the protein, peptide or psm level Default protein" << endl;
  cout << "\t [--separate-searches <string>] \n \t     If the target and decoy searches were run separately, the option then allows the user to specify the location of the decoy search results, the target database search should be provided as required argument." << endl;
  cout << "\t [--fileroot <string>] \n \t     The fileroot string will be added as a prefix to all output file names. Default = none." <<endl;
  cout << "\t [--output-dir <directory>] \n \t     The name of the directory where output files will be created. Default = crux-output." << endl;
  cout << "\t [--overwrite <T/F>] \n \t     Replace existing files (T) or exit if attempting to overwrite (F). Default=F." << endl;
  cout << "\t [--skip-cleanup <T/F>] \n \t     When set to T, prevents the deletion of lookup tables created during the preprocessing step. Default = F." << endl; 
  cout << "\t [--re-run <directory>] \n \t      Re-run Barista analysis using a previously computed set of lookup tables." <<endl;  
  cout << "\t [--use-spec-features <T/F>] \n \t      When set to F, use minimal feature set. Default T." <<endl;  
  cout<<  "\t [--list-of-files <T/F>] \n \t      When set to T,there is a list of files for search result. Default= F." <<endl;  
  cout << endl; 

}

int Barista :: crux_set_command_line_options(int argc, char *argv[])
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
    "optimization",
    "spectrum-parser"
  };
  int num_options = sizeof(option_list)/sizeof(char*);

  const char* argument_list[] = {
    "database",
    "spectra",
    "search results"
  };
  int num_arguments = sizeof(argument_list)/sizeof(char*);

  initialize(
    argument_list, 
    num_arguments, 
    option_list,
    num_options,
    argc,
    argv
  );
  
  string db_source;
  string sqt_source;
  string ms2_source;
  string sqt_decoy_source;
  int separate_search_flag;
  string output_directory;
  
  string enzyme;
  string decoy_prefix;

  string dir_with_tables;
  int found_dir_with_tables;
  bool list_of_files_flag; 
  bool spec_features_flag;

  opt_type = get_string_parameter_pointer("optimization");
  

  fileroot = get_string_parameter_pointer("fileroot");
  if(fileroot != "__NULL_STR")
    fileroot.append(".");
  else
    fileroot = "";
  if(opt_type.compare("psm") == 0)
    qr.set_fileroot(fileroot);
  else if(opt_type.compare("peptide") == 0)
    pr.set_fileroot(fileroot);
  
  overwrite_flag = get_boolean_parameter("overwrite");
  if(opt_type.compare("psm") == 0)
    qr.set_overwrite_flag(overwrite_flag);
  if(opt_type.compare("peptide") == 0)
    pr.set_overwrite_flag(overwrite_flag);

  //options for the parser
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
  feature_file_name << output_directory << "/" << fileroot << "barista.features.txt";
  
  if(found_dir_with_tables){
      parser= new SQTParser();
      //set input and output dirs
      parser->set_output_dir(dir_with_tables);
      if(opt_type.compare("psm") == 0)
	{
	  qr.set_input_dir(dir_with_tables);
	  qr.set_output_dir(output_directory);
	}
      else if(opt_type.compare("peptide") == 0)
	{
	  pr.set_input_dir(dir_with_tables);
	  pr.set_output_dir(output_directory);
	}
      set_input_dir(dir_with_tables);
      set_output_dir(output_directory);

      carp(CARP_INFO, "directory with tables: %s", dir_with_tables.c_str());
      carp(CARP_INFO, "output_directory: %s", output_directory.c_str());
    }else{
      db_source = get_string_parameter_pointer("database");
      ms2_source = get_string_parameter_pointer("spectra");
      sqt_source = get_string_parameter_pointer("search results");
      sqt_decoy_source = get_string_parameter_pointer("separate-searches"); 
      list_of_files_flag=get_boolean_parameter("list-of-files");
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
        FILE_FORMAT_T format=check_file_format(files[i]);
        //sqt_source 
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

        if(sqt_decoy_source != "__NULL_STR")
	  separate_search_flag = 1;
        else
	  separate_search_flag = 0;
      
      parser->write_features_header();
      if(opt_type.compare("protein") == 0)
	parser->set_use_quadratic_features(1);
      if(parser->get_use_quadratic_features())
        parser->add_quadratic_features_header(); 
      d.get_features_header(parser->get_final_features_header());
        //set the output directory for the parser
        if(!parser->set_output_dir(output_directory))
	  return 0;
        //set input and output for the leaning algo (in and out are the same as the out for the parser)
		
	if(opt_type.compare("psm") == 0)
	  {
	    qr.set_input_dir(output_directory);
	    qr.set_output_dir(output_directory);
	  }
	else if(opt_type.compare("peptide") == 0)
	  {
	    pr.set_input_dir(output_directory);
	    pr.set_output_dir(output_directory);
	  }
	set_input_dir(output_directory);
	set_output_dir(output_directory);
	            
        if(!parser->set_database_source(db_source))
	  carp(CARP_FATAL, "could not find the database");
      
        if(separate_search_flag){
	    if(!parser->set_input_sources(ms2_source, sqt_source, sqt_decoy_source))
	      carp(CARP_FATAL, "could not extract features for training");
	    parser->set_num_hits_per_spectrum(1);
        }else{
	   if(!parser->set_input_sources(ms2_source, sqt_source))
	      carp(CARP_FATAL, "could not extract features for training");
         }
        //print some info
        carp(CARP_INFO, "database source: %s", db_source.c_str());
        carp(CARP_INFO, "sqt source: %s", sqt_source.c_str()); 
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
    carp(CARP_FATAL, "Please re-run with database, ms2 input and sqt input.");

  return 1;

}

int Barista::main(int argc, char **argv) {
 //int main(int argc, char **argv){

  if(!crux_set_command_line_options(argc,argv))
    return 1;

  if(opt_type.compare("psm") == 0)
    qr.run();
  else if(opt_type.compare("peptide") == 0)
    pr.run();
  else
    {
      //run();
      run_tries();
      //run_tries_multi_task();
    }
  if(skip_cleanup_flag != 1)
    parser->clean_up(out_dir);
  
  return 0;
}   

bool Barista :: needsOutputDirectory()
{
  return true;
}


string Barista::getName() {
  return "barista";
}

string Barista::getDescription() {
  return "Rank PSMs, peptides and proteins, assigning a confidence measure to "
         "each identification.";
}




COMMAND_T Barista::getCommand(){
  return BARISTA_COMMAND;
}

// returns file extension 
string Barista:: file_extension (string str){ 
  return str.substr(str.rfind('.')+1);
}


//determins which parser to be called 
FILE_FORMAT_T Barista:: check_file_format(string source){
  string ext = file_extension(source);
    if (ext=="sqt")
      return SQT_FORMAT;
    else if(ext=="txt")
      return DELIMITED_FORMAT;
    else 
      carp(CARP_DEBUG,"The file format is invalid");
    return INVALID_FORMAT;
}

void Barista :: get_protein_id(int pepind, vector<string> &prot){
  int * protinds= d.pepind2protinds(pepind) ;
  unsigned prot_length=d.pepind2num_prot(pepind);  
  for (unsigned k=0;k<prot_length;k++){
    string protein_str= d.ind2prot( protinds[k]);
    prot.push_back(protein_str); 
  }
}

void Barista :: print_protein_ids(vector<string> &prots, ofstream &os, int psmind){
  if(file_format_=="txt"){
    for (unsigned int j=0;j<prots.size();j++){
       if(j==prots.size()-1)
      	 os<<prots[j]<<"("<<d.psmind2peptide_position(psmind)<<")\t";
       else
         os<<prots[j]<<"("<<d.psmind2peptide_position(psmind)<<")"<<",";
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
*Local Variables: 
*mode: c
*c-basic-offset: 2
*End:
*/
