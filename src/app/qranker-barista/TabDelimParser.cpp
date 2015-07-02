#include "TabDelimParser.h"

/******************************/

TabDelimParser :: TabDelimParser() 
  : 	num_mixed_labels(0),     
	psmind_to_scan(0),
	psmind_to_charge(0),
	psmind_to_label(0), 
	psmind_to_num_pep(0),
	psmind_to_ofst(0),
	psmind_to_pepind(0),
	psmind_to_neutral_mass(0),
	psmind_to_peptide_mass(0),
	num_features(0),
	num_psm(0),
	num_pos_psm(0),
	num_neg_psm(0),
	num_pep(0),
	num_pep_in_all_psms(0),
	curr_ofst(0),
	psmind(0),
	x(0),
	num_xlink_features(0)
{
  
  //num_psm_features
  num_features = 17;
  //num_xlink_features
  num_xlink_features = 17;

  
  //final_hits_per_spectrum
  fhps = 3;
  //decoy prefix
  decoy_prefix = "random_";
  //max peptide length to be considered
  max_len = 50;
  //min peptide length to be considered
   min_len = 7;
  
}

void TabDelimParser :: clear()
{
  delete[] x; x = (double*)0;
  delete[] psmind_to_scan; psmind_to_scan = (int*)0;
  delete[] psmind_to_charge; psmind_to_charge = (int*)0;
  delete[] psmind_to_label; psmind_to_label = (int*)0;
  delete[] psmind_to_num_pep; psmind_to_num_pep = (int*)0;
  delete[] psmind_to_ofst; psmind_to_ofst = (int*)0;
  delete[] psmind_to_pepind; psmind_to_pepind = (int*)0;
  delete[] psmind_to_neutral_mass; psmind_to_neutral_mass = (double*)0;
  delete[] psmind_to_peptide_mass; psmind_to_peptide_mass = (double*)0;
}



TabDelimParser :: ~TabDelimParser()
{
  clear();
}

/* tokens are
0. scan
1. charge
2. spectrum precursor mz
3. spectrum neutral mass   
4. peptide mass    
5. delta_cn        
6. sp score        
7. xcorr score     
8. xcorr rank      
9. matches/spectrum        
10. sequence        
11. cleavage type   
12. protein id      
13. flanking aa     
14. nzstates        
15. rtime max diff  
16. peptides/spectrum       
17. xcorr sum diff  
18. xcorr max diff
*/



void TabDelimParser :: get_tokens(string &line, vector<string>&tokens, string &delim)
{
  tokens.erase(tokens.begin(), tokens.end());
  string tmp = line;
  string tok;
  size_t pos = tmp.find(delim);
  while(pos != string::npos)
    {
      tok = tmp.substr(0,pos);
      tokens.push_back(tok);
      tmp = tmp.substr(pos+1, tmp.size());
      pos = tmp.find(delim);
    }
  //last token
  tokens.push_back(tmp);
}

void TabDelimParser :: first_pass(ifstream &fin)
{
  string line;
  vector<string> tokens;
  vector<string> peptides;
  string delim1 = "\t";
  string delim2 = ",";
  while(!fin.eof())
    {
      getline(fin,line);
      get_tokens(line, tokens, delim1);
  
      if(tokens.size() > 1)
	{
	  //get all peptides and fill pep_to_ind tables
	  string peps = tokens[10];
	  get_tokens(peps,peptides,delim2);
	  
	  for(unsigned int i = 0; i < peptides.size(); i++)
	    {
	      string pep = peptides[i];
	      int pep_ind = -1;
	      //add peptide to pep_to_ind and ind_to_pep maps
	      if(pep_to_ind.find(pep) == pep_to_ind.end())
		{
		  pep_ind = num_pep;
		  //add a new peptide
		  pep_to_ind[pep] = pep_ind;
		  ind_to_pep[pep_ind] = pep;
		  num_pep++;
	  	}
	      else
		{
		  pep_ind = pep_to_ind[pep];
		  string p = ind_to_pep[pep_ind];
		  if(pep.compare(p) != 0)
		    cout << "warning : did not find peptide in ind_to_pep_table\n"; 
		}
	    }
	  num_pep_in_all_psms += peptides.size();
	  num_psm++;
	}
      
    }
}


void TabDelimParser :: allocate_feature_space()
{
  //space for feature vector
  x = new double[num_features];
  memset(x,0,sizeof(double)*num_features);

  //space for psminfo
  psmind_to_pepind = new int[num_pep_in_all_psms];
  memset(psmind_to_pepind,0,sizeof(int)*num_pep_in_all_psms);
  psmind_to_num_pep = new int[num_psm];
  memset(psmind_to_num_pep,0,sizeof(int)*num_psm);
  psmind_to_ofst = new int[num_psm];
  memset(psmind_to_ofst,0,sizeof(int)*num_psm);
  psmind_to_charge = new int[num_pep_in_all_psms];
  memset(psmind_to_charge,0,sizeof(int)*num_pep_in_all_psms);

  psmind_to_neutral_mass = new double[num_pep_in_all_psms];
  memset(psmind_to_neutral_mass,0,sizeof(double)*num_pep_in_all_psms);
  psmind_to_peptide_mass = new double[num_pep_in_all_psms];
  memset(psmind_to_peptide_mass,0,sizeof(double)*num_pep_in_all_psms);

  psmind_to_scan = new int[num_psm];
  memset(psmind_to_scan,0,sizeof(int)*num_psm);
  psmind_to_label = new int[num_psm];
  memset(psmind_to_label,0,sizeof(int)*num_psm);

}


void TabDelimParser :: extract_psm_features(vector<string> & tokens, double *x)
{
  memset(x,0,sizeof(double)*num_features);
  //log rank by Sp

  //deltaLCN

  //deltaCN

  //xcorr score
  x[3] = atof(tokens[7].c_str());
  //sp score
  
  //matched ions/predicted ions
  
  //observed mass
  
  //peptide length
  
  //charge
  
  //whether n-terminus and c-terminus have proper cleavage sites
  
  // missed cleavages
  
  // number of sequence_comparisons

  //difference between measured and calculated mass

  // absolute value of difference between measured and calculated mass

}



void TabDelimParser :: second_pass(ifstream &fin, int label)
{
  string line;
  vector<string> tokens;
  
  vector<string> peptides;
  vector<string> charge_str;
  vector<string> spectrum_neutral_mass_str;
  vector<string> peptide_mass_str;
  string delim1 = "\t";
  string delim2 = ",";
  while(!fin.eof())
    {
      getline(fin,line);
      get_tokens(line, tokens, delim1);
  
      if(tokens.size() > 1)
	{
	  //extract features
	  extract_psm_features(tokens, x);
	  f_psm.write((char*)x, sizeof(double)*num_features);

	  //fill in tables

	  //get the scan
	  int scan = atoi(tokens[0].c_str());
	  psmind_to_scan[psmind] = scan;

	  //get charges
	  string ch = tokens[1];
	  get_tokens(ch, charge_str,delim2);
	  for(unsigned int i = 0; i < charge_str.size(); i++)
	    psmind_to_charge[curr_ofst+i] = atoi(charge_str[i].c_str());
	  
	  //get spectrum neutral mass
	  string neut_mass = tokens[3];
	  get_tokens(neut_mass, spectrum_neutral_mass_str,delim2);
	  for(unsigned int i = 0; i < spectrum_neutral_mass_str.size(); i++)
	    psmind_to_neutral_mass[curr_ofst+i] = atof(spectrum_neutral_mass_str[i].c_str());
	  
	  //get peptide mass
	  string pep_mass = tokens[4];
	  get_tokens(pep_mass, peptide_mass_str,delim2);
	  vector<double> peptide_mass;
	  peptide_mass.resize(peptide_mass_str.size(),0);
	  for(unsigned int i = 0; i < peptide_mass_str.size(); i++)
	    psmind_to_peptide_mass[curr_ofst+i] = atof(peptide_mass_str[i].c_str());
		  
	  //get all peptides
	  string peps = tokens[10];
	  get_tokens(peps,peptides,delim2);
	  for(unsigned int i = 0; i < peptides.size(); i++)
	    {
	      string pep = peptides[i];
	      int pep_ind = -1;
	      //add peptide to pep_to_ind and ind_to_pep maps
	      if(pep_to_ind.find(pep) == pep_to_ind.end())
		cout << "warning : did not find peptide in ind_to_pep_table\n"; 
	      else
		pep_ind = pep_to_ind[pep];
	      psmind_to_pepind[curr_ofst+i] = pep_ind;
	    }
	  psmind_to_num_pep[psmind] = peptides.size();
	  psmind_to_ofst[psmind] = curr_ofst;
	  
	  psmind_to_label[psmind] = label;
	  if(label == 1)
	    num_pos_psm++;
	  else
	    num_neg_psm++;

	  //augment counters
	  curr_ofst+= peptides.size();
	  psmind++;  
	}
    }
}

/*****************************************************************************************/

void TabDelimParser :: save_data_in_binary(string out_dir)
{
  ostringstream fname;
  //write out data summary
  fname << out_dir << "/summary.txt";
  ofstream f_summary(fname.str().c_str());
  //psm info
  f_summary << num_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  //peptide info
  f_summary << num_pep << endl;
  f_summary.close();
  fname.str("");
  
  //psmind_to_pepind
  fname << out_dir << "/psmind_to_pepind.txt";
  ofstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  f_psmind_to_pepind.write((char*)psmind_to_pepind,sizeof(int)*num_pep_in_all_psms);
  f_psmind_to_pepind.close();
  fname.str("");

  //psmind_to_num_pep
  fname << out_dir << "/psmind_to_num_pep.txt";
  ofstream f_psmind_to_num_pep(fname.str().c_str(),ios::binary);
  f_psmind_to_num_pep.write((char*)psmind_to_num_pep,sizeof(int)*num_psm);
  f_psmind_to_num_pep.close();
  fname.str("");

  //psmind_to_ofst
  fname << out_dir << "/psmind_to_ofst.txt";
  ofstream f_psmind_to_ofst(fname.str().c_str(),ios::binary);
  f_psmind_to_ofst.write((char*)psmind_to_ofst,sizeof(int)*num_psm);
  f_psmind_to_ofst.close();
  fname.str("");

  //psmind_to_scan
  fname << out_dir << "/psmind_to_scan.txt";
  ofstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  f_psmind_to_scan.write((char*)psmind_to_scan,sizeof(int)*num_psm);
  f_psmind_to_scan.close();
  fname.str("");

  //psmind_to_charge
  fname << out_dir << "/psmind_to_charge.txt";
  ofstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  f_psmind_to_charge.write((char*)psmind_to_charge,sizeof(int)*num_pep_in_all_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  ofstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  f_psmind_to_label.write((char*)psmind_to_label,sizeof(int)*num_psm);
  f_psmind_to_label.close();
  fname.str("");
  
  //ind_to_pep
  fname << out_dir << "/ind_to_pep.txt";
  ofstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = ind_to_pep.begin(); it != ind_to_pep.end(); it++)
    f_ind_to_pep << it->first << " " << it->second << "\n";
  f_ind_to_pep.close();
  fname.str("");

  //pep_to_ind
  fname << out_dir << "/pep_to_ind.txt";
  ofstream f_pep_to_ind(fname.str().c_str(),ios::binary);
  for(map<string,int>::iterator it = pep_to_ind.begin(); it != pep_to_ind.end(); it++)
    f_pep_to_ind << it->first << " " << it->second << "\n";
  f_pep_to_ind.close();
  fname.str("");

}

void TabDelimParser :: clean_up(string dir)
{

  ostringstream fname;
      
  //fname << out_dir << "/summary.txt";
  //ofstream f_summary(fname.str().c_str());

  fname << dir << "/psm.txt";
  remove(fname.str().c_str());
  fname.str("");
  
  //psmind_to_pepind
  fname << out_dir << "/psmind_to_pepind.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_num_pep
  fname << out_dir << "/psmind_to_num_pep.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_num_pep
  fname << out_dir << "/psmind_to_ofst.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_scan
  fname << out_dir << "/psmind_to_scan.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_charge
  fname << out_dir << "/psmind_to_charge.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  remove(fname.str().c_str());
  fname.str("");

  //ind_to_pep
  fname << out_dir << "/ind_to_pep.txt";
  remove(fname.str().c_str());
  fname.str("");

  //pep_to_ind
  fname << out_dir << "/pep_to_ind.txt";
  remove(fname.str().c_str());
  fname.str("");

}



/******************************************************************************************/


int TabDelimParser :: run(vector<string> &filenames)
{
  string line;
  for(unsigned int i = 0; i < filenames.size(); i++)
    {
      string fname = filenames[i];
      cout << fname << endl;
      ifstream fin(fname.c_str());
      if(!fin.is_open())
	{
	  cout << "could not open " << fname << " for reading" << endl;
	  return 0;
	}
      getline(fin,line);
      first_pass(fin);
      fin.close();
    }
  cout << num_psm << " " << num_pep << " " << num_pep_in_all_psms << endl;
  allocate_feature_space();

  ostringstream fname;
  fname << out_dir << "/psm.txt";
  f_psm.open(fname.str().c_str());
  fname.str("");

  int label = 0;
  for(unsigned int i = 0; i < filenames.size(); i++)
    {
      string fname = filenames[i];
      cout << fname << endl;
      ifstream fin(fname.c_str());
      if(!fin.is_open())
	{
	  cout << "could not open " << fname << " for reading" << endl;
	  return 0;
	}
      getline(fin,line);
      if(i == 0)
	label = 1;
      else
	label = -1;
      second_pass(fin,label);
      fin.close();
    }
  f_psm.close();
  cout << psmind << " " << num_pos_psm << " " << num_neg_psm << endl;
  save_data_in_binary(out_dir);
  return 1;
}

/*******************************************xlinking functions*************************************/


/* xlink tokens
 * 0. scan    
 * 1. charge  
 * 2. spectrum precursor m/z  
 * 3. spectrum neutral mass   
 * 4. peptide mass mono       
 * 5. peptide mass average    
 * 6. mass error(ppm) 
 * 7. sp score        
 * 8. sp rank 
 * 9. b/y ions matched        
 * 10. b/y ions total  
 * 11. xcorr score     
 * 12. xcorr rank      
 * 13. p-value 
 * 14. matches/spectrum        
 * 15. sequence        
 * 16. protein id(loc) 1       
 * 17. protein id(loc) 2
 */


void TabDelimParser :: first_pass_xlink(ifstream &fin)
{
  string line;
  vector<string> tokens;
  vector<string> subtokens;
  vector<string> peptides;
  string delim1 = "\t";
  string delim2 = ",";
  string delim3 = " ";
  while(!fin.eof())
    {
      getline(fin,line);
      get_tokens(line, tokens, delim1);
      
      if(tokens.size() > 1)
	{
	  //get all peptides and fill tables
	  string pep_and_loc = tokens[15];
	  get_tokens(pep_and_loc,subtokens,delim3);
	  string peps = subtokens[0];
	  get_tokens(peps,peptides,delim2);
	  psmind_to_peptide1[num_psm] = peptides[0];
	  psmind_to_peptide2[num_psm] = peptides[1];
	  psmind_to_loc[num_psm] = subtokens[1];
	  //proteins
	  psmind_to_protein1[num_psm] = tokens[16];
	  psmind_to_protein2[num_psm] = tokens[17];
	  
	  num_psm++;
	}
    }
}


void TabDelimParser :: allocate_feature_space_xlink()
{
  //space for feature vector
  x = new double[num_xlink_features];
  memset(x,0,sizeof(double)*num_xlink_features);

  psmind_to_label = new int[num_psm];
  memset(psmind_to_label,0,sizeof(int)*num_psm);

}


void TabDelimParser :: extract_xlink_features(vector<string> & tokens, double *x)
{
  memset(x,0,sizeof(double)*num_xlink_features);
  //log rank by Sp
  x[0]=atof(tokens[8].c_str());
  //deltaLCN

  //deltaCN

  //xcorr score
  x[3] = atof(tokens[11].c_str());
  //sp score
  x[4] = atof(tokens[7].c_str());
  //matched ions/predicted ions
  
  //observed mass
  
  //peptide length
  
  //charge
  
  //whether n-terminus and c-terminus have proper cleavage sites
  
  // missed cleavages
  
  // number of sequence_comparisons

  //difference between measured and calculated mass

  // absolute value of difference between measured and calculated mass

}



void TabDelimParser :: second_pass_xlink(ifstream &fin, int label)
{
  string line;
  vector<string> tokens;
  
  vector<string> peptides;
  vector<string> charge_str;
  vector<string> spectrum_neutral_mass_str;
  vector<string> peptide_mass_str;
  string delim1 = "\t";
  string delim2 = ",";
  while(!fin.eof())
    {
      getline(fin,line);
      get_tokens(line, tokens, delim1);
  
      if(tokens.size() > 1)
	{
	  //extract features
	  extract_xlink_features(tokens, x);
	  f_psm.write((char*)x, sizeof(double)*num_features);
		  
	  psmind_to_label[psmind] = label;
	  if(label == 1)
	    num_pos_psm++;
	  else
	    num_neg_psm++;

	  //augment counters
	  psmind++;  
	}
    }
}


void TabDelimParser :: save_data_in_binary_xlink(string out_dir)
{
  ostringstream fname;
  //write out data summary
  fname << out_dir << "/summary.txt";
  ofstream f_summary(fname.str().c_str());
  //psm info
  f_summary << num_xlink_features << " " << num_psm << " " << num_pos_psm << " " << num_neg_psm << endl;
  f_summary.close();
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  ofstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  f_psmind_to_label.write((char*)psmind_to_label,sizeof(int)*num_psm);
  f_psmind_to_label.close();
  fname.str("");
  
  //psmind_to_peptide1
  fname << out_dir << "/psmind_to_peptide1.txt";
  ofstream f_psmind_to_peptide1(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_peptide1.begin(); it != psmind_to_peptide1.end(); it++)
    f_psmind_to_peptide1 << it->first << " " << it->second << "\n";
  f_psmind_to_peptide1.close();
  fname.str("");
  
  //psmind_to_peptide2
  fname << out_dir << "/psmind_to_peptide2.txt";
  ofstream f_psmind_to_peptide2(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_peptide2.begin(); it != psmind_to_peptide2.end(); it++)
    f_psmind_to_peptide2 << it->first << " " << it->second << "\n";
  f_psmind_to_peptide2.close();
  fname.str("");

  //psmind_to_loc
  fname << out_dir << "/psmind_to_loc.txt";
  ofstream f_psmind_to_loc(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_loc.begin(); it != psmind_to_loc.end(); it++)
    f_psmind_to_loc << it->first << " " << it->second << "\n";
  f_psmind_to_loc.close();
  fname.str("");

  //psmind_to_peptide1
  fname << out_dir << "/psmind_to_protein1.txt";
  ofstream f_psmind_to_protein1(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_protein1.begin(); it != psmind_to_protein1.end(); it++)
    f_psmind_to_protein1 << it->first << " " << it->second << "\n";
  f_psmind_to_protein1.close();
  fname.str("");
  
  //psmind_to_peptide2
  fname << out_dir << "/psmind_to_protein2.txt";
  ofstream f_psmind_to_protein2(fname.str().c_str(),ios::binary);
  for(map<int,string>::iterator it = psmind_to_protein2.begin(); it != psmind_to_protein2.end(); it++)
    f_psmind_to_protein2 << it->first << " " << it->second << "\n";
  f_psmind_to_protein2.close();
  fname.str("");

}


void TabDelimParser :: clean_up_xlink(string dir)
{

  ostringstream fname;
      
  //fname << out_dir << "/summary.txt";
  //ofstream f_summary(fname.str().c_str());

  fname << dir << "/psm.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_label
  fname << out_dir << "/psmind_to_label.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_peptide1
  fname << out_dir << "/psmind_to_peptide1.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_peptide2
  fname << out_dir << "/psmind_to_peptide2.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_loc
  fname << out_dir << "/psmind_to_loc.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_protein1
  fname << out_dir << "/psmind_to_protein1.txt";
  remove(fname.str().c_str());
  fname.str("");

  //psmind_to_protein2
  fname << out_dir << "/psmind_to_protein2.txt";
  remove(fname.str().c_str());
  fname.str("");

}



int TabDelimParser :: run_on_xlink(vector<string> &filenames)
{
  string line;
  for(unsigned int i = 0; i < filenames.size(); i++)
    {
      string fname = filenames[i];
      cout << fname << endl;
      ifstream fin(fname.c_str());
      if(!fin.is_open())
	{
	  cout << "could not open " << fname << " for reading" << endl;
	  return 0;
	}
      getline(fin,line);
      first_pass_xlink(fin);
      fin.close();
    
    }
  cout << num_psm  << endl;
  allocate_feature_space_xlink();
    
  ostringstream fname;
  fname << out_dir << "/psm.txt";
  f_psm.open(fname.str().c_str());
  fname.str("");

  int label = 0;
  for(unsigned int i = 0; i < filenames.size(); i++)
    {
      string fname = filenames[i];
      cout << fname << endl;
      ifstream fin(fname.c_str());
      if(!fin.is_open())
	{
	  cout << "could not open " << fname << " for reading" << endl;
	  return 0;
	}
      getline(fin,line);
      if(i == 0)
	label = 1;
      else
	label = -1;
      second_pass_xlink(fin,label);
      fin.close();
    }
  f_psm.close();
  cout << psmind << " " << num_pos_psm << " " << num_neg_psm << endl;
  save_data_in_binary_xlink(out_dir);
  return 1;
}


