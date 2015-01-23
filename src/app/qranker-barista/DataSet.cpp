#include "DataSet.h"
#include "SQTParser.h"

Dataset::Dataset() 
  : num_psms(0), num_pos_psms(0), num_neg_psms(0), 
    psmind_to_features((double*)0), 
    psmind_to_label((int*)0),
    psmind_to_pepind((int*)0),
    psmind_to_scan(0),
    psmind_to_charge(0),
    psmind_to_precursor_mass(0),
    psmind_to_fileind(0),
    num_pep(0), num_pos_pep(0), num_neg_pep(0),
    pepind_to_label(0),
    num_prot(0), num_pos_prot(0), num_neg_prot(0),
    protind_to_label(0),
    protind_to_num_all_pep(0),
    psmind_to_sp_rank((int*)0),
    protind_to_length(0),
    psmind_to_xcorr_rank((int*)0),
    psmind_to_matches_spectrum((int*)0),
    psmind_to_by_ions_matched((double*)0),
    psmind_to_by_ions_total((double*)0),
    psmind_to_peptide_position((int*)0)
{
}

Dataset::~Dataset()
{
  delete[] psmind_to_features; psmind_to_features = (double*)0;
  delete[] psmind_to_label; psmind_to_label = (int*)0;
  delete[] psmind_to_pepind; psmind_to_pepind = (int*)0;
  delete[] psmind_to_scan; psmind_to_scan = (int*)0;
  delete[] psmind_to_charge; psmind_to_charge = (int*)0;
  delete[] psmind_to_precursor_mass; psmind_to_precursor_mass = (double*)0;
  delete[] psmind_to_sp_rank; psmind_to_sp_rank=(int*)0;//Sp Rank
  delete[] psmind_to_xcorr_rank; psmind_to_xcorr_rank=(int*)0;//xcorr Rank
  delete[] psmind_to_matches_spectrum; psmind_to_matches_spectrum=(int*)0;//matches/spectrum 
  delete[] psmind_to_by_ions_matched; psmind_to_by_ions_matched=(double*)0;// b/y ions matched
  delete[] psmind_to_by_ions_total; psmind_to_by_ions_total=(double*)0;// b/y ions total
  delete[] psmind_to_peptide_position; psmind_to_peptide_position=(int*)0;//peptide position in protein 
  fileind_to_fname.clear();

  ind_to_pep.clear();
  pepind_to_psminds.clear();
  pepind_to_protinds.clear();
  delete [] pepind_to_label; pepind_to_label = (int*)0;

  delete [] psmind_to_fileind; psmind_to_fileind = (int*)0;
  delete[] protind_to_label; protind_to_label = (int*)0;
  delete[] protind_to_num_all_pep; protind_to_num_all_pep = (int*)0;
  delete[] protind_to_length; protind_to_length = (int*)0;
  protind_to_pepinds.clear();
  ind_to_prot.clear();

}


/****************************************************************************/

void Dataset :: load_data_psm_training()
{

  ostringstream fname;
  fname << in_dir << "/summary";
  ifstream f_summary(fname.str().c_str());
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary.close();
  fname.str("");

  //psm features
  fname << in_dir << "/" << "psm";
  ifstream f_psm_feat(fname.str().c_str(),ios::binary);
  if(!f_psm_feat.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_features = new double[num_psms*num_features];
  f_psm_feat.read((char*)psmind_to_features,sizeof(double)*num_psms*num_features);
  f_psm_feat.close();
  fname.str("");
}

void Dataset :: clear_data_psm_training()
{
  delete [] psmind_to_features; psmind_to_features = (double*)0;
}

void Dataset :: load_labels_psm_training()
{
  ostringstream fname;
  fname << in_dir << "/summary";
  ifstream f_summary(fname.str().c_str());
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary.close();
  fname.str("");

  //psmind_to_label
  fname << in_dir << "/psmind_to_label";
  ifstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  psmind_to_label = new int[num_psms];
  f_psmind_to_label.read((char*)psmind_to_label,sizeof(int)*num_psms);
  f_psmind_to_label.close();
  fname.str("");
}

void Dataset :: clear_labels_psm_training()
{
  delete [] psmind_to_label; psmind_to_label = (int*)0;
}

void Dataset :: load_data_psm_results()
{

  ostringstream fname;
  fname << in_dir << "/summary";
  ifstream f_summary(fname.str().c_str());
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary.close();
  fname.str("");

  //psmind_to_pepind
  fname << in_dir << "/psmind_to_pepind";
  ifstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_pepind.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_pepind = new int[num_psms];
  f_psmind_to_pepind.read((char*)psmind_to_pepind,sizeof(int)*num_psms);
  f_psmind_to_pepind.close();
  fname.str("");
  
  //psmind_to_scan
  fname << in_dir << "/psmind_to_scan";
  ifstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_scan.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_scan = new int[num_psms];
  f_psmind_to_scan.read((char*)psmind_to_scan,sizeof(int)*num_psms);
  f_psmind_to_scan.close();
  fname.str("");

  //psmind_to_charge
  fname << in_dir << "/psmind_to_charge";
  ifstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_charge.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_charge = new int[num_psms];
  f_psmind_to_charge.read((char*)psmind_to_charge,sizeof(int)*num_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //psmind_to_xcorr
  fname << in_dir << "/psmind_to_xcorr";
  ifstream f_psmind_to_xcorr(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_xcorr.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_xcorr = new double[num_psms];
  f_psmind_to_xcorr.read((char*)psmind_to_xcorr,sizeof(double)*num_psms);
  f_psmind_to_xcorr.close();
  fname.str("");

  //psmind_to_deltaCn
  fname << in_dir << "/psmind_to_deltaCn";
  ifstream f_psmind_to_deltaCn(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_deltaCn.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_deltaCn = new double[num_psms];
  f_psmind_to_deltaCn.read((char*)psmind_to_deltaCn,sizeof(double)*num_psms);
  f_psmind_to_deltaCn.close();
  fname.str("");
  
  //psmind_to_spscore
  fname << in_dir << "/psmind_to_spscore";
  ifstream f_psmind_to_spscore(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_spscore.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_spscore = new double[num_psms];
  f_psmind_to_spscore.read((char*)psmind_to_spscore,sizeof(double)*num_psms);
  f_psmind_to_spscore.close();
  fname.str("");

  //psmind_to_calculated_mass
  fname << in_dir << "/psmind_to_calculated_mass";
  ifstream f_psmind_to_calculated_mass(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_calculated_mass.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_calculated_mass = new double[num_psms];
  f_psmind_to_calculated_mass.read((char*)psmind_to_calculated_mass,sizeof(double)*num_psms);
  f_psmind_to_calculated_mass.close();
  fname.str("");

  //psmind_to_precursor_mass
  fname << in_dir << "/psmind_to_precursor_mass";
  ifstream f_psmind_to_precursor_mass(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_precursor_mass.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_precursor_mass = new double[num_psms];
  f_psmind_to_precursor_mass.read((char*)psmind_to_precursor_mass,sizeof(double)*num_psms);
  f_psmind_to_precursor_mass.close();
  fname.str("");
  
  //fileind_to_fname
  fname << in_dir << "/fileind_to_fname";
  ifstream f_fileind_to_fname(fname.str().c_str(),ios::binary);
  if(!f_fileind_to_fname.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  int ind;
  string filename;
  f_fileind_to_fname >> ind;
  f_fileind_to_fname >> filename;
  while(!f_fileind_to_fname.eof())
    {
      fileind_to_fname[ind] = filename;
      f_fileind_to_fname >> ind;
      f_fileind_to_fname >> filename;
    }
  f_fileind_to_fname.close();
  fname.str("");
  
  //psmind_to_filename
  fname << in_dir << "/psmind_to_fileind";
  ifstream f_psmind_to_fileind(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_fileind.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_fileind = new int[num_psms];
  f_psmind_to_fileind.read((char*)psmind_to_fileind,sizeof(int)*num_psms);
  f_psmind_to_fileind.close();
  fname.str("");

  //ind_to_pep
  fname << in_dir << "/ind_to_pep";
  ifstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  string pep;
  f_ind_to_pep >> ind;
  f_ind_to_pep >> pep;
  while(!f_ind_to_pep.eof())
    {
      ind_to_pep[ind] = pep;
      f_ind_to_pep >> ind;
      f_ind_to_pep >> pep;
    }
  f_ind_to_pep.close();
  fname.str("");

  //ind_to_prot
  fname << in_dir << "/ind_to_prot";
  ifstream f_ind_to_prot(fname.str().c_str(),ios::binary);
 
  string prot;
  f_ind_to_prot >> ind;
  f_ind_to_prot >> prot;
  while(!f_ind_to_prot.eof())
    {
      ind_to_prot[ind] = prot;
      f_ind_to_prot >> ind;
      f_ind_to_prot >> prot;
    }
  f_ind_to_prot.close();
  fname.str("");

  //pepind_to_protinds
  fname << in_dir << "/pepind_to_protinds";
  ifstream f_pepind_to_protinds(fname.str().c_str(),ios::binary);
  if(!f_pepind_to_protinds.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  pepind_to_protinds.load(f_pepind_to_protinds);
  f_pepind_to_protinds.close();
  fname.str("");
  
  //psmind_to_Sp_rank 
  fname << in_dir << "/psmind_to_sp_rank";
  ifstream f_psmind_to_sp_rank(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_sp_rank.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }

  psmind_to_sp_rank = new int[num_psms];
  f_psmind_to_sp_rank.read((char*)psmind_to_sp_rank,sizeof(int)*num_psms);
  f_psmind_to_sp_rank.close();
  fname.str("");

  //psmind_to_xcorr_rank 
  fname << in_dir << "/psmind_to_xcorr_rank";
  ifstream f_psmind_to_xcorr_rank(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_xcorr_rank.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_xcorr_rank = new int[num_psms];
  f_psmind_to_xcorr_rank.read((char*)psmind_to_xcorr_rank,sizeof(int)*num_psms);
  f_psmind_to_xcorr_rank.close();
  fname.str("");

  //psmind_to_by_ions_matched 
  fname << in_dir << "/psmind_to_by_ions_matched";
  ifstream f_psmind_to_by_ions_matched(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_by_ions_matched.is_open()){ 
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_by_ions_matched = new double[num_psms];
  f_psmind_to_by_ions_matched.read((char*)psmind_to_by_ions_matched,sizeof(double)*num_psms);
  f_psmind_to_by_ions_matched.close();
  fname.str("");

  //psmind_to_by_ions_total 
  fname << in_dir << "/psmind_to_by_ions_total";
  ifstream f_psmind_to_by_ions_total(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_by_ions_total.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_by_ions_total= new double[num_psms];
  f_psmind_to_by_ions_total.read((char*)psmind_to_by_ions_total,sizeof(double)*num_psms);
  f_psmind_to_by_ions_total.close();
  fname.str("");
  
  //psmind_to_matches_spectrum 
  fname << in_dir << "/psmind_to_matches_spectrum";
  ifstream f_psmind_to_matches_spectrum(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_matches_spectrum.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_matches_spectrum = new int[num_psms];
  f_psmind_to_matches_spectrum.read((char*)psmind_to_matches_spectrum,sizeof(int)*num_psms);
  f_psmind_to_matches_spectrum.close();
  fname.str("");
  
  //psmind_to_peptide_position 
  fname << in_dir << "/psmind_to_peptide_position";
  ifstream f_psmind_to_peptide_position(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_peptide_position.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_peptide_position= new int[num_psms];
  f_psmind_to_peptide_position.read((char*)psmind_to_peptide_position,sizeof(int)*num_psms);
  f_psmind_to_peptide_position.close();
  fname.str("");

}

void Dataset :: clear_data_psm_results()
{
  delete [] psmind_to_label; psmind_to_label = (int*)0; 
  delete [] psmind_to_pepind; psmind_to_pepind = (int*)0;
  delete [] psmind_to_scan; psmind_to_scan = (int*)0;
  delete [] psmind_to_charge; psmind_to_charge = (int*)0;
  delete [] psmind_to_xcorr; psmind_to_xcorr = (double*)0;
  delete [] psmind_to_precursor_mass; psmind_to_precursor_mass = (double*)0;
  fileind_to_fname.clear();
  delete [] psmind_to_fileind; psmind_to_fileind = (int*)0;
  delete [] psmind_to_sp_rank; psmind_to_sp_rank=(int*)0;//sp rank
  delete [] psmind_to_xcorr_rank; psmind_to_xcorr_rank=(int*)0;//xcorr rank
  delete[] psmind_to_matches_spectrum; psmind_to_matches_spectrum=(int*)0;//matches/spectrum
  delete [] psmind_to_by_ions_total; psmind_to_by_ions_total=(double*)0; //b/y ions total 
  delete [] psmind_to_by_ions_matched;  psmind_to_by_ions_matched=(double*)0; //by ions matched 
  delete [] psmind_to_peptide_position; psmind_to_peptide_position=(int*)0; //position of peptide in 
  //the begining of protein 
  ind_to_pep.clear();
  ind_to_prot.clear();
  pepind_to_protinds.clear();
}

string& Dataset :: psmind2fname(int psmind)
{
  int fileind = psmind_to_fileind[psmind];
  return fileind_to_fname[fileind];
}


/************************************************************/
void Dataset :: load_data_prot_training()
{

  ostringstream fname;
  fname << in_dir << "/summary";
  ifstream f_summary(fname.str().c_str());
  if(!f_summary.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary >> num_pep;
  f_summary >> num_pos_pep;
  f_summary >> num_neg_pep;
  f_summary >> num_prot;
  f_summary >> num_pos_prot;
  f_summary >> num_neg_prot;
  f_summary.close();
  fname.str("");
  
  //psm features
  fname << in_dir << "/" << "psm";
  ifstream f_psm_feat(fname.str().c_str(),ios::binary);
  if(!f_psm_feat.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_features = new double[num_psms*num_features];
  f_psm_feat.read((char*)psmind_to_features,sizeof(double)*num_psms*num_features);
  f_psm_feat.close();
  fname.str("");


  //pepind_to_psminds
  fname << in_dir << "/pepind_to_psminds";
  ifstream f_pepind_to_psminds(fname.str().c_str(),ios::binary);
  if(!f_pepind_to_psminds.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  pepind_to_psminds.load(f_pepind_to_psminds);
  f_pepind_to_psminds.close();
  fname.str("");
  
  //protind_to_num_all_pep
  fname << in_dir << "/protind_to_num_all_pep";
  ifstream f_protind_to_num_all_pep(fname.str().c_str(),ios::binary);
  if(!f_protind_to_num_all_pep.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  protind_to_num_all_pep = new int[num_prot];
  f_protind_to_num_all_pep.read((char*)protind_to_num_all_pep,sizeof(int)*num_prot);
  f_protind_to_num_all_pep.close();
  fname.str("");

  //protind_to_pepinds
  fname << in_dir << "/protind_to_pepinds";
  ifstream f_protind_to_pepinds(fname.str().c_str(),ios::binary);
  if(!f_protind_to_pepinds.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  protind_to_pepinds.load(f_protind_to_pepinds);
  f_protind_to_pepinds.close();
  fname.str("");

  //pepind_to_protinds
  fname << in_dir << "/pepind_to_protinds";
  ifstream f_pepind_to_protinds(fname.str().c_str(),ios::binary);
  if(!f_pepind_to_protinds.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  pepind_to_protinds.load(f_pepind_to_protinds);
  f_pepind_to_protinds.close();
  fname.str("");

}

void Dataset :: clear_data_prot_training()
{
  delete [] psmind_to_features; psmind_to_features = (double*)0;
  delete [] protind_to_num_all_pep; protind_to_num_all_pep = (int*)0;
}

void Dataset :: load_labels_prot_training()
{

  ostringstream fname;
  fname << in_dir << "/summary";
  ifstream f_summary(fname.str().c_str());
  if(!f_summary.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary >> num_pep;
  f_summary >> num_pos_pep;
  f_summary >> num_neg_pep;
  f_summary >> num_prot;
  f_summary >> num_pos_prot;
  f_summary >> num_neg_prot;
  f_summary.close();
  fname.str("");

  //psmind_to_label
  fname << in_dir << "/psmind_to_label";
  ifstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_label.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_label = new int[num_psms];
  f_psmind_to_label.read((char*)psmind_to_label,sizeof(int)*num_psms);
  f_psmind_to_label.close();
  fname.str("");
  
  //pepind_to_label
  fname << in_dir << "/pepind_to_label";
  ifstream f_pepind_to_label(fname.str().c_str(),ios::binary);
  if(!f_pepind_to_label.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  pepind_to_label = new int[num_pep];
  f_pepind_to_label.read((char*)pepind_to_label,sizeof(int)*num_pep);
  f_pepind_to_label.close();
  fname.str("");
  
  //protind_to_label
  fname << in_dir << "/protind_to_label";
  ifstream f_protind_to_label(fname.str().c_str(),ios::binary);
  if(!f_protind_to_label.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  protind_to_label = new int[num_prot];
  f_protind_to_label.read((char*)protind_to_label,sizeof(int)*num_prot);
  f_protind_to_label.close();
  fname.str("");

  //ind_to_pep
  fname << in_dir << "/ind_to_pep";
  ifstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  int ind;
  string pep;
  f_ind_to_pep >> ind;
  f_ind_to_pep >> pep;
  while(!f_ind_to_pep.eof()){
      ind_to_pep[ind] = pep;
      f_ind_to_pep >> ind;
      f_ind_to_pep >> pep;
  }
  f_ind_to_pep.close();
  fname.str("");


}

void Dataset :: clear_labels_prot_training()
{
  delete [] psmind_to_label; psmind_to_label = (int*)0;
  delete [] pepind_to_label; pepind_to_label = (int*)0;
  delete [] protind_to_label; protind_to_label = (int*)0;
  ind_to_pep.clear();
}


void Dataset :: load_data_all_results(){

  ostringstream fname;

  //ind_to_pep
  fname << in_dir << "/ind_to_pep";
  ifstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  int ind;
  string pep;
  f_ind_to_pep >> ind;
  f_ind_to_pep >> pep;
  while(!f_ind_to_pep.eof()){
      ind_to_pep[ind] = pep;
      f_ind_to_pep >> ind;
      f_ind_to_pep >> pep;
  }
  f_ind_to_pep.close();
  fname.str("");

  //ind_to_prot
  fname << in_dir << "/ind_to_prot";
  ifstream f_ind_to_prot(fname.str().c_str(),ios::binary);
 
  string prot;
  f_ind_to_prot >> ind;
  f_ind_to_prot >> prot;
  while(!f_ind_to_prot.eof()){
      ind_to_prot[ind] = prot;
      f_ind_to_prot >> ind;
      f_ind_to_prot >> prot;
  }
  f_ind_to_prot.close();
  fname.str("");

  //protind_to_length
  fname << in_dir << "/protind_to_length";
  ifstream f_protind_to_length(fname.str().c_str(),ios::binary);
  if(!f_protind_to_length.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  protind_to_length = new int[num_prot];
  f_protind_to_length.read((char*)protind_to_length,sizeof(int)*num_prot);
  f_protind_to_length.close();
  fname.str("");


  //protind_to_label
  fname << in_dir << "/protind_to_label";
  ifstream f_protind_to_label(fname.str().c_str(),ios::binary);
  if(!f_protind_to_label.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  protind_to_label = new int[num_prot];
  f_protind_to_label.read((char*)protind_to_label,sizeof(int)*num_prot);
  f_protind_to_label.close();
  fname.str("");

  //psm data

  //psmind_to_pepind
  fname << in_dir << "/psmind_to_pepind";
  ifstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_pepind.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_pepind = new int[num_psms];
  f_psmind_to_pepind.read((char*)psmind_to_pepind,sizeof(int)*num_psms);
  f_psmind_to_pepind.close();
  fname.str("");
  
  //psmind_to_scan
  fname << in_dir << "/psmind_to_scan";
  ifstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_scan.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_scan = new int[num_psms];
  f_psmind_to_scan.read((char*)psmind_to_scan,sizeof(int)*num_psms);
  f_psmind_to_scan.close();
  fname.str("");

  //psmind_to_charge
  fname << in_dir << "/psmind_to_charge";
  ifstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_charge.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_charge = new int[num_psms];
  f_psmind_to_charge.read((char*)psmind_to_charge,sizeof(int)*num_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //psmind_to_xcorr
  fname << in_dir << "/psmind_to_xcorr";
  ifstream f_psmind_to_xcorr(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_xcorr.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_xcorr = new double[num_psms];
  f_psmind_to_xcorr.read((char*)psmind_to_xcorr,sizeof(double)*num_psms);
  f_psmind_to_xcorr.close();
  fname.str("");

  //psmind_to_deltaCn
  fname << in_dir << "/psmind_to_deltaCn";
  ifstream f_psmind_to_deltaCn(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_deltaCn.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_deltaCn = new double[num_psms];
  f_psmind_to_deltaCn.read((char*)psmind_to_deltaCn,sizeof(double)*num_psms);
  f_psmind_to_deltaCn.close();
  fname.str("");
  
  //psmind_to_sp_score
  fname << in_dir << "/psmind_to_spscore";
  ifstream f_psmind_to_spscore(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_spscore.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_spscore = new double[num_psms];
  f_psmind_to_spscore.read((char*)psmind_to_spscore,sizeof(double)*num_psms);
  f_psmind_to_spscore.close();
  fname.str("");

  //psmind_to_calculated_mass
  fname << in_dir << "/psmind_to_calculated_mass";
  ifstream f_psmind_to_calculated_mass(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_calculated_mass.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_calculated_mass = new double[num_psms];
  f_psmind_to_calculated_mass.read((char*)psmind_to_calculated_mass,sizeof(double)*num_psms);
  f_psmind_to_calculated_mass.close();
  fname.str("");

  //psmind_to_precursor_mass
  fname << in_dir << "/psmind_to_precursor_mass";
  ifstream f_psmind_to_precursor_mass(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_precursor_mass.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_precursor_mass = new double[num_psms];
  f_psmind_to_precursor_mass.read((char*)psmind_to_precursor_mass,sizeof(double)*num_psms);
  f_psmind_to_precursor_mass.close();
  fname.str("");
  
  //fileind_to_fname
  fname << in_dir << "/fileind_to_fname";
  ifstream f_fileind_to_fname(fname.str().c_str(),ios::binary);
  if(!f_fileind_to_fname.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }

  string filename;
  f_fileind_to_fname >> ind;
  f_fileind_to_fname >> filename;
  while(!f_fileind_to_fname.eof()){
      fileind_to_fname[ind] = filename;
      f_fileind_to_fname >> ind;
      f_fileind_to_fname >> filename;
    }
  f_fileind_to_fname.close();
  fname.str("");
  
  //psmind_to_filename
  fname << in_dir << "/psmind_to_fileind";
  ifstream f_psmind_to_fileind(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_fileind.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_fileind = new int[num_psms];
  f_psmind_to_fileind.read((char*)psmind_to_fileind,sizeof(int)*num_psms);
  f_psmind_to_fileind.close();
  fname.str("");

  //psmind_to_sp_rank
  fname << in_dir << "/psmind_to_sp_rank";
  ifstream f_psmind_to_sp_rank(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_sp_rank.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_sp_rank = new int[num_psms];
  f_psmind_to_sp_rank.read((char*)psmind_to_sp_rank,sizeof(int)*num_psms);
  f_psmind_to_sp_rank.close();
  fname.str("");
  
  //psmind_to_xcorr_rank
  fname << in_dir << "/psmind_to_xcorr_rank";
  ifstream f_psmind_to_xcorr_rank(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_xcorr_rank.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_xcorr_rank = new int[num_psms];
  f_psmind_to_xcorr_rank.read((char*)psmind_to_xcorr_rank,sizeof(int)*num_psms);
  f_psmind_to_xcorr_rank.close();
  fname.str("");
    
  //psmind_to_match_spectrum
  fname << in_dir << "/psmind_to_matches_spectrum";
  ifstream f_psmind_to_matches_spectrum(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_matches_spectrum.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_matches_spectrum = new int[num_psms];
  f_psmind_to_matches_spectrum.read((char*)psmind_to_matches_spectrum,sizeof(int)*num_psms);
  f_psmind_to_matches_spectrum.close();
  fname.str("");
   
  //psmind_to_by_ions_matched
  fname << in_dir << "/psmind_to_by_ions_matched";
  ifstream f_psmind_to_by_ions_matched(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_by_ions_matched.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_by_ions_matched = new double[num_psms];
  f_psmind_to_by_ions_matched.read((char*)psmind_to_by_ions_matched,sizeof(double)*num_psms);
  f_psmind_to_by_ions_matched.close();
  fname.str("");
 
  //psmind_to_by_ions_total
  fname << in_dir << "/psmind_to_by_ions_total";
  ifstream f_psmind_to_by_ions_total(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_by_ions_total.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_by_ions_total = new double[num_psms];
  f_psmind_to_by_ions_total.read((char*)psmind_to_by_ions_total,sizeof(double)*num_psms);
  f_psmind_to_by_ions_total.close();
  fname.str("");
  
  //psmind_to_peptide_position
  fname << in_dir << "/psmind_to_peptide_position";
  ifstream f_psmind_to_peptide_position(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_peptide_position.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_peptide_position = new int[num_psms];
  f_psmind_to_peptide_position.read((char*)psmind_to_peptide_position,sizeof(int)*num_psms);
  f_psmind_to_peptide_position.close();
  fname.str("");
}

void Dataset :: clear_data_all_results()
{
  clear_data_psm_results();
  delete [] protind_to_label; protind_to_label = (int*)0;
  delete [] protind_to_length; protind_to_length = (int*)0;
  protind_to_pepinds.clear();
}




/******************************************************/

void Dataset :: normalize_psms()
{
  for (int i = 0; i < num_features; i++)
    {
      double mean = 0;
      for (int j = 0; j < num_psms; j++)
	mean += psmind_to_features[num_features*j+i];
      mean /= num_psms;
      
      double std = 0;
      for (int j = 0; j < num_psms; j++)
	{
	  psmind_to_features[num_features*j+i] -= mean;
	  std += psmind_to_features[num_features*j+i]*psmind_to_features[num_features*j+i];
	}
      std = sqrt(std/num_psms);

      for (int j = 0; j < num_psms; j++)
	{
	  if(std > 0)
	    psmind_to_features[num_features*j+i] /= std;
	}

      double sm = 0;
      for (int j = 0; j < num_psms; j++)
	{
	  sm += psmind_to_features[num_features*j+i]*psmind_to_features[num_features*j+i];
	}
      //cout << i << " " << sm/num_psms << endl;
    }
}


int Dataset :: print_features(string &filename)
{
  ofstream os(filename.c_str());
  if(!os.is_open())
    return 0;

  delete [] psmind_to_label; psmind_to_label = (int*)0;
  delete [] psmind_to_scan; psmind_to_scan = (int*)0;

  ostringstream fname;
  //psmind_to_label
  fname << in_dir << "/psmind_to_label";
  ifstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  psmind_to_label = new int[num_psms];
  f_psmind_to_label.read((char*)psmind_to_label,sizeof(int)*num_psms);
  f_psmind_to_label.close();
  fname.str("");
  
  //psmind_to_scan
  fname << in_dir << "/psmind_to_scan";
  ifstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  psmind_to_scan = new int[num_psms];
  f_psmind_to_scan.read((char*)psmind_to_scan,sizeof(int)*num_psms);
  f_psmind_to_scan.close();
  fname.str("");
  //print features header
  os<<"scan\t"<<"label\t";
  for(unsigned i=0;i<features_header_.size()-1;i++)
    os<<features_header_[i]<<"\t";   
  os<<features_header_[features_header_.size()-1]<<endl;

  for (int j = 0; j < num_psms; j++){
      os << psmind_to_scan[j] << "\t";
      os << psmind_to_label[j] << "\t"; 
      for (int i = 0; i < num_features-1; i++){
	os << psmind_to_features[num_features*j+i] << "\t"; 
      }
      os << psmind_to_features[num_features*j+(num_features-1)] << endl;
    }
  os.close();

  delete [] psmind_to_label; psmind_to_label = (int*)0;
  delete [] psmind_to_scan; psmind_to_scan = (int*)0;
  return 1;
}

/*****************************************************/
void Dataset :: load_data_pep_training()
{

  ostringstream fname;
  fname << in_dir << "/summary";
  ifstream f_summary(fname.str().c_str());
  if(!f_summary.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary >> num_pep;
  f_summary >> num_pos_pep;
  f_summary >> num_neg_pep;
  f_summary >> num_prot;
  f_summary >> num_pos_prot;
  f_summary >> num_neg_prot;
  f_summary.close();
  fname.str("");
  
  //psm features
  fname << in_dir << "/" << "psm";
  ifstream f_psm_feat(fname.str().c_str(),ios::binary);
  if(!f_psm_feat.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_features = new double[num_psms*num_features];
  f_psm_feat.read((char*)psmind_to_features,sizeof(double)*num_psms*num_features);
  f_psm_feat.close();
  fname.str("");


  //pepind_to_psminds
  fname << in_dir << "/pepind_to_psminds";
  ifstream f_pepind_to_psminds(fname.str().c_str(),ios::binary);
  if(!f_pepind_to_psminds.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  pepind_to_psminds.load(f_pepind_to_psminds);
  f_pepind_to_psminds.close();
  fname.str("");
  
}


void Dataset :: clear_data_pep_training()
{
  delete [] psmind_to_features; psmind_to_features = (double*)0;
}

void Dataset :: load_labels_pep_training()
{

  ostringstream fname;
  fname << in_dir << "/summary";
  ifstream f_summary(fname.str().c_str());
  if(!f_summary.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  f_summary >> num_features;
  f_summary >> num_psms;
  f_summary >> num_pos_psms;
  f_summary >> num_neg_psms;
  f_summary >> num_pep;
  f_summary >> num_pos_pep;
  f_summary >> num_neg_pep;
  f_summary >> num_prot;
  f_summary >> num_pos_prot;
  f_summary >> num_neg_prot;
  f_summary.close();
  fname.str("");

  //psmind_to_label
  fname << in_dir << "/psmind_to_label";
  ifstream f_psmind_to_label(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_label.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_label = new int[num_psms];
  f_psmind_to_label.read((char*)psmind_to_label,sizeof(int)*num_psms);
  f_psmind_to_label.close();
  fname.str("");
  
  //pepind_to_label
  fname << in_dir << "/pepind_to_label";
  ifstream f_pepind_to_label(fname.str().c_str(),ios::binary);
  if(!f_pepind_to_label.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  pepind_to_label = new int[num_pep];
  f_pepind_to_label.read((char*)pepind_to_label,sizeof(int)*num_pep);
  f_pepind_to_label.close();
  fname.str("");
}

void Dataset :: clear_labels_pep_training()
{
  delete [] psmind_to_label; psmind_to_label = (int*)0;
  delete [] pepind_to_label; pepind_to_label = (int*)0;
}

void Dataset :: load_data_pep_results()
{

  ostringstream fname;

  //ind_to_pep
  fname << in_dir << "/ind_to_pep";
  ifstream f_ind_to_pep(fname.str().c_str(),ios::binary);
  int ind;
  string pep;
  f_ind_to_pep >> ind;
  f_ind_to_pep >> pep;
  while(!f_ind_to_pep.eof()){
      ind_to_pep[ind] = pep;
      f_ind_to_pep >> ind;
      f_ind_to_pep >> pep;
  }
  f_ind_to_pep.close();
  fname.str("");

  //pepind_to_protinds
  fname << in_dir << "/pepind_to_protinds";
  ifstream f_pepind_to_protinds(fname.str().c_str(),ios::binary);
  if(!f_pepind_to_protinds.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  pepind_to_protinds.load(f_pepind_to_protinds);
  f_pepind_to_protinds.close();
  fname.str("");

  //ind_to_prot
  fname << in_dir << "/ind_to_prot";
  ifstream f_ind_to_prot(fname.str().c_str(),ios::binary);
 
  string prot;
  f_ind_to_prot >> ind;
  f_ind_to_prot >> prot;
  while(!f_ind_to_prot.eof()){
      ind_to_prot[ind] = prot;
      f_ind_to_prot >> ind;
      f_ind_to_prot >> prot;
  }
  f_ind_to_prot.close();
  fname.str("");

  //protind_to_length
  fname << in_dir << "/protind_to_length";
  ifstream f_protind_to_length(fname.str().c_str(),ios::binary);
  if(!f_protind_to_length.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  protind_to_length = new int[num_prot];
  f_protind_to_length.read((char*)protind_to_length,sizeof(int)*num_prot);
  f_protind_to_length.close();
  fname.str("");


  //protind_to_label
  fname << in_dir << "/protind_to_label";
  ifstream f_protind_to_label(fname.str().c_str(),ios::binary);
  if(!f_protind_to_label.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  protind_to_label = new int[num_prot];
  f_protind_to_label.read((char*)protind_to_label,sizeof(int)*num_prot);
  f_protind_to_label.close();
  fname.str("");



  //psm data

  //psmind_to_pepind
  fname << in_dir << "/psmind_to_pepind";
  ifstream f_psmind_to_pepind(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_pepind.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_pepind = new int[num_psms];
  f_psmind_to_pepind.read((char*)psmind_to_pepind,sizeof(int)*num_psms);
  f_psmind_to_pepind.close();
  fname.str("");
  
  //psmind_to_scan
  fname << in_dir << "/psmind_to_scan";
  ifstream f_psmind_to_scan(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_scan.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_scan = new int[num_psms];
  f_psmind_to_scan.read((char*)psmind_to_scan,sizeof(int)*num_psms);
  f_psmind_to_scan.close();
  fname.str("");

  //psmind_to_charge
  fname << in_dir << "/psmind_to_charge";
  ifstream f_psmind_to_charge(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_charge.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_charge = new int[num_psms];
  f_psmind_to_charge.read((char*)psmind_to_charge,sizeof(int)*num_psms);
  f_psmind_to_charge.close();
  fname.str("");

  //psmind_to_xcorr
  fname << in_dir << "/psmind_to_xcorr";
  ifstream f_psmind_to_xcorr(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_xcorr.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_xcorr = new double[num_psms];
  f_psmind_to_xcorr.read((char*)psmind_to_xcorr,sizeof(double)*num_psms);
  f_psmind_to_xcorr.close();
  fname.str("");

  //psmind_to_deltaCn
  fname << in_dir << "/psmind_to_deltaCn";
  ifstream f_psmind_to_deltaCn(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_deltaCn.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_deltaCn = new double[num_psms];
  f_psmind_to_deltaCn.read((char*)psmind_to_deltaCn,sizeof(double)*num_psms);
  f_psmind_to_deltaCn.close();
  fname.str("");
  
  //psmind_to_sp_score
  fname << in_dir << "/psmind_to_spscore";
  ifstream f_psmind_to_spscore(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_spscore.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_spscore = new double[num_psms];
  f_psmind_to_spscore.read((char*)psmind_to_spscore,sizeof(double)*num_psms);
  f_psmind_to_spscore.close();
  fname.str("");

  //psmind_to_calculated_mass
  fname << in_dir << "/psmind_to_calculated_mass";
  ifstream f_psmind_to_calculated_mass(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_calculated_mass.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_calculated_mass = new double[num_psms];
  f_psmind_to_calculated_mass.read((char*)psmind_to_calculated_mass,sizeof(double)*num_psms);
  f_psmind_to_calculated_mass.close();
  fname.str("");

  //psmind_to_precursor_mass
  fname << in_dir << "/psmind_to_precursor_mass";
  ifstream f_psmind_to_precursor_mass(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_precursor_mass.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_precursor_mass = new double[num_psms];
  f_psmind_to_precursor_mass.read((char*)psmind_to_precursor_mass,sizeof(double)*num_psms);
  f_psmind_to_precursor_mass.close();
  fname.str("");
  
  //fileind_to_fname
  fname << in_dir << "/fileind_to_fname";
  ifstream f_fileind_to_fname(fname.str().c_str(),ios::binary);
  if(!f_fileind_to_fname.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }

  string filename;
  f_fileind_to_fname >> ind;
  f_fileind_to_fname >> filename;
  while(!f_fileind_to_fname.eof()){
      fileind_to_fname[ind] = filename;
      f_fileind_to_fname >> ind;
      f_fileind_to_fname >> filename;
    }
  f_fileind_to_fname.close();
  fname.str("");
  
  //psmind_to_filename
  fname << in_dir << "/psmind_to_fileind";
  ifstream f_psmind_to_fileind(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_fileind.is_open())
    {
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_fileind = new int[num_psms];
  f_psmind_to_fileind.read((char*)psmind_to_fileind,sizeof(int)*num_psms);
  f_psmind_to_fileind.close();
  fname.str("");

  //psmind_to_sp_rank
  fname << in_dir << "/psmind_to_sp_rank";
  ifstream f_psmind_to_sp_rank(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_sp_rank.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
    }
  psmind_to_sp_rank = new int[num_psms];
  f_psmind_to_sp_rank.read((char*)psmind_to_sp_rank,sizeof(int)*num_psms);
  f_psmind_to_sp_rank.close();
  fname.str("");
  
  //psmind_to_xcorr_rank
  fname << in_dir << "/psmind_to_xcorr_rank";
  ifstream f_psmind_to_xcorr_rank(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_xcorr_rank.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_xcorr_rank = new int[num_psms];
  f_psmind_to_xcorr_rank.read((char*)psmind_to_xcorr_rank,sizeof(int)*num_psms);
  f_psmind_to_xcorr_rank.close();
  fname.str("");
    
  //psmind_to_match_spectrum
  fname << in_dir << "/psmind_to_matches_spectrum";
  ifstream f_psmind_to_matches_spectrum(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_matches_spectrum.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_matches_spectrum = new int[num_psms];
  f_psmind_to_matches_spectrum.read((char*)psmind_to_matches_spectrum,sizeof(int)*num_psms);
  f_psmind_to_matches_spectrum.close();
  fname.str("");
   
  //psmind_to_by_ions_matched
  fname << in_dir << "/psmind_to_by_ions_matched";
  ifstream f_psmind_to_by_ions_matched(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_by_ions_matched.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_by_ions_matched = new double[num_psms];
  f_psmind_to_by_ions_matched.read((char*)psmind_to_by_ions_matched,sizeof(double)*num_psms);
  f_psmind_to_by_ions_matched.close();
  fname.str("");
 
  //psmind_to_by_ions_total
  fname << in_dir << "/psmind_to_by_ions_total";
  ifstream f_psmind_to_by_ions_total(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_by_ions_total.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_by_ions_total = new double[num_psms];
  f_psmind_to_by_ions_total.read((char*)psmind_to_by_ions_total,sizeof(double)*num_psms);
  f_psmind_to_by_ions_total.close();
  fname.str("");
  
  //psmind_to_peptide_position
  fname << in_dir << "/psmind_to_peptide_position";
  ifstream f_psmind_to_peptide_position(fname.str().c_str(),ios::binary);
  if(!f_psmind_to_peptide_position.is_open()){
      cout << "could not open file " << fname.str() <<  " for reading data\n";
      return;
  }
  psmind_to_peptide_position = new int[num_psms];
  f_psmind_to_peptide_position.read((char*)psmind_to_peptide_position,sizeof(int)*num_psms);
  f_psmind_to_peptide_position.close();
  fname.str("");
}

void Dataset :: clear_data_pep_results()
{
  clear_data_psm_results();
}



/*
 * Local Variable 
 * mode: c
 * c-basic-offset: 2
 * End: 
 */

/*
int main()
{
  Dataset* d = new Dataset();
  
  d->set_input_dir("yeast");
  d->load_prot_data();
  //d->normalize_psms();
  return 0;
}
*/
