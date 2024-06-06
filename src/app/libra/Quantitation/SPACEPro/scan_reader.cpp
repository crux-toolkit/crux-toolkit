#include "scan_reader.h"

using namespace std;
using namespace mzParser;

bool scan_reader::getTIC(peptide_lists& my_peptide_lists, my_parameters& my_params){
  BasicSpectrum mySpec;
  MzParser myfile(&mySpec);

  if (!myfile.load(my_params.mzml.c_str())) {
    string local = my_params.mzml;
    size_t ret = local.find_last_of('\\');
    if (ret == string::npos) local.find_last_of('/');
    if (ret == string::npos) return false;
    local = local.substr(ret + 1, local.size());
    cout << "WARNING: " << my_params.mzml << " not found. Attempting: " << local << endl;
    if (!myfile.load(local.c_str())) return false;
  }

  size_t i=0;
  while (myfile.readSpectrum()) { //this line checks that the current spectrum is valid. The moment it is not, the end of file was reached.
    if (mySpec.getMSLevel() != 2 ) continue;

    while(i<my_peptide_lists.all_psm.size() && my_peptide_lists.all_psm[i].scan<mySpec.getScanNum()) i++;
    if(i==my_peptide_lists.all_psm.size()) break;
    if(my_peptide_lists.all_psm[i].scan == mySpec.getScanNum()) my_peptide_lists.all_psm[i].tic=mySpec.getTotalIonCurrent();
  }
  return true;

}

bool scan_reader::mzml(peptide_lists& my_peptide_lists, my_parameters& my_params) {

	//READ MZML FILE FOR BOTH FULLY TRYPTIC AND MISCLEAVED PEPTIDES (MZML)
	BasicSpectrum mySpec;
  MzParser myfile(&mySpec);

	if(!myfile.load(my_params.mzml.c_str())) {
    string local= my_params.mzml;
    size_t ret = local.find_last_of('\\');
    if(ret==string::npos) local.find_last_of('/');
    if(ret==string::npos) return false;
    local=local.substr(ret+1,local.size());
    cout << "WARNING: " << my_params.mzml << " not found. Attempting: " << local << endl;
    if (!myfile.load(local.c_str())) return false;
  }

	//MH: We need to get through precursor peak extraction fast.
	//To do so, we need to do it a) using a single pass through all arrays, and
	//b) using fast lookups instead of full iteration.
	//Step 1 is to sort all PSMs by retention time
	sort(my_peptide_lists.all_psm.begin(), my_peptide_lists.all_psm.end(), compareRTime);
	size_t pepIndex = 0; // we will start from the first sorted peptide

	//MH: Start reading our mzML file.
	while (myfile.readSpectrum()) { //this line checks that the current spectrum is valid. The moment it is not, the end of file was reached.
    if(mySpec.getMSLevel()>1) continue;

		size_t i = pepIndex;
		while (i < my_peptide_lists.all_psm.size() && (my_peptide_lists.all_psm[i].xml_rtime / 60) < (mySpec.getRTime() - my_params.ret_time)) i++;
		if (i == my_peptide_lists.all_psm.size()) break; //if we've checked every peptide, stop now.
		pepIndex = i; //mark our new start point for the next iteration

		while (i < my_peptide_lists.all_psm.size() && (my_peptide_lists.all_psm[i].xml_rtime / 60) < (mySpec.getRTime() + my_params.ret_time)) {
      //match compensation voltage
      if(fabs(my_peptide_lists.all_psm[i].compensation_voltage-mySpec.getCompensationVoltage())>0.1) {
        i++; //go to the next PSM
        continue;
      }

			//compute the desired m/z value to find, and a tolerance around that value
			double mz = (my_peptide_lists.all_psm[i].pre_neutral_mass + my_peptide_lists.all_psm[i].charge * 1.007276466) / my_peptide_lists.all_psm[i].charge;
			double tolerance = my_params.ppm * mz / 1000000;
			int ret = findPeakMZ(mySpec, mz, tolerance);
			dsXIC dsx;
			dsx.rTime = mySpec.getRTime();
			dsx.tot = 0;
			if (ret > -1) {
				dsx.intensity = (float)mySpec[ret].intensity;
        //printf("Match: %.8lf +/- %.8lf\t%.8lf\t%.8lf\t%d\n",mz,tolerance,mySpec[ret].mz,mz- mySpec[ret].mz,ret);
			} else {
				dsx.intensity = 0;
			}
			my_peptide_lists.all_psm[i].XIC.push_back(dsx);
			i++;//MH: go to the next PSM
		}
		//myfile.readFile(NULL, mySpec); //read the next spectrum
	}

	return true;

}



//MH: Function to binary search for an mz value within a specific tolerance
int scan_reader::findPeakMZ(mzParser::BasicSpectrum& spec, double mz, double tol){
  size_t sz = spec.size();
  size_t lower = 0;
  size_t mid = sz / 2;
  size_t upper = sz;

  //our boundaries
  double LB=mz-tol;
  double UB=mz+tol;

  //stop now if the spectrum is empty.
  if(sz<1) return -1;

  while(lower<upper){
    if(spec[mid].mz>UB){ //too high
      upper = mid;
      mid = (lower + upper) / 2;
    } else if(spec[mid].mz<LB) { //too low
      lower = mid + 1;
      mid = (lower + upper) / 2;
    } else { //we have a match, now we must step left and right to make sure it's the max
      double max=spec[mid].intensity;
      size_t maxI=mid;
      size_t i=mid;
      while(i>0 && spec[i-1].mz>LB){
        i--;
        if(spec[i].intensity>max){
          max=spec[i].intensity;
          maxI=i;
        }
      }
      i=mid;
      while(i<sz-1 && spec[i+1].mz<UB){
        i++;
        if (spec[i].intensity > max) {
          max = spec[i].intensity;
          maxI = i;
        }
      }
      return (int)maxI; 
    }
  }
  return -1; //-1 means can't find peak;
}
