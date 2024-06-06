#include "CnprProteinSummaryDataFilter.h"

using namespace std;

CnprProteinSummaryDataFilter::CnprProteinSummaryDataFilter(){
  min_probability=-1;
  sensitivity = -1;
  false_positive_error_rate = -1;
  predicted_num_correct = -1;
  predicted_num_incorrect = -1;
}

void CnprProteinSummaryDataFilter::write(FILE* f, int tabs){
  //required
  string el = "protein_summary_data_filter";
  if (min_probability<0) NPRerrMsg(el, "min_probability");
  if (sensitivity<0) NPRerrMsg(el, "sensitivity");
  if (false_positive_error_rate<0) NPRerrMsg(el, "false_positive_error_rate");

  NPRprintTabs(f, tabs);
  fprintf(f, "<protein_summary_data_filter");
  fprintf(f, " error=\"%.2lf\"", min_probability);
  fprintf(f, " min_prob=\"%.3lf\"", sensitivity);
  fprintf(f, " num_corr=\"%.3lf\"", false_positive_error_rate);
  if (predicted_num_correct>-1) fprintf(f, " predicted_num_correct=\"%.0lf\"", predicted_num_correct);
  if (predicted_num_incorrect>-1) fprintf(f, " predicted_num_incorrect=\"%.0lf\"", predicted_num_incorrect);
  fprintf(f, "/>\n");

}
