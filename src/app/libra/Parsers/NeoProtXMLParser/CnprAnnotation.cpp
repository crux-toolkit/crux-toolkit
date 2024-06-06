#include "CnprAnnotation.h"

using namespace std;

void CnprAnnotation::write(FILE* f, int tabs){
  //required
  string el = "annotation";
  if (protein_description.empty()) NPRerrMsg(el, "protein_description");

  NPRprintTabs(f, tabs);
  fprintf(f, "<annotation");
  fprintf(f, " protein_description=\"%s\"", protein_description.c_str());
  if (!ipi_name.empty()) fprintf(f, " ipi_name=\"%s\"", ipi_name.c_str());
  if (!refseq_name.empty()) fprintf(f, " refseq_name=\"%s\"", refseq_name.c_str());
  if (!swissprot_name.empty()) fprintf(f, " swissprot_name=\"%s\"", swissprot_name.c_str());
  if (!ensembl_name.empty()) fprintf(f, " ensembl_name=\"%s\"", ensembl_name.c_str());
  if (!trembl_name.empty()) fprintf(f, " trembl_name=\"%s\"", trembl_name.c_str());
  if (!locus_link_name.empty()) fprintf(f, " locus_link_name=\"%s\"", locus_link_name.c_str());
  if (!flybase.empty()) fprintf(f, " flybase=\"%s\"", flybase.c_str());
  fprintf(f, "/>\n");

}
