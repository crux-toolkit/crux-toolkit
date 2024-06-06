#include "CnprModAminoacidMass.h"

using namespace std;

CnprModAminoacidMass::CnprModAminoacidMass(){
  position=0;
  mass=0;
}

void CnprModAminoacidMass::write(FILE* f, int tabs){
  //required
  string el = "mod_aminoacid_mass";
  if (position==0) NPRerrMsg(el, "position");
  if (mass=0) NPRerrMsg(el, "mass");

  NPRprintTabs(f, tabs);
  fprintf(f, "<mod_aminoacid_mass");
  fprintf(f, " position=\"%d\"", position);
  fprintf(f, " mass=\"%.6lf\"", mass);
  fprintf(f, "/>\n");

}
