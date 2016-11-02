#include "GlobalParams.h"
#include "util/Params.h"
#include "parameter.h"
#include "util/StringUtils.h"

using namespace std;


MASS_TYPE_T GlobalParams::isotopic_mass_;
int GlobalParams::missed_cleavages_;
int GlobalParams::max_aas_modified_;
FLOAT_T GlobalParams::min_mass_;
FLOAT_T GlobalParams::max_mass_;
WINDOW_TYPE_T GlobalParams::precursor_window_type_;
FLOAT_T GlobalParams::precursor_window_;
int GlobalParams::min_length_;
int GlobalParams::max_length_;
string GlobalParams::xlink_prevents_cleavage_;
string GlobalParams::max_ion_charge_;
ION_TYPE_T GlobalParams::primary_ions_;
MASS_TYPE_T GlobalParams::fragment_mass_;
bool GlobalParams::precursor_ions_;
ENZYME_T GlobalParams::enzyme_;
DIGEST_T GlobalParams::digestion_;
FLOAT_T GlobalParams::remove_precursor_tolerance_;
OBSERVED_PREPROCESS_STEP_T GlobalParams::stop_after_;
bool GlobalParams::xlink_include_inter_;
bool GlobalParams::xlink_include_intra_;
bool GlobalParams::xlink_include_inter_intra_;
bool GlobalParams::xlink_include_deadends_;
bool GlobalParams::xlink_include_linears_;
bool GlobalParams::xlink_include_selfloops_;
int GlobalParams::max_xlink_mods_;
int GlobalParams::mod_precision_;
int GlobalParams::xlink_top_n_;
vector<int> GlobalParams::isotope_windows_;
FLOAT_T GlobalParams::fraction_to_fit_;
bool GlobalParams::xlink_use_ion_cache_;
MASS_FORMAT_T GlobalParams::mod_mass_format_;

void GlobalParams::set() {
  isotopic_mass_ = get_mass_type_parameter("isotopic-mass");
  missed_cleavages_ = Params::GetInt("missed-cleavages");
  max_aas_modified_ = Params::GetInt("max-aas-modified");
  min_mass_ = Params::GetDouble("min-mass");
  max_mass_ = Params::GetDouble("max-mass");
  precursor_window_type_ = string_to_window_type(Params::GetString("precursor-window-type"));
  precursor_window_ = Params::GetDouble("precursor-window");
  min_length_ = Params::GetInt("min-length");
  max_length_ = Params::GetInt("max-length");
  xlink_prevents_cleavage_ = Params::GetString("xlink-prevents-cleavage");
  max_ion_charge_ = Params::GetString("max-ion-charge");
  string_to_ion_type(Params::GetString("primary-ions"), &primary_ions_);
  fragment_mass_ = get_mass_type_parameter("fragment-mass");
  precursor_ions_ = Params::GetBool("precursor-ions");
  enzyme_ = get_enzyme_type_parameter("enzyme");
  digestion_ = get_digest_type_parameter("digestion");
  remove_precursor_tolerance_ = Params::GetDouble("remove-precursor-tolerance");
  stop_after_ = string_to_observed_preprocess_step(Params::GetString("stop-after"));
  xlink_include_inter_ = Params::GetBool("xlink-include-inter");
  xlink_include_intra_ = Params::GetBool("xlink-include-intra");
  xlink_include_inter_intra_ = Params::GetBool("xlink-include-inter-intra");
  xlink_include_deadends_ = Params::GetBool("xlink-include-deadends");
  xlink_include_selfloops_ = Params::GetBool("xlink-include-selfloops");
  xlink_include_linears_ = Params::GetBool("xlink-include-linears");
  max_xlink_mods_ = Params::GetInt("max-xlink-mods");
  mod_precision_ = Params::GetInt("mod-precision");
  xlink_top_n_ = Params::GetInt("xlink-top-n");
  isotope_windows_ = StringUtils::Split<int>(Params::GetString("isotope-windows"), ',');
  fraction_to_fit_ = Params::GetDouble("fraction-top-scores-to-fit");
  xlink_use_ion_cache_ = Params::GetBool("xlink-use-ion-cache");
  mod_mass_format_ = get_mass_format_type_parameter("mod-mass-format");
}

const MASS_TYPE_T& GlobalParams::getIsotopicMass() {
  return isotopic_mass_;
}

const int& GlobalParams::getMissedCleavages() {
  return missed_cleavages_;
}

const int& GlobalParams::getMaxAasModified() {
  return max_aas_modified_;
}

const FLOAT_T& GlobalParams::getMinMass() {
  return min_mass_;
}

const FLOAT_T& GlobalParams::getMaxMass() {
  return max_mass_;
}

const WINDOW_TYPE_T& GlobalParams::getPrecursorWindowType() {
  return precursor_window_type_;
}

const FLOAT_T& GlobalParams::getPrecursorWindow() {
  return precursor_window_;
}

const int& GlobalParams::getMinLength() {
  return min_length_;
}

const int& GlobalParams::getMaxLength() {
  return max_length_;
}

const string& GlobalParams::getXLinkPreventsCleavage() {
  return xlink_prevents_cleavage_;
}

const string& GlobalParams::getMaxIonCharge() {
  return max_ion_charge_;
}

const ION_TYPE_T& GlobalParams::getPrimaryIons() {
  return primary_ions_;
}

const MASS_TYPE_T& GlobalParams::getFragmentMass() {
  return fragment_mass_;
}

const bool& GlobalParams::getPrecursorIons() {
  return precursor_ions_;
}

const ENZYME_T& GlobalParams::getEnzyme() {
  return enzyme_;
}

const DIGEST_T& GlobalParams::getDigestion() {
  return digestion_;
}

const FLOAT_T& GlobalParams::getRemovePrecursorTolerance() {
  return remove_precursor_tolerance_;
}

const OBSERVED_PREPROCESS_STEP_T& GlobalParams::getStopAfter() {
  return stop_after_;
}

const bool& GlobalParams::getXLinkIncludeInter() {
  return xlink_include_inter_;
}
  
const bool& GlobalParams::getXLinkIncludeIntra() {
  return xlink_include_intra_;
}

const bool& GlobalParams::getXLinkIncludeInterIntra() {
  return xlink_include_inter_intra_;
}

const bool& GlobalParams::getXLinkIncludeDeadends() {
  return xlink_include_deadends_;
}
 
const bool& GlobalParams::getXLinkIncludeSelfloops() {
  return xlink_include_selfloops_;
}

const bool& GlobalParams::getXLinkIncludeLinears() {
  return xlink_include_linears_;
}
const int& GlobalParams::getMaxXLinkMods() {
  return max_xlink_mods_;
}

const int& GlobalParams::getModPrecision() {
  return mod_precision_;
}

const int& GlobalParams::getXLinkTopN() {
  return xlink_top_n_;
}

const vector<int>& GlobalParams::getIsotopeWindows() {
  return isotope_windows_;
}

const FLOAT_T& GlobalParams::getFractionToFit() {
  return fraction_to_fit_;
}

const bool& GlobalParams::getXLinkUseIonCache() {
  return xlink_use_ion_cache_;
}

const MASS_FORMAT_T& GlobalParams::getModMassFormat() {
  return mod_mass_format_;
}

