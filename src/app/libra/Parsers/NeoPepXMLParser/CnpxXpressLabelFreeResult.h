#ifndef _CNPXXPRESSLABELFREERESULT_H
#define _CNPXXPRESSLABELFREERESULT_H

#include <string>
#include <vector>

class CnpxXpressLabelFreeResult {
public:

  CnpxXpressLabelFreeResult();
  CnpxXpressLabelFreeResult(bool b);

  bool present();
  void write(FILE* f);

  int first_scan;
  int last_scan;
  int peak_intensity_scan;
  int charge;

  std::string first_scan_RT_seconds;
  std::string last_scan_RT_seconds;
  std::string precursor_mz;
  std::string peak_area;
  std::string peak_intensity;
  std::string peak_intensity_RT_seconds;

private:
  bool active;

};

#endif
