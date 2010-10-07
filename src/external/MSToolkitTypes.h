#ifndef _MSTOOLKITTYPES_H
#define _MSTOOLKITTYPES_H

namespace MSToolkit{
enum MSFileType {
  MS1,
  MS2,
  MS3,
  ZS,
  UZS,
  IonSpec,
	SRM,
  Unknown
};

enum MSFileFormat {
	bms1,
	bms2,
	cms1,
	cms2,
	ms1,
	ms2,
	mzXML,
	zs,
	uzs,
	dunno
};

enum MSTag {
  no,
  D,
  H,
  I,
  S,
  Z
};

struct MSHeader {
	char header[16][128];
};

struct MSScanInfo {
	int scanNumber[2];
	float rTime;
	int numDataPoints;
	int numZStates;
};

struct Peak_T {
  double mz;
  float intensity;
};

struct ZState {
  int z;
  double mz;
};

} // namespace
#endif


