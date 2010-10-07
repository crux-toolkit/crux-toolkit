#ifndef _MSREADER_H
#define _MSREADER_H

#include "Spectrum.h"
#include "MSObject.h"
#include <cstring>
#include "zlib.h"
#include "ramp.h"
#include "base64.h"

//Macro for general large file support - RAMP toolkit style
#ifndef _LARGEFILE_SOURCE // use MSFT API for 64 bit file pointers
typedef fpos_t f_off;
#else // a real OS with real file handling
#define _FILE_OFFSET_BITS 64
typedef off_t f_off;
#define fopen(a,b) fopen64(a,b)
#endif 

using namespace std;

namespace MSToolkit{

class MSReader {
 public:
  //Constructors & Destructors
  MSReader();
	~MSReader();

  //Functions
	void appendFile(char* c, bool text, Spectrum& s);
	void appendFile(char* c, bool text, MSObject& m);
  MSHeader& getHeader();
	//Spectrum readBinaryFile(char* c, Spectrum& s, int scNum=0);
  //Spectrum readMSFile(char* c,int scNum=0);
	
  MSFileType getFileType();
  int getPercent();
	void setPrecision(int i, int j);
	void setPrecisionInt(int i);
	void setPrecisionMZ(int i);
  void writeFile(char* c, bool text, MSObject& m);

	bool readFile(char* c, bool text, Spectrum& s, int scNum=0);
	bool readFile(char* c, MSFileFormat f, Spectrum& s, int scNum=0);
	void setFilter(MSFileType m);

	//File compression
	void setCompression(bool b);


 protected:

 private:
  //Data Members
  FILE *fileIn;
  MSHeader header;
	int headerIndex;
  MSFileType fileType;
  f_off lEnd;
  f_off lPivot;
  f_off lFWidth;
	int iIntensityPrecision;
	int iMZPrecision;

	//File compression
	bool compressMe;

	//mzXML support variables;
  ramp_fileoffset_t  *pScanIndex;
	RAMPFILE  *rampFileIn;
	bool rampFileOpen;
	int rampLastScan;
	int rampIndex;
	MSFileType filter;

  //Functions
  void closeFile();
  int openFile(char* c, bool text=false);
  bool findSpectrum(int i);
	void readCompressSpec(FILE* fileIn, MSScanInfo& ms, Spectrum& s);

	void writeBinarySpec(FILE* fileOut, Spectrum& s);
	void writeCompressSpec(FILE* fileOut, Spectrum& s);
	void writeTextSpec(FILE* fileOut, Spectrum& s);
	void writeSpecHeader(FILE* fileOut, bool text, Spectrum& s);
};

} // namespace
#endif

