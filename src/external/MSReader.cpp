#include "MSReader.h"
#include <iostream>
using namespace std;
using namespace MSToolkit;

MSReader::MSReader(){
	fileIn=NULL;
	iIntensityPrecision=1;
	iMZPrecision=4;
	filter=Unknown;
	rampFileOpen=false;
	compressMe=false;
};

MSReader::~MSReader(){
	closeFile();
	if(rampFileOpen) {
		rampCloseFile(rampFileIn);
		free(pScanIndex);
	};

};

void MSReader::closeFile(){
  if(fileIn!=NULL) fclose(fileIn);
};

MSHeader& MSReader::getHeader(){
  return header;
};

/* 0 = File opened correctly
   1 = Could not open file
*/
int MSReader::openFile(char *c,bool text){
	int i;

	//for finding the end of a large file with MS
	#ifndef _LARGEFILE_SOURCE
	if(text) i=_open(c,_O_TEXT);
	else i=_open(c,_O_BINARY);
	if(i>-1) {
		_lseeki64(i,0,2);
		//printf("%I64d\n",_telli64(i));
		lEnd=_telli64(i);
		//printf("%I64d %I64d\n",lEnd,_telli64(i));
		_close(i);
	};
	#endif

	if(text) fileIn=fopen(c,"rt");
	else fileIn=fopen(c,"rb");

  if(fileIn==NULL) {
		for(i=0;i<16;i++) strcpy(header.header[i],"\0");
    headerIndex=0;
    fileType=Unknown;
    return 1;
  } else {
    fileType=Unknown;

		//if we don't have the eof position, get it here.
		#ifdef _LARGEFILE_SOURCE
    fseeko(fileIn,0,2);
    lEnd = ftello(fileIn);
		//fgetpos(fileIn,&lEnd);
		#endif
    
		lPivot = 0;
    lFWidth = lEnd/2;
    
		#ifdef _LARGEFILE_SOURCE
		fseeko(fileIn,0,0);
		#else
		fsetpos(fileIn,&lPivot);
	  #endif
	
		if(text){
			for(i=0;i<16;i++) strcpy(header.header[i],"\0");
			headerIndex=0;
		} else {
			fread(&header,sizeof(MSHeader),1,fileIn);
		};

	  return 0;
  };
};

MSFileType MSReader::getFileType(){
  return fileType;
};

/*
Spectrum MSReader::readBinaryFile(char *c, Spectrum& s, int scNum){
	MS1ScanInfo ms;
	Peak_T p;
	int i;

	s.clear();

	if(c!=NULL){
		closeFile();
		if(openFile(c,false)==1) return s;
	} else if(fileIn==NULL) {
		return s;
	};


	fread(&ms,sizeof(MS1ScanInfo),1,fileIn);
	if(scNum!=0){
		while(ms.scanNumber[0]!=scNum){
			fseek(fileIn,ms.numDataPoints*12,1);
			fread(&ms,sizeof(MS1ScanInfo),1,fileIn);
			if(feof(fileIn)) return s;
		};
	};
	if(feof(fileIn)) return s;

	s.setScanNumber(ms.scanNumber[0]);
	s.setRTime(ms.rTime);
	for(i=0;i<ms.numDataPoints;i++){
		fread(&p.mz,8,1,fileIn);
		fread(&p.intensity,4,1,fileIn);
		s.add(p);
	};

	return s;

};
*/

bool MSReader::readFile(char *c, bool text, Spectrum& s, int scNum){
	MSScanInfo ms;
	Peak_T p;
	ZState z;
	int i;

	//variables for text reading only
	bool firstScan = false;
  bool bScan = true;
  bool bDoneHeader = false;
  char tstr[256];
  char ch;
	char *tok;
	f_off fpoint;

	//variables for compressed files
	uLong mzLen, intensityLen;

	//clear any spectrum data
	s.clear();

	//check for valid file and if we can access it
	if(c!=NULL){
		closeFile();
		if(openFile(c,text)==1) return false;
	} else if(fileIn==NULL) {
		return false;
	};

	//Handle binary and text files differently
	if(!text){

		//if binary file, read scan info sequentially, skipping to next scan if requested
		fread(&ms,sizeof(MSScanInfo),1,fileIn);
		//cout << "F2:" << ftell(fileIn) << " " << ms.scanNumber[0] << endl;

		if(scNum!=0){

			fpoint=sizeof(MSHeader);
			#ifdef _LARGEFILE_SOURCE
			fseeko(fileIn,sizeof(MSHeader),0);
			#else
			fsetpos(fileIn,&fpoint);
			#endif

			fread(&ms,sizeof(MSScanInfo),1,fileIn);

			while(ms.scanNumber[0]!=scNum){
				#ifdef _LARGEFILE_SOURCE
				fseeko(fileIn,ms.numZStates*12,1);
				if(compressMe){
					fread(&mzLen,sizeof(uLong),1,fileIn);
					fread(&intensityLen,sizeof(uLong),1,fileIn);
					fseeko(fileIn,mzLen+intensityLen,1);
				} else {	
					fseeko(fileIn,ms.numDataPoints*12,1);
				};
				#else
				fgetpos(fileIn,&fpoint);
				fpoint+=ms.numZStates*12;
				if(compressMe){
					fsetpos(fileIn,&fpoint);
					fread(&mzLen,sizeof(uLong),1,fileIn);
					fread(&intensityLen,sizeof(uLong),1,fileIn);
					fgetpos(fileIn,&fpoint);
					fpoint+=(mzLen+intensityLen);
				} else {
					fpoint+=ms.numDataPoints*12;
				};
				fsetpos(fileIn,&fpoint);
				#endif

				fread(&ms,sizeof(MSScanInfo),1,fileIn);
				if(feof(fileIn)) return false;
			};
		};
		if(feof(fileIn)) return false;

		//read any charge states (for MS2 files)
		for(i=0;i<ms.numZStates;i++){
			fread(&z.z,4,1,fileIn);
			fread(&z.mz,8,1,fileIn);
			s.addZState(z);
		};

		s.setScanNumber(ms.scanNumber[0]);
		s.setRTime(ms.rTime);

		//read compressed data to the spectrum object
		if(compressMe) {
			
			readCompressSpec(fileIn,ms,s);

		//or read binary data to the spectrum object
		} else {	
			for(i=0;i<ms.numDataPoints;i++){
				fread(&p.mz,8,1,fileIn);
				fread(&p.intensity,4,1,fileIn);
				s.add(p);
			};
		};

		//return success
		return true;

	} else {

		//if reading text files, some parsing is required.
		while(true){
			if(feof(fileIn)) break;
    
			//scan next character in the file
			ch=fgetc(fileIn);
			ungetc(ch,fileIn);
   
			switch(ch){
			case 'D':
				//D lines are ignored
				fgets(tstr,256,fileIn);
				break;

			case 'H':
				//Header lines are recorded as strings up to 16 lines at 256 characters each
				fgets(tstr,256,fileIn);
				if(!bDoneHeader) {
					tok=strtok(tstr," \t\n\r");
					tok=strtok(NULL,"\n\r");
					strcat(tok,"\n");
					if(headerIndex<16) strcpy(header.header[headerIndex++],tok);
					else cout << "Header too big!!" << endl;
				};
				break;

			case 'I':
				//I lines are recorded only if they contain retention times
				fgets(tstr,256,fileIn);
				tok=strtok(tstr," \t\n\r");
				tok=strtok(NULL," \t\n\r");
				if(strcmp(tok,"RTime")==0) {
					tok=strtok(NULL," \t\n\r,");
					s.setRTime(atof(tok));
				};
				break;

			case 'S':
				//Scan numbers are recorded and mark all following data is spectrum data
				//until the next tag

				//Reaching an S tag also indicates there are no more header lines
				bDoneHeader=true;

				if(firstScan) {
					//if we are here, a scan was read and we just reached the next scan tag
					//therefore, stop reading further.
					return true;

				} else {
					fgets(tstr,256,fileIn);
					tok=strtok(tstr," \t\n\r");
					tok=strtok(NULL," \t\n\r");
					tok=strtok(NULL," \t\n\r");
					s.setScanNumber(atoi(tok));
					tok=strtok(NULL," \t\n\r");
					if(tok!=NULL)	s.setMZ(atof(tok));
					if(scNum != 0){
						if(s.getScanNumber() != scNum) {
							if(s.getScanNumber()<scNum) bScan=findSpectrum(1);
							else bScan=findSpectrum(-1);
							s.setScanNumber(0);
							if(bScan==false) return false;
							break;
						};
					};
					firstScan=true;
				};
				break;

			case 'Z':
				//Z lines are recorded for MS2 files
				fgets(tstr,256,fileIn);
				tok=strtok(tstr," \t\n\r");
				tok=strtok(NULL," \t\n\r");
				z.z=atoi(tok);
				tok=strtok(NULL," \t\n\r");
				z.mz=atof(tok);
				s.addZState(z);
				break;

			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				//lines beginning with numbers are data; if they belong to a scan we are not
				//interested in, we ignore them.
				if(scNum != 0){
					if(s.getScanNumber()!=scNum) {
						fgets(tstr,256,fileIn);
						break;
					};
				};
				//otherwise, read in the line
				fscanf(fileIn,"%lf %f\n",&p.mz,&p.intensity);
				s.add(p);
				break;
			
			default:
				//if the character is not recognized, ignore the entire line.
				fscanf(fileIn,"%s\n",tstr);
				break;
			};
		};

	};
  
  return true;

};


bool MSReader::findSpectrum(int i){

  if(i==0){
    lPivot = lEnd/2;
    lFWidth = lPivot/2;
  } else if(i==-1){
    lPivot -= lFWidth;
    lFWidth /= 2;
  } else {
    lPivot += lFWidth;
    lFWidth /= 2;
  };
  
	#ifdef _LARGEFILE_SOURCE
  fseeko(fileIn,lPivot,0);
	#else
	fsetpos(fileIn,&lPivot);
	#endif
  return (lFWidth>0 && lPivot>0 && lPivot<lEnd);

};


int MSReader::getPercent(){
	f_off pos;
  if(fileIn!=NULL){
		#ifdef _LARGEFILE_SOURCE
    return (int)((double)ftello(fileIn)/lEnd*100);
		#else
		fgetpos(fileIn,&pos);
		return (int)((double)pos/lEnd*100);
		#endif
  };
	if(rampFileIn!=NULL){
		return (int)((double)rampIndex/rampLastScan*100);
	};
  return -1;
};

void MSReader::writeFile(char* c, bool text, MSObject& m){

  FILE* fileOut;
  int i;

  //if a filename isn't specified, check to see if the
  //MSObject has a filename.
  if(c == NULL) {
		return;
  } else {
		if(text) fileOut=fopen(c,"wt");
		else fileOut=fopen(c,"wb");
	};

  //output file header lines;
	if(text){
		for(i=0;i<16;i++){
			if(m.getHeader().header[i][0]!='\0') {
				fputs("H\t",fileOut);
				fputs(m.getHeader().header[i],fileOut);
			};
		};
	} else {
		fwrite(&m.getHeader(),sizeof(MSHeader),1,fileOut);
	};
  
	//output spectra;
  for(i=0;i<m.size();i++){

		//output spectrum header
		writeSpecHeader(fileOut,text,m.at(i));

		//output scan
		if(text){
			writeTextSpec(fileOut,m.at(i));
		} else if(compressMe){
			writeCompressSpec(fileOut,m.at(i));
		} else {
			writeBinarySpec(fileOut,m.at(i));
		};

  };
    
	fclose(fileOut);
};

void MSReader::appendFile(char* c, bool text, Spectrum& s){
	FILE* fileOut;

	if(c == NULL) return;
	
	if(text)fileOut=fopen(c,"at");
	else fileOut=fopen(c,"ab");

	//output spectrum header
	writeSpecHeader(fileOut,text,s);

  //output spectrum
	if(text){
		writeTextSpec(fileOut,s);
	} else if(compressMe){
		writeCompressSpec(fileOut,s);
	} else {
		writeBinarySpec(fileOut,s);
	};
    
	fclose(fileOut);


};

void MSReader::appendFile(char* c, bool text, MSObject& m){

  FILE* fileOut;
  int i;

  //if a filename isn't specified, check to see if the
  //MSObject has a filename.
  if(c == NULL) {
		return;
  } else {
		if(text) fileOut=fopen(c,"at");
		else fileOut=fopen(c,"ab");
	};

  //output spectra;
  for(i=0;i<m.size();i++){

		//output spectrum header
		writeSpecHeader(fileOut,text,m.at(i));

		//output spectrum
		if(text){
			writeTextSpec(fileOut,m.at(i));
		} else if(compressMe){
			writeCompressSpec(fileOut,m.at(i));
		} else {
			writeBinarySpec(fileOut,m.at(i));
		};

	};
    
	fclose(fileOut);
};

void MSReader::setPrecision(int i, int j){
	iIntensityPrecision=i;
	iMZPrecision=j;
};

void MSReader::setPrecisionInt(int i){
	iIntensityPrecision=i;
};

void MSReader::setPrecisionMZ(int i){
	iMZPrecision=i;
};

bool MSReader::readFile(char* c, MSFileFormat f, Spectrum& s, int scNum){

	//Redirect functions to appropriate places, if possible.
	switch(f){
		case ms1:
		case ms2:
			return readFile(c,true,s,scNum);
			break;
		case bms1:
		case bms2:
			setCompression(false);
			return readFile(c,false,s,scNum);
			break;
		case cms1:
		case cms2:
			setCompression(true);
			return readFile(c,false,s,scNum);
			break;
		case mzXML:
			break;
		default:
			return false;
			break;
	};

	//if we got here, it's because we're reading mzXML format

	ramp_fileoffset_t indexOffset;
	ScanHeaderStruct scanHeader;
	RAMPREAL *pPeaks;
	int i,j;
	
	if(c!=NULL) {
		//open the file if new file was requested
		if(rampFileOpen) {
			rampCloseFile(rampFileIn);
			rampFileOpen=false;
			free(pScanIndex);
		};
		rampFileIn = rampOpenFile(c);
		if (rampFileIn == NULL) {
      cout << "Error reading input file " << c << endl;
      return false;
		};
		rampFileOpen=true;
		//read the index
		indexOffset = getIndexOffset(rampFileIn);
		pScanIndex = readIndex(rampFileIn,indexOffset,&rampLastScan);
		rampIndex=0;

	} else {
		//if no new file requested, check to see if one is open already
		if (rampFileIn == NULL) return false;
	};


	//clear any spectrum data
	s.clear();

	//read scan header
	if(scNum!=0) {
		rampIndex=0;
		for(i=1;i<rampLastScan;i++){
			readHeader(rampFileIn, pScanIndex[i], &scanHeader);
			if(scanHeader.acquisitionNum==scNum) {
				rampIndex=i;
				break;
			};
		};
		if(rampIndex==0) return false;

		readHeader(rampFileIn, pScanIndex[rampIndex], &scanHeader);
		switch(filter){
		case MS1:
			if(scanHeader.msLevel!=1)	return false;
			break;
		case MS2:
			if(scanHeader.msLevel!=2)	return false;
			break;
		case MS3:
			if(scanHeader.msLevel!=3)	return false;
			break;
		default:
			//no filter
			break;
		};
		s.setScanNumber(scanHeader.acquisitionNum);
		s.setRTime((float)scanHeader.retentionTime);
		pPeaks = readPeaks(rampFileIn, pScanIndex[rampIndex]);
		j=0;
		for(i=0;i<scanHeader.peaksCount;i++){
			s.add((double)pPeaks[j],(float)pPeaks[j+1]);
			j+=2;
		};

	} else {
		
		//read next index
		while(true){
			rampIndex++;

			//reached end of file
			if(rampIndex>rampLastScan) {
				//rampCloseFile(rampFileIn);
				//rampFileIn = NULL;
				return false;
			};
			readHeader(rampFileIn, pScanIndex[rampIndex], &scanHeader);
			
			switch(filter){
			case MS1:
				if(scanHeader.msLevel!=1)	continue;
				break;
			case MS2:
				if(scanHeader.msLevel!=2)	continue;
				break;
			case MS3:
				if(scanHeader.msLevel!=3)	continue;
				break;
			default:
				//no filter
				break;
			};

			//if we got here, we passed the filter.
			break;
		};

		s.setScanNumber(scanHeader.acquisitionNum);
		s.setRTime((float)scanHeader.retentionTime);
		pPeaks = readPeaks(rampFileIn, pScanIndex[rampIndex]);
		j=0;
		for(i=0;i<scanHeader.peaksCount;i++){
			s.add((double)pPeaks[j],(float)pPeaks[j+1]);
			j+=2;
		};

	};

	free(pPeaks);
	return true;

};

void MSReader::setFilter(MSFileType m){
	filter=m;
};

void MSReader::setCompression(bool b){
	compressMe=b;
};

void MSReader::writeCompressSpec(FILE* fileOut, Spectrum& s){

	int j;

	//file compression
	int err;
	uLong len;
	Byte *comprM, *comprI;
  uLong comprLenM, comprLenI;
	double *pD;
	float *pF;
	uLong sizeM;
	uLong sizeI;

	//Build arrays to hold scan prior to compression
	// Ideally, we would just use the scan vectors, but I don't know how yet.
	pD = new double[s.size()];
	pF = new float[s.size()];
	for(j=0;j<s.size();j++){
		pD[j]=s.at(j).mz;
		pF[j]=s.at(j).intensity;
	};

	//compress mz
	len = (uLong)s.size()*sizeof(double);
	sizeM = len;
	comprLenM = compressBound(len);
	comprM = (Byte*)calloc((uInt)comprLenM, 1);
	err = compress(comprM, &comprLenM, (const Bytef*)pD, len);

	//compress intensity
	len = (uLong)s.size()*sizeof(float);
	sizeI = len;
	comprLenI = compressBound(len);
	comprI = (Byte*)calloc((uInt)comprLenI, 1);
	err = compress(comprI, &comprLenI, (const Bytef*)pF, len);

	fwrite(&comprLenM,sizeof(uLong),1,fileOut);
	fwrite(&comprLenI,sizeof(uLong),1,fileOut);
	fwrite(comprM,comprLenM,1,fileOut);
	fwrite(comprI,comprLenI,1,fileOut);

	//clean up memory
	free(comprM);
	free(comprI);
	delete [] pD;
	delete [] pF;

};

void MSReader::readCompressSpec(FILE* fileIn, MSScanInfo& ms, Spectrum& s){

	int i;
	Peak_T p;

	//variables for compressed files
	uLong uncomprLen;
	uLong mzLen, intensityLen;
	Byte *compr;
	double *mz;
	float *intensity;

	fread(&mzLen,sizeof(uLong),1,fileIn);
	fread(&intensityLen,sizeof(uLong),1,fileIn);
		
	compr = new Byte[mzLen];
	mz = new double[ms.numDataPoints];
	uncomprLen=ms.numDataPoints*sizeof(double);
	fread(compr,mzLen,1,fileIn);
	uncompress((Bytef*)mz, &uncomprLen, compr, mzLen);
	delete [] compr;
			
	compr = new Byte[intensityLen];
	intensity = new float[ms.numDataPoints];
	uncomprLen=ms.numDataPoints*sizeof(float);
	fread(compr,intensityLen,1,fileIn);
	uncompress((Bytef*)intensity, &uncomprLen, compr, intensityLen);
	delete [] compr;
		
	for(i=0;i<ms.numDataPoints;i++){
		p.mz = mz[i];
		p.intensity = intensity[i];
		s.add(p);
	};			

	delete [] mz;
	delete [] intensity;

};

void MSReader::writeTextSpec(FILE* fileOut, Spectrum& s) {

	int j,k;
	char t[64];

	for(j=0;j<s.size();j++){
		sprintf(t,"%.*f",iIntensityPrecision,s.at(j).intensity);
		k=strlen(t);
		if(k>2 && iIntensityPrecision>0){
			if(t[0]=='0'){
				fprintf(fileOut,"%.*f 0\n",iMZPrecision,s.at(j).mz);
			} else if(t[k-1]=='0'){
				fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision-1,s.at(j).intensity);
			} else {
				fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
			};
		} else {
			fprintf(fileOut,"%.*f %.*f\n",iMZPrecision,s.at(j).mz,iIntensityPrecision,s.at(j).intensity);
		};
	};
 
};

void MSReader::writeBinarySpec(FILE* fileOut, Spectrum& s) {
	int j;

	for(j=0;j<s.size();j++){
		fwrite(&s.at(j).mz,8,1,fileOut);
		fwrite(&s.at(j).intensity,4,1,fileOut);
	};

};

void MSReader::writeSpecHeader(FILE* fileOut, bool text, Spectrum& s) {

	MSScanInfo ms;
	MSFileType mft;
	int j;

	//output scan info
	if(text){
		mft=s.getFileType();
		if(mft==MS2 || mft==MS3 || mft==SRM){
			fprintf(fileOut,"S\t%d\t%d\t%.*f\n",s.getScanNumber(),s.getScanNumber(),2,s.getMZ());
		} else {
			fprintf(fileOut,"S\t%d\t%d\n",s.getScanNumber(),s.getScanNumber());
		};
		if(s.getRTime()>0) fprintf(fileOut,"I\tRTime\t%.*f\n",4,s.getRTime());
		for(j=0;j<s.sizeZ();j++){
			fprintf(fileOut,"Z\t%d\t%.*f\n",s.atZ(j).z,2,s.atZ(j).mz);
		};
	} else {
		ms.scanNumber[0]=ms.scanNumber[1]=s.getScanNumber();
		ms.rTime=s.getRTime();
		ms.numDataPoints=s.size();
		ms.numZStates=s.sizeZ();
		fwrite(&ms,sizeof(MSScanInfo),1,fileOut);
		for(j=0;j<ms.numZStates;j++){
			fwrite(&s.atZ(j).z,4,1,fileOut);
			fwrite(&s.atZ(j).mz,8,1,fileOut);
		};
	};

};


