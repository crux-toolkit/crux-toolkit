// cached_ramp.h
//
// DESCRIPTION:
//
// uses #define to redirect RAMP API to a smarter mzdata/mzxml reader - 
// for use in TPP where the authors tend to reopen the same file a lot, so we cache
//
// Copyright (c) 2006,2007 Insilicos, LLC
//
// This library is free software; you can redistribute it and/or 
// modify it under the terms of the GNU Lesser General Public 
// License as published by the Free Software Foundation; either 
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public 
// License along with this library; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
// 
// Brian Pratt 
// Insilicos LLC 
// www.insilicos.com
//
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CACHED_RAMP__H__03CA0A06_4827_4CC4_BE2B_72178B5D0E71__INCLUDED_)
#define AFX_CACHED_RAMP__H__03CA0A06_4827_4CC4_BE2B_72178B5D0E71__INCLUDED_

#include "ramp.h" // the real RAMP API - make sure it's included before our #defines

//
// RAMP-like function set that caches file opens and file reads
//
RAMPFILE *cached_ramp_rampOpenFile(const char *filename);
void cached_ramp_rampCloseFile(RAMPFILE *pFI);
ramp_fileoffset_t cached_ramp_getIndexOffset(RAMPFILE *pFI);
ramp_fileoffset_t *cached_ramp_readIndex(RAMPFILE *pFI,
                ramp_fileoffset_t indexOffset,
                int *iLastScan);
void cached_ramp_readHeader(RAMPFILE *pFI,
                ramp_fileoffset_t lScanIndex, // read from this file position
                struct ScanHeaderStruct *scanHeader);
const struct ScanHeaderStruct *cached_ramp_readHeader(RAMPFILE *pFI,
                ramp_fileoffset_t lScanIndex); // read from this file position
int  cached_ramp_readMsLevel(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex);
double cached_ramp_readStartMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex);
double cached_ramp_readEndMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex);
int cached_ramp_readPeaksCount(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex);
RAMPREAL *cached_ramp_readPeaks(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex);
const RAMPREAL *cached_ramp_readPeaks_const(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex);
void cached_ramp_readRunHeader(RAMPFILE *pFI,
                   ramp_fileoffset_t *pScanIndex,
                   struct RunHeaderStruct *runHeader,
                   int iLastScan);
void cached_ramp_readMSRun(RAMPFILE *pFI,
                   struct RunHeaderStruct *runHeader);

#ifndef CACHED_RAMP_HOME
//
// redirect the RAMP API to our cacheing code with minimal calling source code impact
//
#define rampOpenFile cached_ramp_rampOpenFile
#define rampCloseFile cached_ramp_rampCloseFile
#define getIndexOffset cached_ramp_getIndexOffset
#define readIndex cached_ramp_readIndex
#define readHeader cached_ramp_readHeader
#define readMsLevel cached_ramp_readMsLevel
#define readStartMz cached_ramp_readStartMz
#define readEndMz cached_ramp_readEndMz
#define readPeaksCount cached_ramp_readPeaksCount
#define readPeaks cached_ramp_readPeaks
#define readRunHeader cached_ramp_readRunHeader
#define readMSRun cached_ramp_readMSRun
#endif // CACHED_RAMP_HOME



#endif // !defined(AFX_CACHED_RAMP__H__03CA0A06_4827_4CC4_BE2B_72178B5D0E71__INCLUDED_)
