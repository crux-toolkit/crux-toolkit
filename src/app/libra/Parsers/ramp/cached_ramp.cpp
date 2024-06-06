// cached_ramp.cpp
//
// DESCRIPTION:
//
// smarter mzdata/mzxml reader - for use in TPP where the authors tend 
// to reopen the same file a lot, so we cache
//
// NOTES:
//
// this class is meant to be a singleton, so we're just going to go 
// with static globals instead of junking up the class declaration
//
// VERSION:
//
// $Id: cached_ramp.cpp 8868 2023-02-24 18:46:05Z real_procopio $
//
// Copyright (c) 2006 Insilicos, LLC
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

#include "string.h"
#include <string>
#include <stdio.h>
#include "Common/sysdepend.h"
#define CACHED_RAMP_HOME // code lives here
#include "cached_ramp.h"
#include <map> // if you're using VC6 you need STLPort 
#include "Common/Array.h"
#include "Common/util.h"

class cachedscan {
public:
   cachedscan() : m_data(NULL),m_index(0) {}
   ~cachedscan() {
     free(m_data); 
   }

  RAMPREAL *m_data;
  int peakCapacity;
  ramp_fileoffset_t m_index;
};

class mzCache {
public:
   mzCache() {
      m_filePtr = NULL;
      m_filename = NULL;
      m_scanIndex = NULL;
      m_indexOffset = 0;
      m_scanCount = 0; // 1 .. totScan
      memset(&m_runHeader,0,sizeof(RunHeaderStruct));
//#define GETSTATS
#ifdef GETSTATS // for cache miss investigations
      m_nMisses = 0;
      m_nReads = 0;
      m_nCachedReads = 0;
#endif
      // prepare the cache
      m_cachedScans.setSize(500); // cache the last n spectra read
      m_cachedScans.nullify();
   };
   ~mzCache() {
#ifdef GETSTATS // for cache miss investigations
      printf("%s %d reads, %d missed, %d saved\n",m_filename,m_nReads,m_nMisses,m_nCachedReads);
#endif
      free(m_scanIndex);
      rampCloseFile(m_filePtr);
      delete[] m_filename;
      for (int i=m_cachedScans.size();i--;) {
         delete m_cachedScans[i];
      }
   }
   RAMPFILE *m_filePtr;
   char *m_filename;
   ramp_fileoffset_t *m_scanIndex;
   ramp_fileoffset_t m_indexOffset;
   RunHeaderStruct m_runHeader;
   std::map<ramp_fileoffset_t, struct ScanHeaderStruct *> m_scanHeaders;
   Array<cachedscan *> m_cachedScans; // cache the last n spectra read
   int m_scanCount; // 1 .. totScan
#ifdef GETSTATS // for cache miss investigations
   Array<unsigned long> m_scansRead; 
   int m_nCachedReads;
   int m_nReads;
   int m_nMisses;
#endif
};

typedef std::map<RAMPFILE *,mzCache *> mzmap_t;

static mzmap_t *the_map; // just one per app
static int n_constructions = 0;


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

static void cleanup_cached_ramp() {
  if (n_constructions) {
    mzmap_t::const_iterator iter;
    for (iter=the_map->begin(); iter != the_map->end(); ++iter) {
      delete(iter->second);
    }
    delete the_map;
    n_constructions = 0;
  }
}

static void init_cached_ramp() {
  if (!n_constructions) { // this should never exceed 1
    // redirect ramp users to our function table
    the_map = new mzmap_t;
    atexit(cleanup_cached_ramp);
    n_constructions++;
  }
}

RAMPFILE *cached_ramp_rampOpenFile(const char *filename) {
    init_cached_ramp(); // in case this is first access
    mzmap_t::const_iterator iter;
    for (iter=the_map->begin(); iter != the_map->end(); ++iter) {
       if (!strcmp(filename,(iter->second->m_filename))) {
          return iter->first; // already open
       }
    }
    RAMPFILE *pFI = rampOpenFile(filename);
    mzCache *cache = new mzCache; 
    (*the_map)[pFI] = cache;
    cache->m_filename = strCopy(filename);
    cache->m_filePtr = pFI;
    return pFI;
}

void cached_ramp_rampCloseFile(RAMPFILE *pFI) {
  // don't actually close it - they'll probably just reopen, and make us
  // re-read the index - but do let go of the scan and header caches
  mzmap_t::iterator iter;
  for (iter=the_map->begin(); iter != the_map->end(); ++iter) {
    if (pFI == iter->first) { // found our entry
      for (int i=iter->second->m_cachedScans.size();i--;) {
	delete iter->second->m_cachedScans[i];
	iter->second->m_cachedScans[i] = NULL;
      }
      if (!iter->second->m_scanHeaders.empty()) {
	for (map <ramp_fileoffset_t, struct ScanHeaderStruct *>::const_iterator it = iter->second->m_scanHeaders.begin(); it != iter->second->m_scanHeaders.end(); it++)
	  delete it->second;
      }
      break;
    }
  }
}

ramp_fileoffset_t cached_ramp_getIndexOffset(RAMPFILE *pFI) {
   mzCache *cache = (*the_map)[pFI];
   if (!cache->m_indexOffset) {
      cache->m_indexOffset = getIndexOffset(pFI);
   }
   return cache->m_indexOffset;
}

ramp_fileoffset_t *cached_ramp_readIndex(RAMPFILE *pFI,
                ramp_fileoffset_t indexOffset,
                int *iLastScan) {
   mzCache *cache = (*the_map)[pFI];
   if (!cache->m_scanIndex) { // nab the index
      cache->m_scanIndex = readIndex(pFI,indexOffset,iLastScan);
      cache->m_scanCount = *iLastScan;
   } else {
      *iLastScan = cache->m_scanCount;
   }
   // now copy the index
   ramp_fileoffset_t *result = (ramp_fileoffset_t *)malloc((cache->m_scanCount+1)*sizeof(ramp_fileoffset_t));
   memcpy(result,cache->m_scanIndex,(cache->m_scanCount+1)*sizeof(ramp_fileoffset_t));
   return result;
}

// local helper func
const struct ScanHeaderStruct *cached_ramp_getHeader(RAMPFILE *pFI,
                                          ramp_fileoffset_t lScanIndex) {
   mzCache *cache = (*the_map)[pFI];
   struct ScanHeaderStruct *hdr = cache->m_scanHeaders[lScanIndex];
   if (!hdr) { // first access?
      hdr = cache->m_scanHeaders[lScanIndex] = new(struct ScanHeaderStruct);
      readHeader(pFI,lScanIndex, hdr); 
   }
   return hdr;
}

void cached_ramp_readHeader(RAMPFILE *pFI,
                ramp_fileoffset_t lScanIndex, // read from this file position
                struct ScanHeaderStruct *scanHeader) {
   const struct ScanHeaderStruct *hdr = cached_ramp_getHeader(pFI,lScanIndex);
   memcpy(scanHeader,hdr,sizeof(struct ScanHeaderStruct));
}

const struct ScanHeaderStruct *cached_ramp_readHeader(RAMPFILE *pFI,
	      ramp_fileoffset_t lScanIndex) { // read from this file position
   return cached_ramp_getHeader(pFI,lScanIndex);
}

int  cached_ramp_readMsLevel(RAMPFILE *pFI,
                           ramp_fileoffset_t lScanIndex) {
   const struct ScanHeaderStruct *hdr = cached_ramp_getHeader(pFI,lScanIndex);
   return hdr->msLevel;
}
double cached_ramp_readStartMz(RAMPFILE *pFI,
                             ramp_fileoffset_t lScanIndex) {
   const struct ScanHeaderStruct *hdr = cached_ramp_getHeader(pFI,lScanIndex);
   return hdr->lowMZ;
}
double cached_ramp_readEndMz(RAMPFILE *pFI,
		   ramp_fileoffset_t lScanIndex) {
   const struct ScanHeaderStruct *hdr = cached_ramp_getHeader(pFI,lScanIndex);
   return hdr->highMZ;
}
int cached_ramp_readPeaksCount(RAMPFILE *pFI,
                 ramp_fileoffset_t lScanIndex) {
   const struct ScanHeaderStruct *hdr = cached_ramp_getHeader(pFI,lScanIndex);
   return hdr->peaksCount;
}
void cached_ramp_readRunHeader(RAMPFILE *pFI,
                   ramp_fileoffset_t *pScanIndex,
                   struct RunHeaderStruct *runHeader,
                   int iLastScan) {
   mzCache *cache = (*the_map)[pFI];
   if (!cache->m_runHeader.scanCount) {
      readRunHeader(pFI,pScanIndex,runHeader,iLastScan);
      memcpy(&(cache->m_runHeader),runHeader,sizeof(struct RunHeaderStruct));
   } else {
      memcpy(runHeader,&(cache->m_runHeader),sizeof(struct RunHeaderStruct));
   }
}
void cached_ramp_readMSRun(RAMPFILE *pFI,
                         struct RunHeaderStruct *runHeader) {
   mzCache *cache = (*the_map)[pFI];
   if (!cache->m_runHeader.highMZ) {
      readMSRun(pFI,runHeader);
      memcpy(&(cache->m_runHeader),runHeader,sizeof(struct RunHeaderStruct));
   } else {
      memcpy(runHeader,&(cache->m_runHeader),sizeof(struct RunHeaderStruct));
   }
}

//IMPORTANT!!! When you ask for a bDeepCopy you must deallocate the memory when you are done otherwise it creates a memory leak!
static RAMPREAL *cached_ramp_readPeaksHelper(RAMPFILE *pFI,
                              ramp_fileoffset_t lScanIndex,
                              bool bDeepCopy) {
   mzCache *cache = (*the_map)[pFI];
#ifdef GETSTATS // for cache miss investigations
   cache->m_nReads++;
#endif
   // do we have this in the cache?
   int slot;
   for (slot=cache->m_cachedScans.size();slot--;) {
      if (!cache->m_cachedScans[slot]) {
         break; // didn't find anything, but we have an available slot
      }
      if (cache->m_cachedScans[slot]->m_index == lScanIndex) {
#ifdef GETSTATS // for cache miss investigations
         cache->m_nCachedReads++;
#endif
         break; // we read that recently
      }
   }
   if ((slot<0) || !cache->m_cachedScans[slot]) { // read it
      if (slot < 0) { // reuse oldest slot (which is slot 0)
        delete cache->m_cachedScans[slot = 0];
      } 
#ifdef GETSTATS // for cache miss investigations
      int m;
      for (m=cache->m_scansRead.size();m--;) {
         if (cache->m_scansRead[m]==lScanIndex) { 
            cache->m_nMisses++;
            break;
         }
      }
      if (m<0) { // first access
         cache->m_scansRead.insertAtEnd(lScanIndex);
      }
#endif
      // read the scan and cache it
      cachedscan *scan = new cachedscan; 
      cache->m_cachedScans[slot] = scan;

      RAMPREAL *tmp = readPeaks(pFI,lScanIndex); //readPeaks no longer allocates memory when not needed must make a copy here!

      if ((pFI->fileType == 3 || pFI->fileType == 1) && pFI->mzML->getIonMobility()) {
	scan->m_data =   (RAMPREAL *) malloc((pFI->bs->size()+1) * 3 * sizeof(RAMPREAL) + 1);
	memcpy(scan->m_data,tmp,(pFI->bs->size()+1)*3*sizeof(RAMPREAL));
      }
      else  {
	scan->m_data =   (RAMPREAL *) malloc((pFI->bs->size()+1) * 2 * sizeof(RAMPREAL) + 1);
      	memcpy(scan->m_data,tmp,(pFI->bs->size()+1)*2*sizeof(RAMPREAL));
      }

       
      scan->m_index = lScanIndex;
   }
   // copy the scan data
   struct ScanHeaderStruct *hdr = cache->m_scanHeaders[lScanIndex];
   RAMPREAL *result;
   
   
   if (bDeepCopy) {
     if ((pFI->fileType == 3 || pFI->fileType == 1) && pFI->mzML->getIonMobility()) {
       result = (RAMPREAL *)malloc((hdr->peaksCount+1)*3*sizeof(RAMPREAL)); // catch that last -1
     }
     else {
       result = (RAMPREAL *)malloc((hdr->peaksCount+1)*2*sizeof(RAMPREAL)); // catch that last -1
     }
     
     if (hdr->peaksCount) {
       if ((pFI->fileType == 3 || pFI->fileType == 1) && pFI->mzML->getIonMobility()) {
	 memcpy(result,cache->m_cachedScans[slot]->m_data,(hdr->peaksCount+1)*3*sizeof(RAMPREAL));
       }
       else {
	 memcpy(result,cache->m_cachedScans[slot]->m_data,(hdr->peaksCount+1)*2*sizeof(RAMPREAL));
       }
     }
     else {
       *result = -1; // empty scan
     }
   }
   else {
     result = cache->m_cachedScans[slot]->m_data;
   }
   
   // now update the cache aging
   cache->m_cachedScans.moveLast(slot);
   return result;
}

RAMPREAL *cached_ramp_readPeaks(RAMPFILE *pFI,
                              ramp_fileoffset_t lScanIndex) {
   return cached_ramp_readPeaksHelper(pFI, lScanIndex, true); // true=make a copy
}

const RAMPREAL *cached_ramp_readPeaks_const(RAMPFILE *pFI,
                              ramp_fileoffset_t lScanIndex) {
   return cached_ramp_readPeaksHelper(pFI, lScanIndex, false); // false=don't make a copy
}




