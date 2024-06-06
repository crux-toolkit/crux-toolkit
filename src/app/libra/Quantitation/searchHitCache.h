/*

Program       : searchHitCache.h                                                 

cacheing for programs that read peptide info from pepXML files

Copyright (C) 2007 Labkey

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/

#ifndef SEARCHHITCACHE_H
#define SEARCHHITCACHE_H

#include "Parsers/Parser/Tag.h"
typedef charptr_int_map peplist; // borrowing the string+index collection class from Tag
// create a multimap to locate peptide in stored search hits
// we'll cast those void ptrs to type T in searchHitCache
#ifdef _MSC_VER // assuming VC8 STL
#include <hash_map> // if you're using VC6 you need STLPort 
#ifdef _STLP_HASH_MAP 
typedef stdext::hash_multimap<const char *,void *, std::hash<const char*>, tageqstr> searchHitMap;
#else // VC8/Dinkumware
typedef stdext::hash_multimap<const char *,void *, stdext::hash_compare<const char*, tagltstr>> searchHitMap;
#endif
#else
#include <ext/hash_map>
using namespace __gnu_cxx;
typedef hash_multimap<const char *,void *, __gnu_cxx::hash<const char*>, tageqstr> searchHitMap;
#endif
typedef std::pair<searchHitMap::const_iterator,searchHitMap::const_iterator> hitRange;

// helper class for cacheing pepXML reads
template <class T> class searchHitCache {
public:
   searchHitCache() {
      bReady = false;
   };
   ~searchHitCache() {
      std::vector<char *>strs_to_free; // can't free the keys while map exists
      for (searchHitMap::iterator iter = hits_.begin();iter!=hits_.end();iter++) {
         strs_to_free.push_back((char *)(iter->first)); // we allocated it, so cast away const
         delete ((T *)iter->second); // we allocated this too
      }
      for (int i = (int)strs_to_free.size();i--;) {
         free(strs_to_free[i]);
      }
   }

   void addHit(const char *peptideName,const T &peptideInfo) {
      T* p = new T(peptideInfo);
      hits_.insert(std::pair<const char *,void *>(strdup(peptideName),(void *)p));
   }
         
   // return a sorted vector of entries that match the listed peptides
   std::vector<const T*> getMatches(const peplist &peptides) const {
      std::vector<const T*> hits;
      for (int p=0;p<peptides.size();p++) {
         hitRange range = 
            hits_.equal_range(peptides.getPair(p)->charptr_);
         while (range.first != range.second) {
            hits.push_back((const T *)(range.first->second));
            range.first++;
         }
      }
      return hits;
   }
   void setReady() {
      bReady = true; // if it's empty, its just that pepxml was uninteresting
   }
   bool isReady() const {
      return bReady; // if it's empty, its just that pepxml was uninteresting
   }

private:
   bool bReady;
   searchHitMap hits_; // note (not necessarily unique) peptide name and associated info

};


#endif


