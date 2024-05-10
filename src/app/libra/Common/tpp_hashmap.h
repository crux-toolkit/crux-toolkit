//
// tpp_hashmap.h
//
// help with portability issues around STL hash_map
//
// Copyright (c) 2006, 2007 Insilicos LLC and LabKey Software. All rights reserved.
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
// DESCRIPTION:
//
// NOTES:
//
// SVN INFO:
// $Id$
//
//////////////////////////////////////////////////////////////////////

#pragma warning(disable: 4786)

#if !defined(AFX_TPP_HASHMAP_H__A003F5B4_8E62_49E2_BBA4_F7F3FA90E2B2__INCLUDED_)
#define AFX_TPP_HASHMAP_H__A003F5B4_8E62_49E2_BBA4_F7F3FA90E2B2__INCLUDED_

#ifdef _MSC_VER
#include <hash_map> // if you're using VC6 you need STLPort 
#ifdef _STLP_HASH_MAP  // STLPort
// a bit of ugly cut and paste to deal with VC template borkeness
inline size_t __stl_hash_string(const char* __s)
{
  _STLP_FIX_LITERAL_BUG(__s)
    unsigned long __h = 0; 
  for ( ; *__s; ++__s)
    __h = 5*__h + *__s;

  return size_t(__h);
}
struct cchash
{
  size_t operator()(const char* __s) const { _STLP_FIX_LITERAL_BUG(__s) return __stl_hash_string(__s); }
};
struct stdstringhash
{
  size_t operator()(const std::string & __s) const { _STLP_FIX_LITERAL_BUG(__s) return __stl_hash_string(__s.c_str()); }
};
#define TPP_HASH_MULTIMAP_T std::hash_multimap
#define TPP_HASHMAP_T std::hash_map
#define TPP_HASH_T std::hash
#else // assume VC8 STL
#define TPP_HASHMAP_VC8
#define TPP_HASH_MULTIMAP_T stdext::hash_multimap
#define TPP_HASHMAP_T stdext::hash_map
#define TPP_HASH_T stdext::hash
#endif
#else // GCC
#if (__GNUC__ < 4)
#include <hash_map.h>
#define TPP_HASH_MULTIMAP_T hash_multimap
#define TPP_HASHMAP_T hash_map
#define TPP_HASH_T hash
#else
#include <ext/hash_map>
#define TPP_HASH_MULTIMAP_T __gnu_cxx::hash_multimap
#define TPP_HASHMAP_T __gnu_cxx::hash_map
#define TPP_HASH_T __gnu_cxx::hash
#endif
#define cchash TPP_HASH_T<const char *> 
struct stdstringhash {
  size_t operator()(const std::string& s) const {
    cchash h;
    return h(s.c_str());
  }
};
#endif
#include <string.h>
struct equalstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2)< 0;
  }
};
struct equalstdstring
{
  bool operator()(const std::string &s1, const std::string &s2) const
  {
    return strcmp(s1.c_str(), s2.c_str()) == 0;
  }
};

struct ltstdstring
{
  bool operator()(const std::string &s1, const std::string &s2) const
  {
    return strcmp(s1.c_str(), s2.c_str()) < 0;
  }
};
#ifdef _DEBUG_SORT // alpha sorted may be easier for debug, but it's slow
#define TPP_CONSTCHARP_HASHMAP(T) std::map<const char *, T,  ltstr>
#define TPP_STDSTRING_HASHMAP(T) std::map<std::string, T,  ltsrdstring>
#define TPP_STDSTRING_HASH_MULTIMAP(T) std::multimap<std::string, T,  ltsrdstring>
#else
#ifdef TPP_HASHMAP_VC8
#define TPP_CONSTCHARP_HASHMAP(T) TPP_HASHMAP_T<const char *, T, stdext::hash_compare<const char *,ltstr>>
#define TPP_STDSTRING_HASHMAP(T) TPP_HASHMAP_T<std::string, T, stdext::hash_compare<std::string,ltstdstring>>
#define TPP_STDSTRING_HASH_MULTIMAP(T) TPP_HASH_MULTIMAP_T<std::string, T, stdext::hash_compare<std::string,ltstdstring>>
#else
#define TPP_CONSTCHARP_HASHMAP(T) TPP_HASHMAP_T<const char *, T, cchash, equalstr> 
#define TPP_STDSTRING_HASHMAP(T) TPP_HASHMAP_T<std::string, T, stdstringhash, equalstdstring> 
#define TPP_STDSTRING_HASH_MULTIMAP(T) TPP_HASHMAP_T<std::string, T, stdstringhash, equalstdstring> 
#endif
#endif

#endif // AFX_TPP_HASHMAP_H__A003F5B4_8E62_49E2_BBA4_F7F3FA90E2B2__INCLUDED_
