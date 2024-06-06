// mzp_cached_ramp.cpp
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

#define CACHED_RAMP_HOME // code lives here
#include "Parsers/mzParser/cached_ramp.h"
#ifdef _RAMP_H
#error "wrong interface - got ramp.h, wanted mzParser.h"
#endif
using namespace std;
using namespace mzParser;
#include "../ramp/cached_ramp.cpp"
