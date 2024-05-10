/*

Copyright (C) 2014 Institute for Systems Biology

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

Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA

*/

/*
   Version and build information are compiled into TPP programs using
   this source file.  The version information is expected to be defined
   as TPP_VERSIONINFO and should contain the version and any additional
   build indentification.

   Hint, use -DTPP_VERSIONINFO="" to define when compiling.
*/
#include "TPPVersion.h"

#ifndef TPP_VERSIONINFO
#warning "TPP_VERSIONINFO was not defined"
#define TPP_VERSIONINFO_STR "TPP version information unavailable"
#else
#define TPP_VERSIONINFO_STR TPP_VERSIONINFO
#endif

const char *szTPPVersionInfo = TPP_VERSIONINFO_STR;
