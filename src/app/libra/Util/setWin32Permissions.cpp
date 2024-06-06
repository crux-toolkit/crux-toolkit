//
// DESCRIPTION:
//
// tiny wrapper for setACL.exe for use in TPP win32 (MinGW or VC8) installer
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
// General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// NOTES:
//
//
// TODO:
//

#include <stdlib.h>
#include <stdio.h>

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't bark about "unsafe" functions
#endif

int main(int argc,char *argv[]) {
   char buf[4096];
   printf("TPP Installer:  setting file permissions in %s... this may take a few minutes, please do not close this window.\n",argv[1]);
   // cut loose from parent permissions ( -actn setprot -op "dacl:p_nc;sacl:p_nc")
   sprintf(buf,"setACL -on \"%s\" -ot file -actn clear -clr \"dacl\" -actn rstchldrn -rst \"dacl\"",argv[1]);
   // set generous permissions for:
   // BUILTIN\ADMINISTRATORS S-1-5-32-544
   // BUILTIN\USERS SID=S-1-5-32-545
   // BUILTIN\GUESTS SID=S-1-5-32-546
   // Everyone SID=S-1-1-0
   // Anonymous SID=S-1-5-7
   sprintf(buf,"setACL -on \"%s\" -ot file -actn ace -ace \"n:S-1-5-32-544;p:change,read_ex,write,full;s:y\" -ace \"n:S-1-1-0;p:change,read_ex,write,full;s:y\" -ace \"n:S-1-5-32-545;p:change,read_ex,write,full;s:y\" -ace \"n:S-1-5-32-546;p:change,read_ex,write,full;s:y\" -ace \"n:S-1-5-7;p:change,read_ex,write,full;s:y\"",argv[1]);

   return system(buf);
}
