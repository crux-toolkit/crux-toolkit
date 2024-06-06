/*
Program       : ASAPRatio
Author        : Xiao-jun Li <xli@systemsbiology.org>
Date          : 09.17.02
SVN Info      : $Id: ASAPRatio_txtFns.h 7996 2019-12-25 00:16:42Z real_procopio $

Header file of ASAPRatio_txtFns.c

Copyright (C) 2002 Xiao-jun Li

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

Xiao-jun Li
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
xli@systemsbiology.org
*/

#ifndef _ASAPRATIO_TXTFNS_H_
#define _ASAPRATIO_TXTFNS_H_


#include "Common/sysdepend.h" // defn

//
// Constants
//
#define _MXSTRLEN_ 20000 // maximium string length


//
// Structures
//

// field structure in .html
typedef struct {
  char * name;
  char * value;
} htmlFieldStrct;


//
// Functions
//

// This function extracts a segment from a string between "startTag[]" and "endTag[]".
// The "preTag[]" is on the upstream of and used to specify without ambiguity on "startTag[]". 
// Tags can be set to NULL if not needed and are excluded from the returning segment. 
char *getSegment(const char *string, const char *preTag, const char *startTag, const char *endTag);

// This function gets rid of any space at the beginning or end of a string.
void getRidOfSpace(char *string); 

// This function gets consecutive sections of a string, separated by "sep".
char **getStrSects(int *sectNum, const char *string, char sep);

// This function gets the queryString passed by CGI.
char *getQueryString(void);

// This function gets a field value, specified by its corresponding name, from a queryString.
char *getHtmlFieldValue(const char *fieldName, const char *queryString);

// This function gets the length of a htmlString when all tags are removed.
int getHtmlStrLngth(const char *htmlString);

// This function removes any html tag.
void rmHtmlTag(char *string);

// This function reads from a "file" a "char ** stringList".
// It gets rid of blank lines. If "getComment" is set to 0, it also gets rid of any comments,
// which are proceeded by "#". Then, one must escape "#" as "\#" to have "#" within a string.
// The function also returns the total "int *stringNum" by point. 
char ** getDataStrings(const char * file, int * stringNum, int getComment);


#endif /* _ASAPRATIO_TXTFNS_H_ */
