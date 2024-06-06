/*
Program       : ASAPRatio
Author        : Xiao-jun Li <xli@systemsbiology.org>
Date          : 09.17.02
SVN Info      : $Id: ASAPRatio_txtFns.cpp 7996 2019-12-25 00:16:42Z real_procopio $

Functions for text handling in ASAPRatio

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "ASAPRatio_txtFns.h"

// This function extracts a segment from a string between "startTag[]" and "endTag[]".
// The "preTag[]" is on the upstream of and used to specify without ambiguity on "startTag[]".
// Tags can be set to NULL if not needed and are excluded from the returning segment. 
char *getSegment(const char *string, const char *preTag, const char *startTag, const char *endTag)
{
  char *segment;
  const char *tmpStr;

  segment = (char *) calloc(strlen(string)+1, sizeof(char));

  // preTag
  if(preTag == NULL) {
    strcpy(segment, string);
  }
  else {
    if((tmpStr = strstr(string, preTag)) == NULL){
      free(segment);
      return NULL;
    }
    else {
      strcpy(segment, tmpStr+strlen(preTag));
    }
  } 
    
  // startTag
  if(startTag != NULL) {
    if((tmpStr = strstr(segment, startTag)) == NULL) {
      free(segment);
      return NULL;
    }
    else {
      memmove(segment, segment+strlen(startTag), strlen(segment+strlen(startTag))+1);
    }
  } 
    
  // endTag
  if(endTag != NULL) {
    if((tmpStr = strstr(segment, endTag)) == NULL){
      free(segment);
      return NULL;
    }
    else{
      segment[strlen(segment)-strlen(tmpStr)] = '\0';
    }
  } 
  
  return segment;
}



// This function gets rid of any space at the beginning or end of a string.
void getRidOfSpace(char *string) 
{
  int lngth = (int)strlen(string);
  int i;

  // start
  for (i = 0; i < lngth; ++i) {
    if (isspace(string[i]) == 0) // not a space
      break;
  }
  strcpy(string, string+i);

  // end
  lngth = (int)strlen(string);
  for (i = lngth-1; i >= 0; --i) {
    if (isspace(string[i]) == 0) // not a space
      break;
  }
  string[i+1] = '\0';

  return;
}


// This function gets consecutive sections of a string, separated by "sep".
char **getStrSects(int *sectNum, char *string, char sep)
{
  int lngth = (int)strlen(string);
  char **sects;
  int startIndx, endIndx;
  int i;

  // sectNum
  *sectNum = 1;
  for (i = 0; i < lngth; ++i) {
    if (string[i] == sep) {
      ++(*sectNum);
    }
  }

  // sects
  sects = (char **) calloc(*sectNum, sizeof(char *));
  *sectNum = 0;
  startIndx = 0;
  for (i = 0; i < lngth; ++i) {
    if (string[i] == sep) {
      endIndx = i;
      sects[*sectNum] = (char *) calloc(endIndx-startIndx+1, sizeof(char));
      strncpy(sects[*sectNum], string+startIndx, endIndx-startIndx);
      getRidOfSpace(sects[*sectNum]);
      ++(*sectNum);
      startIndx = endIndx + 1;
    }
  }
  endIndx = i;
  sects[*sectNum] = (char *) calloc(endIndx-startIndx+1, sizeof(char));
  strncpy(sects[*sectNum], string+startIndx, endIndx-startIndx);
  getRidOfSpace(sects[*sectNum]);
  ++(*sectNum);
  
  return sects;
}


// This function gets the queryString passed by CGI.
char *getQueryString(void)
{
  char x2c2(char *what);

  char *queryStr; 
  char *queryLenStr;
  long queryLngth;
  char *tmpStr;
  int indx, queryIndx;

  // get queryStr
  if(strcmp(getenv("REQUEST_METHOD"), "GET") == 0) {
    tmpStr = getenv("QUERY_STRING");
    if(tmpStr && ((queryLngth = (int)strlen(tmpStr)) > 0)){
      queryStr = (char *) calloc(queryLngth+1, sizeof(char));
      strcpy(queryStr, tmpStr);
      queryStr[queryLngth] = '\0';
    }
    else
      return NULL;
  } 
  else if(strcmp(getenv("REQUEST_METHOD"), "POST") == 0
	  && (queryLenStr = getenv("CONTENT_LENGTH")) != NULL 
	  && sscanf(queryLenStr, "%ld", &queryLngth) == 1
	  && queryLngth > 0) {
    queryStr = (char *) calloc(queryLngth+1, sizeof(char));
    size_t r = fread(queryStr, sizeof(char), queryLngth, stdin);
    queryStr[queryLngth] = '\0';
  }
  else {
    return NULL;
  }
  //printf("DDS: %s\n", queryStr);

  // decoding
  indx = 0;
  queryIndx = 0;
  while(indx < queryLngth) {
    if(queryStr[indx] == '+')
      queryStr[queryIndx] = ' ';
    else if(queryStr[indx] == '%') {
      queryStr[queryIndx] = x2c2(queryStr+indx+1);
      indx +=2; 
    }     
    else
      queryStr[queryIndx] = queryStr[indx];
    ++indx;
    ++queryIndx;
  } //while(indx <= queryLngth) {
  queryStr[queryIndx] = '\0';
  
  return queryStr;
}

/*
  This function if for decoding queryString.
*/
char x2c2(char *what)
{
  char digit;
  digit = what[0] >= 'A' ? ((what[0] & 0xdf)-'A')+10 : what[0]-'0';
  digit *= 16;
  digit += what[1] >= 'A' ? ((what[1] & 0xdf)-'A')+10 : what[1]-'0';
  return digit;
}


// This function gets a field value, specified by its corresponding name, from a queryString.
char *getHtmlFieldValue(const char *fieldName, const char *queryString)
{
  char *fieldValue;
  const char *tmpValue;
  int qLngth = (int)strlen(queryString);
  int fLngth = (int)strlen(fieldName);
  int lngth;
  int strIndx = 0;
  int i;

  while((tmpValue = strstr(queryString+strIndx, fieldName)) != NULL) { 
    lngth = (int)strlen(tmpValue);
    if ( (tmpValue == queryString || *(tmpValue-1) == '&') ) { // don't mistake "cidIndx" for "Indx"
      for (i = fLngth; i < lngth; ++i){
	if(isspace(tmpValue[i]) == 0){
	  if(tmpValue[i] == '='){ // match
	    if((fieldValue = getSegment(tmpValue, fieldName, "=", "&")) != NULL
	       || (strchr(tmpValue+i+1, '=') == NULL
		   && (fieldValue = getSegment(tmpValue, fieldName, "=", NULL)) != NULL)){
	      getRidOfSpace(fieldValue); 	    
	      if(strlen(fieldValue) > 0)
		return fieldValue;
	      else {
		free(fieldValue);
		return NULL;
	      }
	    } // if((fieldValue = getSegment(tmpValue, fieldName, "=", "&")) != NULL
	    else 
	      return NULL;
	  } // if(tmpValue[i] == '='){ // match
	  else {
	    break;
	  }
	} // if(isspace(tmpValue[i]) == 0){
      } // for (i = strlen(fieldName); i < strlen(tmpValue); ++i){
    } //  if ( (tmpValue == queryString || *(tmpValue-1) == '&') )
    strIndx = qLngth - lngth + fLngth;
  } //   while((tmpValue = strstr(queryString+strIndx, fieldName)) != NULL) { 

  return NULL;
}


// This function removes any html tag.
void rmHtmlTag(char *string)
{
  char *tmpString;
  char *tag;
  int lngth;

  tmpString = (char *) calloc(strlen(string)+1, sizeof(char));
  while((tag = strchr(string, '<')) != NULL){
    // front tag
    lngth = (int)(strlen(string) - strlen(tag));
    strncpy(tmpString, string, lngth); 
    tmpString[lngth] = '\0';
    // end tag
    strcpy(string, tag);
    if((tag = strchr(string, '>')) != NULL){
      lngth += (int)strlen(tag+1);
      strcat(tmpString, tag+1);
      tmpString[lngth] = '\0';
    }
    strcpy(string, tmpString);
  }
  free(tmpString);

  getRidOfSpace(string);

  return;
}

// This function gets the length of a htmlString when all tags are removed.
int getHtmlStrLngth(const char *htmlString)
{
  char *tmpString;
  int lngth = (int)strlen(htmlString);
  
  tmpString = (char *) calloc(lngth+1, sizeof(char));
  strcpy(tmpString, htmlString);
  rmHtmlTag(tmpString); 
  lngth = (int)strlen(tmpString);
  free(tmpString);

  return lngth;
}


// This function reads from a "file" a "char ** stringList". It gets rid of blank lines.
// If "getComment" is set to 0, it also gets rid of any comments, which are proceeded by "#".
// Then, one must escape "#" as "\#" to have "#" within a string. 
// The function also returns the total "int *stringNum" by point. 
char ** getDataStrings(const char * file, int * stringNum, int getComment)
{
  typedef struct get_data_strings_Strct { // to create a linked list of strings
    char * string;  // string
    struct get_data_strings_Strct * next; // link
  } get_data_strings_LinkedList;
  int get_data_strings_IsEmptyString(char * string, int getComment);

  FILE * fin; // pt for input datafile
  get_data_strings_LinkedList *dataList, *currentData;
  char tempString[_MXSTRLEN_];
  char ** stringList;
  long lngth;
  int index;

  // open file
  if((fin = fopen(file,"r")) == NULL) {
    printf("Error in reading \"%s\"! \n", file);
    exit(0);
  }

  // initialize
  *stringNum = 0;
  dataList = (get_data_strings_LinkedList *) 
    calloc(1, sizeof(get_data_strings_LinkedList));  
  currentData = dataList;
  // read string list
  while (!feof(fin)) {
    if (fgets(tempString, _MXSTRLEN_, fin) != NULL 
	&& get_data_strings_IsEmptyString(tempString, getComment) != 1) {
      lngth = (int)strlen(tempString) + 1;
      currentData->string
	= (char *) calloc(lngth, sizeof(char));      
      strcpy(currentData->string, tempString);
      ++ (*stringNum);
      currentData->next = (get_data_strings_LinkedList *) 
	calloc(1, sizeof(get_data_strings_LinkedList));  
      currentData = currentData->next;
    } // if (fgets(tempString, GET_DATA_STRINGS_MXSTRLEN, fin) != NULL 
  } // while (!feof(fin)) {
  currentData->next = NULL;
  fclose(fin);

  // convert linked list into array of strings
  stringList = (char **) calloc(*stringNum, sizeof(char *));
  index = 0; // string index
  while(index < *stringNum) {
    lngth = (int)strlen(dataList->string) + 1;
    stringList[index] = (char *) calloc(lngth, sizeof(char));
    strcpy(stringList[index], dataList->string);
    ++ index;
    currentData = dataList->next;
    free(dataList->string);
    free(dataList);
    dataList = currentData;
  }
  free(dataList);

  return stringList;
}

// This function checks whether or not a string is empty, 
//  or a blank line. If yes, it returns 1. If "getComment" is set
//  to 0, it also gets rid of comments, if any. 
int get_data_strings_IsEmptyString(char * string, int getComment)
{
  int lngth = (int)strlen(string);
  int indx = 0;
  int i;

  if (lngth == 0) return 1; // empty

  if(getComment == 0 && string[0] == '#') return 1; // comment

  // get rid of any space at the end
  for (i = lngth-1; i >= 0; --i) {
    if (isspace(string[i]) != 0) // space
      string[i] = '\0';
    else
      break;
  }

  // check for blank line
  lngth = (int)strlen(string);
  if (lngth == 0) // blank
    return 1;
    
  // for non-empty string, get rid of any comment 
  if(getComment == 0) {
    for (i = 1; i < lngth; ++i) {
      if (string[i] == '#') {
	if (string[i-1] == '\\') { // escape format for '#'
	  strcpy(string+i-1, string+i);
	  --i;
	  --lngth;
	}
	else { // comment
	  lngth = i;
	  string[i] = '\0';
	  break;
	}
      }
    }
  }
  return 0;
}
