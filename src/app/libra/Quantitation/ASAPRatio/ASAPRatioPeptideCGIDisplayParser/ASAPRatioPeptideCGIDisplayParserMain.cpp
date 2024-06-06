/*
Program       : ASAPRatio
Author        : Xiao-jun Li <xli@systemsbiology.org>
Date          : 09.17.02
SVN Info      : $Id: ASAPRatioPeptideCGIDisplayParserMain.cpp 8445 2021-04-20 01:01:01Z real_procopio $

CGI program for displaying and modifying ASAPRatio peptide relative abundance

Copyright (C) 2002 Xiao-jun Li
UI enhancements by L.Mendoza (2006,2020)

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
#include <math.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <sys/param.h>
#endif
#include <sys/stat.h> 
#include "Common/constants.h"
#include "Common/util.h"

#include "Common/TPPVersion.h" // contains version number, name, revision

#include "Quantitation/ASAPRatio/ASAPRatio_Fns/ASAPRatio_txtFns.h"
#include "Quantitation/ASAPRatio/ASAPRatio_Fns/ASAPRatio_numFns.h"
#include "Quantitation/ASAPRatio/ASAPRatio_Fns/ASAPRatio_pepFns.h"
#include "ASAPRatioPeptideCGIDisplayParser.h"
#include "ASAPRatioPeptideUpdateParser.h"
#include "Parsers/Parser/Tag.h"
#include "Common/util.h"

using namespace mzParser;

Array<Tag*>* generateXML(pepDataStrct data_, int index, Boolean cover, Boolean data) {
  Array<Tag*>* output = new Array<Tag*>;

  if(data_.indx == -1)
    return output; // done

  Tag* next;
  char text[100];
  const char* lcpeakNames[] = {"asapratio_lc_lightpeak", "asapratio_lc_heavypeak"};

  if(cover) {
    if(data)
      next = new Tag("asapratio_result", True, False);
    else
      next = new Tag("asapratio_result", True, True);
    // fill it up
    if(index) {
      sprintf(text, "%d", index);
      next->setAttributeValue("index", text);
    }
    if(data_.pepRatio[0] == -2) {
      next->setAttributeValue("mean", "-1");
      sprintf(text, "%0.2f", data_.pepRatio[1]);
      next->setAttributeValue("error", text);
      next->setAttributeValue("heavy2light_mean", "-1");
      next->setAttributeValue("heavy2light_error", "-1");
    }
    else if(data_.pepRatio[0] == -1) {
      next->setAttributeValue("mean", "999.");
      sprintf(text, "%0.2f", data_.pepRatio[1]);
      next->setAttributeValue("error", text);
      next->setAttributeValue("heavy2light_mean", "0.00");
      next->setAttributeValue("heavy2light_error", "0.00");
    }
    else {
      sprintf(text, "%0.2f", data_.pepRatio[0]);
      next->setAttributeValue("mean", text);
      sprintf(text, "%0.2f", data_.pepRatio[1]);
      next->setAttributeValue("error", text);
      sprintf(text, "%0.2f", data_.pepH2LRatio[0]);
      next->setAttributeValue("heavy2light_mean", text);
      sprintf(text, "%0.2f", data_.pepH2LRatio[1]);
      next->setAttributeValue("heavy2light_error", text);
    }

    //next->write(cout);
    output->insertAtEnd(next);

    // now the light mass
    //    sprintf(text, "%f", data_.msLight);
    //    next->setAttributeValue("light_mass", text);

  }
  if(data) {
    next = new Tag("asapratio_peptide_data", True, False);
    if(index) {
      sprintf(text, "%d", index);
      next->setAttributeValue("index", text);
    }
    sprintf(text, "%d", data_.indx);
    next->setAttributeValue("status", text);
    sprintf(text, "%d", data_.cidIndx);
    next->setAttributeValue("cidIndex", text);
    sprintf(text, "%0.4f", data_.msLight);
    next->setAttributeValue("light_mass", text);
    sprintf(text, "%0.4f", data_.msHeavy);
    next->setAttributeValue("heavy_mass", text);
    sprintf(text, "%d", data_.areaFlag);
    next->setAttributeValue("area_flag", text);
    output->insertAtEnd(next);
    for(int k = 0; k < _ASAPRATIO_MXQ_; k++) {
      next = new Tag("asapratio_contribution", True, False);
      sprintf(text, "%0.4f", data_.pkRatio[k]);
      next->setAttributeValue("ratio", text);
      sprintf(text, "%0.4f", data_.pkError[k]);
      next->setAttributeValue("error", text);
      sprintf(text, "%d", k+1);
      next->setAttributeValue("charge", text);
      sprintf(text, "%d", data_.pkCount[k]);
      next->setAttributeValue("use", text);
      output->insertAtEnd(next);
      for(int j = 0; j < 2; j++) {
	next = new Tag(lcpeakNames[j], True, True);
	sprintf(text, "%d", data_.peaks[k][j].indx);
	next->setAttributeValue("status", text);
	sprintf(text, "%d", data_.peaks[k][j].valley[0]);
	next->setAttributeValue("left_valley", text);
	sprintf(text, "%d", data_.peaks[k][j].valley[1]);
	next->setAttributeValue("right_valley", text);
	sprintf(text, "%0.2e", data_.peaks[k][j].bckgrnd);
	next->setAttributeValue("background", text);
	sprintf(text, "%0.2e", data_.peaks[k][j].area[0]);
	next->setAttributeValue("area", text);
	sprintf(text, "%0.2e", data_.peaks[k][j].area[1]);
	next->setAttributeValue("area_error", text);
	sprintf(text, "%0.4f", data_.peaks[k][j].time[0]);
	next->setAttributeValue("time", text);
	sprintf(text, "%0.4f", data_.peaks[k][j].time[1]);
	next->setAttributeValue("time_width", text);
	sprintf(text, "%d", data_.peaks[k][j].peak);
	next->setAttributeValue("is_heavy", text);
	output->insertAtEnd(next);
      } // light/heavy
      next = new Tag("asapratio_contribution", False, True);
      output->insertAtEnd(next);
    } // next precursor ion charge

    next = new Tag("asapratio_peptide_data", False, True);
    output->insertAtEnd(next);

    if(cover) {
      next = new Tag("asapratio_result", False, True);
      output->insertAtEnd(next);
    } // if also cover
  } // if data
  return output;
}


// This function gets a pepDataStrct from a queryString.
pepDataStrct *getPepDataStrctFromQueryString(char *queryString)
{
  pepDataStrct *peptide;
  char *tmpValue;
  char tmpField[100];
  int i, j, k;

  peptide = (pepDataStrct *) calloc(1, sizeof(pepDataStrct));

  // indx
  if ((tmpValue = getHtmlFieldValue("peptide_indx", queryString)) != NULL) {
    sscanf(tmpValue, "%d", &(peptide->indx));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // scan
  if ((tmpValue = getHtmlFieldValue("peptide_scan", queryString)) != NULL) {
    sscanf(tmpValue, "%ld", &(peptide->scan));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // chrg
  if ((tmpValue = getHtmlFieldValue("peptide_chrg", queryString)) != NULL) {
    sscanf(tmpValue, "%d", &(peptide->chrg));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // cidIndx
  if ((tmpValue = getHtmlFieldValue("peptide_cidIndx", queryString)) 
      != NULL) {
    sscanf(tmpValue, "%d", &(peptide->cidIndx));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // msLight
  if ((tmpValue = getHtmlFieldValue("peptide_msLight", queryString)) 
      != NULL) {
    sscanf(tmpValue, "%lf", &(peptide->msLight));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // msHeavy
  if ((tmpValue = getHtmlFieldValue("peptide_msHeavy", queryString)) 
      != NULL) {
    sscanf(tmpValue, "%lf", &(peptide->msHeavy));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // eltn
  if ((tmpValue = getHtmlFieldValue("peptide_eltn", queryString)) != NULL) {
    sscanf(tmpValue, "%d", &(peptide->eltn));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // areaFlag
  if ((tmpValue = getHtmlFieldValue("peptide_areaFlag", queryString)) 
      != NULL) {
    if(strcmp(tmpValue, "raw") == 0) 
      peptide->areaFlag = 1;
    else if(strcmp(tmpValue, "fitting") == 0) 
      peptide->areaFlag = 2;
    else
      peptide->areaFlag = 0;
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  // peaks
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {
    for (j = 0; j < 2; ++j) {
      // indx
      sprintf(tmpField, "peptide_peaks_%d_%d_indx", i, j);
      if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
	sscanf(tmpValue, "%d", &(peptide->peaks[i][j].indx));
	free(tmpValue);
      }
      else {
	free(peptide);
	return NULL;
      }
      // peak
      sprintf(tmpField, "peptide_peaks_%d_%d_peak", i, j);
      if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
	sscanf(tmpValue, "%d", &(peptide->peaks[i][j].peak));
	free(tmpValue);
      }
      else {
	free(peptide);
	return NULL;
      }
      // valley
      for(k = 0; k < 2; ++k) {
	sprintf(tmpField, "peptide_peaks_%d_%d_valley_%d", i, j, k);
	if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
	  sscanf(tmpValue, "%d", &(peptide->peaks[i][j].valley[k]));
	  free(tmpValue);
	}
	else {
	  free(peptide);
	  return NULL;
	}
      }
      // bckgrnd
      sprintf(tmpField, "peptide_peaks_%d_%d_bckgrnd", i, j);
      if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
	sscanf(tmpValue, "%lf", &(peptide->peaks[i][j].bckgrnd));
	free(tmpValue);
      }
      else {
	free(peptide);
	return NULL;
      }
      // area
      for(k = 0; k < 2; ++k) {
	sprintf(tmpField, "peptide_peaks_%d_%d_area_%d", i, j, k);
	if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
	  sscanf(tmpValue, "%lf", &(peptide->peaks[i][j].area[k]));
	  free(tmpValue);
	}
	else {
	  free(peptide);
	  return NULL;
	}
      }
      // time
      for(k = 0; k < 2; ++k) {
	sprintf(tmpField, "peptide_peaks_%d_%d_time_%d", i, j, k);
	if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
	  sscanf(tmpValue, "%lf", &(peptide->peaks[i][j].time[k]));
	  free(tmpValue);
	}
	else {
	  free(peptide);
	  return NULL;
	}
      }
    }
  } //   for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {

  // pkRatio
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {
    sprintf(tmpField, "peptide_pkRatio_%d", i);
    if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
      sscanf(tmpValue, "%lf", &(peptide->pkRatio[i]));
      free(tmpValue);
    }
    else {
      free(peptide);
      return NULL;
    }
  }

  // pkError
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {
    sprintf(tmpField, "peptide_pkError_%d", i);
    if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
      sscanf(tmpValue, "%lf", &(peptide->pkError[i]));
      free(tmpValue);
    }
    else {
      free(peptide);
      return NULL;
    }
  }

  // pkCount
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {
    sprintf(tmpField, "peptide_pkCount_%d", i);
    if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
      sscanf(tmpValue, "%d", &(peptide->pkCount[i]));
      free(tmpValue);
    }
    else {
      free(peptide);
      return NULL;
    }
  }

  // pepRatio
  for (i = 0; i < 2; ++i) {
    sprintf(tmpField, "peptide_pepRatio_%d", i);
    if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
      sscanf(tmpValue, "%lf", &(peptide->pepRatio[i]));
      free(tmpValue);
    }
    else {
      free(peptide);
      return NULL;
    }
  }

  // pepTime
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      sprintf(tmpField, "peptide_pepTime_%d_%d", i, j);
      if ((tmpValue = getHtmlFieldValue(tmpField, queryString)) != NULL) {
	sscanf(tmpValue, "%lf", &(peptide->pepTime[i][j]));
	free(tmpValue);
      }
      else {
	free(peptide);
	return NULL;
      }
    }
  }
  
  // pepArea
  if ((tmpValue = getHtmlFieldValue("peptide_pepArea", queryString)) 
      != NULL) {
    sscanf(tmpValue, "%lf", &(peptide->pepArea));
    free(tmpValue);
  }
  else {
    free(peptide);
    return NULL;
  }

  return peptide;
}


// This function displays a pepDataStrct in cgi program.
void displayPepDataStrctInCgi(pepDataStrct peptide, char *cgiAction, 
			      char *xmlFile, char *baseName,
			      int index, int ratioType,
			      double *accRatio, char *pngFileLink, 
			      char *pngFilePath, char* spectrumName) {
  void displayPepDataStrctInCgi(pepDataStrct peptide, char *cgiAction, 
				char *xmlFile, char *baseName,
				int index, int ratioType,
				double *accRatio, char *pngFileLink, 
				char *pngFilePath, char* spectrumName, 
				Boolean quantHighBG, Boolean zeroBG, double mzBound);
  displayPepDataStrctInCgi( peptide, cgiAction, 
			    xmlFile, baseName,
			    index,  ratioType,
			    accRatio, pngFileLink, 
			    pngFilePath,  spectrumName, 
			    False,  False, -1);
}

void displayPepDataStrctInCgi(pepDataStrct peptide, char *cgiAction, 
			      char *xmlFile, char *baseName,
			      int index, int ratioType,
			      double *accRatio, char *pngFileLink, 
			      char *pngFilePath, char* spectrumName, 
			      Boolean quantHighBG, Boolean zeroBG, double mzBound)
{  
  void displayPepDataStrctInCgi(pepDataStrct peptide, char *cgiAction, 
				char *xmlFile, char *baseName,
				int index, int ratioType,
				double *accRatio, char *pngFileLink, 
				char *pngFilePath, char* spectrumName, 
				Boolean quantHighBG, Boolean zeroBG, double mzBound, bool wavelet);
  displayPepDataStrctInCgi( peptide, cgiAction, 
			    xmlFile, baseName,
			    index,  ratioType,
			    accRatio, pngFileLink, 
			    pngFilePath,  spectrumName, 
			    quantHighBG,  zeroBG, mzBound, false);
  
}

void displayPepDataStrctInCgi(pepDataStrct peptide, char *cgiAction, 
			      char *xmlFile, char *baseName,
			      int index, int ratioType,
			      double *accRatio, char *pngFileLink, 
			      char *pngFilePath, char* spectrumName, 
			      Boolean quantHighBG, Boolean zeroBG, double mzBound, bool wavelet)
{
  const int cgiType = 2; // set to 1 to use gnuplot/png files; 2 uses flot.js
  char *mzXMLFile;
  char **ratioStrings;
  char tmpField[1000];
  double tmpRatio[2];
  time_t currTime;
  int i, j, k;
  int randn;

  /*
   * create .png files
   */

  // pngFileBase
  (void) time(&currTime);
  srand48((long) currTime);
  randn = lrand48();
  sprintf(pngFilePath+strlen(pngFilePath), "_%d", randn);
  sprintf(pngFileLink+strlen(pngFileLink), "_%d", randn);

  // plot data
  int len;
  mzXMLFile = (char *) calloc(len=strlen(baseName)+strlen(xmlFile)+10, sizeof(char));
  rampConstructInputPath(mzXMLFile,len,xmlFile,baseName); // .mzXML or .mzData
  // JS: unimplemented function in mzParser
  // rampValidateOrDeriveInputFilename(mzXMLFile,len,spectrumName);

  if(peptide.indx >= 0) { // if valid data
    getPepDataStrct(&peptide, mzXMLFile, pngFilePath, cgiType, quantHighBG, zeroBG, mzBound, wavelet);
  }
  free(mzXMLFile);


  // form
  printf("<form method='POST' action='%s'>\n\n", cgiAction);
  // hidden "default" submit button - gets submitted when user presses 'enter' anywhere on form
  // does not work on IE.  Too bad...
  printf("<div class='hideit'><input name='submit' value='Evaluate_Ratio' type='submit'/></div>\n");


  printf("<div class='tppbanner' banner-bg-text='ASAPRatio'>&nbsp;&nbsp;ASAPRatio Peptide :: <b class='tppOrange'>%s</b>\n<br/>\n",spectrumName);
  printf("</div>\n");


  // display manual
  printf("<div class='ratiobox tppgraytitle'>\n");
  fflush(stdout);

  // interim ratio and run parameters
  ratioStrings = ASAPRatio_ratioOutput(peptide.pepRatio, ratioType);
  printf("<span class='fields'>");
  printf("<span style='font-size:large;'>Interim Ratio ");
  if(ratioType == 0) 
    printf("(L/H)");
  else if(ratioType == 1) 
    printf("(H/L)");
  printf(": </span><br>\n");

  printf("QuantHighBG :<br>\n");
  printf("ZeroBG :<br>\n");
  printf("mzBound :<br>\n");
  printf("</span>\n");

  printf("<span class='values'>");
  printf("<span style='font-size:large;'><b>%s &plusmn; %s</b> (%s%%)</span><br>\n", ratioStrings[0], ratioStrings[1], ratioStrings[2]);
  if (quantHighBG)
    printf("True<br>\n");
  else
    printf("False<br>\n");

  if (zeroBG)
    printf("True<br>\n");
  else
    printf("False<br>\n");

  printf("%f<br>\n",mzBound);
  printf("</span>\n");

  fflush(stdout);
  freeMtrx((void **)ratioStrings, 3);


  printf("<span class='controls'>\n");
  printf("Set Accepted Ratio to:<br>\n");
  printf("<input name='submit' value='Interim_Ratio' type='submit'>\n");
  printf("<input style='margin-left:10px;' name='submit' value='0:1' type='submit'>\n");
  printf("<input style='margin-left:10px;' name='submit' value='1:0' type='submit'>\n");
  printf("<input style='margin-left:10px;' name='submit' value='Unknown' type='submit'>\n");
  printf("</span>\n");


  // accepted ratio
  ratioStrings = ASAPRatio_ratioOutput(accRatio, ratioType);
  printf("<span class='ratiocontainer'>");
  printf("Accepted <b>Peptide</b> Ratio ");
  if(ratioType == 0)
    printf("(L/H):<br>");
  else if(ratioType == 1)
    printf("(H/L):<br>");
  //else if(ratioType == 2) 
  //  printf("");
  printf("<span style='font-size:xx-large;' class='tpporange'>");
  printf("<b>%s &plusmn; %s</b> (%s%%)\n", ratioStrings[0], ratioStrings[1], ratioStrings[2]);
  printf("</span></span>\n");

  fflush(stdout);
  freeMtrx((void **)ratioStrings, 3);


  // close div
  printf("</div>\n<br/><br/>\n\n");
  fflush(stdout);


  // display individual peaks

  // nav tabs and animation checkboxes
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    printf("<span class='tab'>\n");
    printf("<input name='UIanim_chg' id='UIanim_chg%d' value='%d' type='checkbox'",i,i);
    if(peptide.pkCount[i] == 1)
      printf(" checked='checked'>");
    else
      printf(">");
    printf("<br/>\n");

    printf("<span class='nav' id='plus%d_nav' onclick='display(\"plus%d\")'>",i+1,i+1);

    if(peptide.peaks[i][0].indx >= 0 ||
       peptide.peaks[i][1].indx >= 0 )
      printf("+%d",i+1);
    else
      printf("[ +%d ]",i+1); // no data

    // add a little asterisk to charge state where CID was made
    if(i+1 == peptide.chrg)
      printf("*");

    printf("</span></span>\n");
  }
  printf("<span id='showall_nav' onclick='displayAll();'>Display All</span>\n");

  printf("<span id='animate_nav'>Animate: \n");
  printf("<input id='UIanim_fast' type='button' onClick=\"anim_ctrl('fast');\" value='Fast'>\n");
  printf("<input id='UIanim_slow' type='button' onClick=\"anim_ctrl('slow');\" value='Slow'>\n");
  printf("<input id='UIanim_stop' type='button' onClick=\"anim_ctrl('stop');\" value='Stop' disabled='disabled'>\n");
  printf("</span><br/>\n");

  // chromatogram table
  printf("<table id='chrtable'>\n");

  // Evaluate Ratio button and areaFlag select
  printf("<tr>\n<td style='text-align:right;padding:5px 30px;'><span style='margin-right:20px;' id='yscale_text'></span>\n");
  printf("<input id='y_peak' type='button' title='Independent scales' onclick=\"adjust_y('peak');\" value='Peak' disabled='disabled'>\n");
  printf("<input id='y_norm' type='button' title='Normalize y-scales' onclick=\"adjust_y('norm');\" value='Norm'>\n");
  printf("<input style='margin-left:70px;' name='submit' value='Evaluate_Ratio' type='submit'/></td>\n");

  printf("<td>Area type: ");
  printf("<select onchange='this.className=\"inputchanged\"' name='peptide_areaFlag'>");
  if(peptide.areaFlag == 1) {
    printf("<option selected>raw</option>");
    printf("<option>average</option>");
    printf("<option>fitting</option>");
  }
  else if(peptide.areaFlag == 2) {
    printf("<option selected>fitting</option>");
    printf("<option>average</option>");
    printf("<option>raw</option>");
  }
  else {
    printf("<option selected>average</option>");
    printf("<option>raw</option>");
    printf("<option>fitting</option>");
  }
  printf("</select></td>\n</tr>\n\n");
  fflush(stdout);

  for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {
    if(peptide.peaks[i][0].indx >= 0 ||
       peptide.peaks[i][1].indx >= 0 ) {
      // pngFile
      sprintf(tmpField, "%s_%d.png", pngFileLink, i+1);

      printf("<tr class='showit' id='plus%d_tr'>\n<td style='padding-left:20px;'>", i+1);
      printf("<div class='chr'>\n");

      if (cgiType==1)
	printf("<img src='%s'/>\n", makeTmpPNGFileSrcRef(tmpField).c_str());
      else {  // assume 2
	printf("<div style='text-align:right'>Light +%d, m/z <b>%.4f</b></div>\n",i+1,(peptide.msLight+i*_ASAPRATIO_HM_)/(double)(i+1));
	printf("<div id='light%d_chr' style='height:300px;width:700px;'></div>\n",i+1);
	printf("<div style='text-align:right'>Heavy +%d, m/z <b>%.4f</b></div>\n",i+1,(peptide.msHeavy+i*_ASAPRATIO_HM_)/(double)(i+1));
	printf("<div id='heavy%d_chr' style='height:300px;width:700px;'></div>\n",i+1);
	printf("<script type='text/javascript' language='javascript'>\n");
	if (peptide.chrg == i+1)
	  printf("defaultcharge=%d;\n",i+1);

	printf("if (raw_data_%d && fit_data_%d) {\n",i+1,i+1);
	printf("$(document).ready(function() {\n");
	for (j = 0; j < 2; ++j) {
	  string lh;
	  if(j == 0)
	    lh = "light" + to_string(i+1);
	  else
	    lh = "heavy" + to_string(i+1);

	  printf("allplots['%s_plot'] = $.plot($('#%s_chr'), [\n",lh.c_str(),lh.c_str());

	  printf("{ data: raw_data_%d[%d], color: '#ff0000', label:'Raw', lines: { show:true }, points: { show:true } },\n",i+1,j);
	  printf("{ data: fit_data_%d[%d], color: '#006dff', label:'Fit', lines: { show:true }, points: { show:false} },\n",i+1,j);
	  printf("{ data: peaks_data_%d[%d], color: '#22eceb', label:'Area', lines: { show:true, lineWidth: 1 }, points: { show:false} },\n",i+1,j);
	  printf("],\n");
	  printf("{\n");
	  printf("legend: { position:'nw' },\n");
	  printf("selection: { mode:'x', color: '#ff5f00' },\n");
	  printf("yaxis:  { min: 0, tickFormatter: function (val, axis, precision) { return val.toExponential(1); } },\n");
	  printf("grid:   { markings: [\n");
	  printf("          { color: '#e5ffff', xaxis:{from:%d, to:%d} },\n",peptide.peaks[i][j].valley[0],peptide.peaks[i][j].valley[1]);
	  printf("          { color: '#e0e0e0', yaxis:{from:%f, to:0 } },\n",peptide.peaks[i][j].bckgrnd);
	  printf("          { color: '#82003b', yaxis:{from:%f, to:%f} },\n",peptide.peaks[i][j].bckgrnd,peptide.peaks[i][j].bckgrnd);
	  if (peptide.chrg == i+1 && peptide.cidIndx == j) // CID
	    printf("          { color: '#ff5f00', xaxis:{from:%d, to:%d} }\n",peptide.scan,peptide.scan);
	  else
	    printf("          { color: '#909090', xaxis:{from:%d, to:%d} }\n",peptide.scan,peptide.scan);
	  printf("         ] }\n");

	  printf("});\n");

	  printf("$('#%s_chr').bind( 'plotselected',  function( event, ranges ) { if (mode == 'bg') bg_level('%s', ranges.yaxis.to); else scan_range('%s', ranges.xaxis.from, ranges.xaxis.to, true); } );\n",lh.c_str(),lh.c_str(),lh.c_str());
	  printf("$('#%s_chr').bind( 'plotselecting', function( event, ranges ) { if(ranges){ if (mode == 'bg') bg_level('%s', ranges.yaxis.to); else scan_range('%s', ranges.xaxis.from, ranges.xaxis.to, true); } } );\n",lh.c_str(),lh.c_str(),lh.c_str());
	  printf("$('#%s_chr').bind( 'plotunselected',function( event, ranges ) { if (mode == 'bg') { $('#%s_bg').val(%.2e); $('#%s_bg').removeClass(\"inputchanged\"); } else { $('#%s_from').val(%d); $('#%s_to').val(%d); $('#%s_from').removeClass(\"inputchanged\"); $('#%s_to').removeClass(\"inputchanged\");} } );\n",
		 lh.c_str(),lh.c_str(),peptide.peaks[i][j].bckgrnd,lh.c_str(),lh.c_str(),peptide.peaks[i][j].valley[0],lh.c_str(),peptide.peaks[i][j].valley[1],lh.c_str(),lh.c_str());

	}
	printf("});\n}\n");
	printf("</script>\n");
      }

      printf("</div></td>\n");


      // table
      // charge state
      printf("<td valign='top' id='plus%d' ", i+1);
      if(peptide.pkCount[i] == 1)
	printf("class='accepted'");
      else
	printf("class='rejected'");
      printf(">\n<table class='chgtable'>\n<tr><th style='padding: 0px;' colspan='2'><span class='banner_chg'>+%d</span></th></tr>\n", i+1);

      // ratio
      tmpRatio[0] = peptide.pkRatio[i];
      tmpRatio[1] = peptide.pkError[i];

      ratioStrings = ASAPRatio_ratioOutput(tmpRatio, ratioType);
      printf("<tr class='banner1'><td>Ratio:</td>");
      printf("<td class='ratio'>%s &plusmn; %s (%s%%)</td></tr>\n", 
	     ratioStrings[0], ratioStrings[1], ratioStrings[2]);
      freeMtrx((void **)ratioStrings, 3);

      // weight
      printf("<tr><td>Weight:</td>");
      printf("<td>%g</td></tr>\n", 
	     peptide.peaks[i][0].area[0]+peptide.peaks[i][1].area[0]);

      // acceptance
      printf("<tr><td><b>Acceptance:</b></td>");
      printf("<td><label><input type='radio' name='peptide_pkCount_%d' onClick='highlight(\"plus%d\",\"yes\")' value='1'", i, i+1);

      if(peptide.pkCount[i] == 1)
	printf(" checked='checked'");

      printf("/><span>Yes</span></label> <label><input type='radio' name='peptide_pkCount_%d' onClick=\"highlight('plus%d','no')\" value='0'", i, i+1);

      if(peptide.pkCount[i] != 1)
	printf(" checked='checked'");

      printf("/><span>No</span></label></td>");
      printf("</tr>\n");

      printf("<tr><td><br><br></td><td></td></tr>\n"); 

      // light
      if(peptide.chrg == i+1 && peptide.cidIndx == 0)
	printf("<tr><th class='banner_cid' ");
      else
	printf("<tr><th class='banner2' ");

      printf("colspan='2'>Light +%d</th></tr>\n", i+1);      

      // light scan
      printf("<tr><td>Scan Range: </td><td>");
      printf("<input onchange='this.className=\"inputchanged\"' type='text' id='light%d_from' name='peptide_peaks_%d_0_valley_0' size='5' value='%d'/> ",i+1,i, peptide.peaks[i][0].valley[0]);
      printf("<input onchange='this.className=\"inputchanged\"' type='text' id='light%d_to'   name='peptide_peaks_%d_0_valley_1' size='5' value='%d'/></td></tr>\n",i+1,i, peptide.peaks[i][0].valley[1]);

      // light bckgrnd
      printf("<tr><td>Background: </td>");
      printf("<td><input onchange='this.className=\"inputchanged\"' type='text' id='light%d_bg' name='peptide_peaks_%d_0_bckgrnd' size='8' value='%.2e'/><label title='set bg to zero' onclick='bg_level(\"light%d\",0);'>&#10060;</label></td></tr>\n", i+1,i, peptide.peaks[i][0].bckgrnd, i+1);

      // light elution time
      printf("<tr><td>Time: </td><td>%.2f &plusmn; %.2f (min)</td></tr>\n", 
	     peptide.peaks[i][0].time[0], peptide.peaks[i][0].time[1]);

      // id table
      printf("<tr id='light%d_ids'></tr>\n\n",i+1);

      // spacer
      printf("<tr><td><br><br></td><td></td></tr>\n");

      // heavy
      if(peptide.chrg == i+1 && peptide.cidIndx == 1)
	printf("<tr><th class='banner_cid' ");
      else
	printf("<tr><th class='banner2' ");
      printf("colspan='2'>Heavy +%d</th></tr>\n", i+1);      

      // heavy scan
      printf("<tr><td>Scan Range: </td><td>");
      printf("<input onchange='this.className=\"inputchanged\"' type='text' id='heavy%d_from' name='peptide_peaks_%d_1_valley_0' size='5' value='%d'/> ",i+1,i, peptide.peaks[i][1].valley[0]);
      printf("<input onchange='this.className=\"inputchanged\"' type='text' id='heavy%d_to'   name='peptide_peaks_%d_1_valley_1' size='5' value='%d'/></td></tr>\n",i+1,i, peptide.peaks[i][1].valley[1]);

      // heavy bckgrnd
      printf("<tr><td>Background: </td> ");
      printf("<td><input onchange='this.className=\"inputchanged\"' type='text' id='heavy%d_bg' name='peptide_peaks_%d_1_bckgrnd' size='8' value='%.2e'/><label title='set bg to zero' onclick='bg_level(\"heavy%d\",0);'>&#10060;</label></td></tr>\n", i+1,i, peptide.peaks[i][1].bckgrnd,i+1);

      // heavy elution time
      printf("<tr><td>Time: </td><td>%.2f &plusmn; %.2f (min)</td></tr>\n", 
	     peptide.peaks[i][1].time[0], peptide.peaks[i][1].time[1]);

      // id table
      printf("<tr id='heavy%d_ids'></tr>\n\n",i+1);

      printf("</table>\n\n");

      printf("<br><input style='float:right' type='button' title='Retrieve and overlay all peptide IDs in scan range' onclick=\"get_ids(%d,%f,%f,2.0,'%s',this);\" value='Get IDs in Neighborhood'>\n",
	     i+1,(peptide.msLight+i*_ASAPRATIO_HM_)/(double)(i+1),(peptide.msHeavy+i*_ASAPRATIO_HM_)/(double)(i+1),xmlFile);

      printf("</td>\n</tr>\n\n");
    }
    else {
      printf("<tr class='hideit' id='plus%d_tr'>\n", i+1);
      printf("<td class='rejected' style='text-align:center;min-width:640px;'><b>-- No +%d Data --</b></td>\n",i+1);
      printf("<td style='min-width:275px;vertical-align:top' class='rejected' id='plus%d'>\n", i+1);
      printf("<table class='chgtable'>\n<tr><th style='padding: 0px;' colspan='2'><span class='banner_chg'>+%d</span></th></tr>\n", i+1);
      printf("<tr class='banner1'><td>Ratio:</td>");
      printf("<td class='ratio'>No Ratio</td></tr>\n");
      printf("<tr><td>Weight:</td><td>------</td></tr>\n");
      printf("<tr><td><b>Acceptance:</b></td>");
      printf("<td>&nbsp;&nbsp;<b>No</b></td>");
      printf("</tr>\n");
      printf("<tr><td><br/><br/></td><td></td></tr>\n");
      printf("</table>\n</td>\n</tr>\n\n");
    }
    fflush(stdout);
  }

  printf("<tr>\n<td style='text-align:right; padding:5px 30px;'><span style='margin-right:20px;' id='mode_text'></span>\n");
  printf("<input id='ctrl_sr' type='button' title='Adjust Scan Range' onclick=\"adjust_ctrl('sr');\" value='Scan' disabled='disabled'>\n");
  printf("<input id='ctrl_bg' type='button' title='Adjust Background' onclick=\"adjust_ctrl('bg');\" value='Bg'>\n");
  printf("<input style='margin-left:70px;' name='submit' value='Evaluate_Ratio' type='submit'/></td>\n");
  printf("<td>&nbsp;</td>\n</tr>\n");
  printf("</table>\n\n\n");
  fflush(stdout);


  // hidden fields
  printf("<input type='hidden' name='Xmlfile' value='%s' />\n", xmlFile);
  printf("<input type='hidden' name='Spectrum' value='%s' />\n", spectrumName);
  printf("<input type='hidden' name='Basename' value='%s' />\n", baseName);
  printf("<input type='hidden' name='Indx' value='%d' />\n", index);
  printf("<input type='hidden' name='ratioType' value='%d' />\n", ratioType);
  printf("<input type='hidden' name='quantHighBG' value='%d' />\n", (int)quantHighBG);
  printf("<input type='hidden' name='zeroBG' value='%d' />\n", (int)zeroBG);
  printf("<input type='hidden' name='mzBound' value='%f' />\n", mzBound);
  printf("<input type='hidden' name='peptide_pepArea' value='%f' />\n", peptide.pepArea);
  printf("<input type='hidden' name='peptide_indx' value='%d' />\n", peptide.indx);
  printf("<input type='hidden' name='peptide_scan' value='%ld' />\n", peptide.scan);
  printf("<input type='hidden' name='peptide_chrg' value='%d' />\n", peptide.chrg);
  printf("<input type='hidden' name='peptide_cidIndx' value='%d' />\n", peptide.cidIndx);
  printf("<input type='hidden' name='peptide_msLight' value='%f' />\n", peptide.msLight);
  printf("<input type='hidden' name='peptide_msHeavy' value='%f' />\n", peptide.msHeavy);
  printf("<input type='hidden' name='peptide_eltn' value='%d' />\n", peptide.eltn);
  // printf("<input type='hidden' name='AsapIndex' value='%d' />\n", pepIndx);

  for (i = 0; i < 2; ++i) {
    printf("<input type='hidden' name='accRatio_%d' value='%f' />\n", i, accRatio[i]);
    printf("<input type='hidden' name='peptide_pepRatio_%d' value='%f' />\n", i, peptide.pepRatio[i]);
    for (j = 0; j < 2; ++j)
      printf("<input type='hidden' name='peptide_pepTime_%d_%d' value='%f' />\n", i, j, peptide.pepTime[i][j]);
  }
  
  fflush(stdout);

  // peaks
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    for (j = 0; j < 2; ++j) {
      // indx and peak
      printf("<input type='hidden' name='peptide_peaks_%d_%d_indx' value='%d'/>\n", i, j, peptide.peaks[i][j].indx);
      printf("<input type='hidden' name='peptide_peaks_%d_%d_peak' value='%d'/>\n", i, j, peptide.peaks[i][j].peak);
      // valley and bckgrnd
      if(peptide.peaks[i][0].indx < 0 &&
	 peptide.peaks[i][1].indx < 0 ) {
	for(k = 0; k < 2; ++k) 
	  printf("<input type='hidden' name='peptide_peaks_%d_%d_valley_%d' value='%d'/>\n", i, j, k, peptide.peaks[i][j].valley[k]);
	printf("<input type='hidden' name='peptide_peaks_%d_%d_bckgrnd' value='%f'/>\n", i, j, peptide.peaks[i][j].bckgrnd);
      } // if(peptide.peaks[i][0].indx < 0

      // area and time
      for(k = 0; k < 2; ++k) {
	printf("<input type='hidden' name='peptide_peaks_%d_%d_area_%d' value='%f'/>\n", i, j, k, peptide.peaks[i][j].area[k]);
	printf("<input type='hidden' name='peptide_peaks_%d_%d_time_%d' value='%f'/>\n", i, j, k, peptide.peaks[i][j].time[k]);
      }
    }
    fflush(stdout);
  }


  // pkRatio, pkError, and pkCount
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    printf("<input type='hidden' name='peptide_pkRatio_%d' value='%f' />\n", i, peptide.pkRatio[i]);
    printf("<input type='hidden' name='peptide_pkError_%d' value='%f' />\n", i, peptide.pkError[i]);

    if(peptide.peaks[i][0].indx < 0 &&
       peptide.peaks[i][1].indx < 0 ) {
      printf("<input type='hidden' name='peptide_pkCount_%d' value='%d' />\n", i, peptide.pkCount[i]);
    }
    fflush(stdout);
  }

  printf("\n</form>\n\n");

  if (cgiType==1)
    printf("<br/>Graphs are generated by <a target='_new' href='http://www.gnuplot.info'>gnuplot</a>.<br/>\n");

  fflush(stdout);
  return;
}


int main(int argc, char **argv)
{
  hooks_tpp handler(argc,argv); // set up install paths etc
 
  // cgi variables
  char *queryString;
  char cgiAction[1000];
  char szWebserverRoot[1000];
  const char *pStr;
  // parameters
  char *xmlFile=NULL;  // pepXML file
  char *baseName=NULL; // base name for .mzML file
  char *timeStamp; 
  int index;
  int ratioType;
  char* spectrumName=NULL; 
  int quantHighBG = 0;
  int zeroBG = 0;
  int wv = 0;
  bool wavelet = false;
  double mzBound = -1;

  // files
  char pngFilePath[1000]; 
  char pngFileLink[1000];

  pStr=getWebserverRoot();
  if (pStr==NULL) {
    printf("<PRE> Environment variable WEBSERVER_ROOT does not exist.\n\n");
#ifdef __MINGW__
    printf(" For Windows users, you can set this environment variable\n");
    printf(" through the Advanced tab under System Properties when you\n");
    printf(" right-mouse-click on your My Computer icon.\n\n");
    printf(" Set this environment variable to your webserver's document\n");
    printf(" root directory such as c:\\inetpub\\wwwroot for IIS or\n");
    printf(" c:\\website\\htdocs or WebSite Pro.\n\n");
#endif
    printf(" Exiting.\n");
    exit(0);
  }
  else {
    strcpy(szWebserverRoot, pStr); 
  }
  // Check if szWebserverRoot is present
  if (access(szWebserverRoot, F_OK)) {
    printf(" Cannot access the webserver's root directory:\n");
    printf("    %s\n", szWebserverRoot);
    printf(" This was set as the environment variable WEBSERVER_ROOT\n\n");
    printf(" For Windows users, you can check this environment variable\n");
    printf(" through the Advanced tab under System Properties when you\n");
    printf(" right-mouse-click on your My Computer icon.\n\n");
    printf(" Exiting.\n");
    exit(1);
  }


  // variables
  pepDataStrct *peptide=NULL;
  char *mzXMLFile;
  char tmpString[1000];
  int wrtIndx;
  double accRatio[2];
  char *tmpValue;

  int i;

  // html header, style-sheet, and javascript
  printf("Content-type: text/html\n\n");
  printf("<html>\n<head>\n");
  printf("<title>ASAPRatio: Peptide ratio</title>\n");
  printf("<link rel='stylesheet' type='text/css' href='%scss/tpp.css'>\n", getHtmlUrl());
  printf("<link rel='stylesheet' type='text/css' href='%scss/asap.css'>\n", getHtmlUrl());
  printf("<script type='text/javascript' src='%sjs/tpp.js'></script>\n", getHtmlUrl());
  printf("<script type='text/javascript' src='%sjs/jquery.min.js'></script>\n", getHtmlUrl());
  printf("<script type='text/javascript' src='%sjs/jquery.flot.js'></script>\n", getHtmlUrl());
  printf("<script type='text/javascript' src='%sjs/jquery.flot.selection.js'></script>\n", getHtmlUrl());
  printf("<script type='text/javascript' src='%sjs/jquery.flot.errorbars.js'></script>\n", getHtmlUrl());
  printf("<script type='text/javascript' src='%sjs/jquery.flot.axislabels.js'></script>\n", getHtmlUrl());
  printf("<script type='text/javascript' src='%sjs/asap.js'></script>\n", getHtmlUrl());

  printf("<script language='JavaScript'>\n");
  printf("    var charges = new Array();\n");
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i)
    printf("    charges[%d] = 'plus%d';\n",i,i+1);
  printf("\n");

  printf("    function animate(){\n");
  printf("      if (++animi == %d) animi=0; \n", _ASAPRATIO_MXQ_);
  printf("      tt = 0;\n");
  printf("      if (animation && document.getElementById('UIanim_chg'+animi).checked) {\n");
  printf("          display(charges[animi]);\n");
  printf("	    tt = animf;\n");
  printf("      }\n");
  printf("      t = setTimeout('animate();', tt);\n");
  printf("    }\n");

  printf("</script>\n");
  printf("</head>\n\n");

  printf("<body onload='init();'>\n<div id='tppWrapper'>\n");
  printf("<div id='tppFiller'><b>ASAPRatio Peptide</b> (%s) :: <b>loading...</b><br><br>\n<span id='tppSpinner' style='margin-left:150px;'><div class='tppspinner1'></div><div class='tppspinner2'></div><div class='tppspinner3'></div></span>",szTPPVersionInfo);
  fflush(stdout);


  // collect information

  // initial steps

  // get queryString
  queryString = getQueryString();
  if(queryString == NULL) {
    printf("<font color='red'>Error in passing parameters from web.</font><br/>\n");
    printf("</body></html>\n");
    fflush(stdout);
    return 1;
  }

  // get cgiAction
  if((tmpValue = getenv("SCRIPT_NAME")) != NULL) {
    sprintf(cgiAction, "%s", tmpValue);
  }
  else {
    printf("<font color='red'>Cannot find SCRIPT_NAME. </font><br/>\n");
    printf("</body></html>\n");
    fflush(stdout);
    free(queryString);
    return 1;
  }


  // collect parameters from web file

  // collect information from link in XML file 
  if(strcmp(getenv("REQUEST_METHOD"), "GET") == 0) {
    // xmlFile
    if((xmlFile = getHtmlFieldValue("Xmlfile", queryString)) == NULL){
      printf("<font color='red'>No input for pepXML file!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      return 1;
    }
    fixPath(xmlFile,1); // tidy up path sep chars etc - expect existence

    // baseName
    if((baseName = getHtmlFieldValue("Basename", queryString)) == NULL){
      printf("<font color='red'>No input for .mzML file!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      return 1;
    }

    // timeStamp
    if((timeStamp = getHtmlFieldValue("Timestamp", queryString)) == NULL){
      printf("<font color='red'>No input for time stamp!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      free(baseName);
      return 1;
    }

    // spectrumName
    if((spectrumName = getHtmlFieldValue("Spectrum", queryString)) == NULL){
      printf("<font color='red'>No input for spectrum!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      free(baseName);
      return 1;
    }
    // quantHighBG
    if((tmpValue = getHtmlFieldValue("quantHighBG", queryString)) != NULL) {
      sscanf(tmpValue, "%d", &quantHighBG);
      free(tmpValue);

    }
    // zeroBG
    if((tmpValue = getHtmlFieldValue("zeroBG", queryString)) != NULL) {
      sscanf(tmpValue, "%d", &zeroBG);
      free(tmpValue);

    }
    // wavelet
    if((tmpValue = getHtmlFieldValue("wavelet", queryString)) != NULL) {
      sscanf(tmpValue, "%d", &wv);
      free(tmpValue);
      if (wv) {
	wavelet = true;
      }
    }
    // mzBound
    if((tmpValue = getHtmlFieldValue("mzBound", queryString)) != NULL) {
      sscanf(tmpValue, "%lf", &mzBound);
      free(tmpValue);

    }
    // Indx
    if((tmpValue = getHtmlFieldValue("Indx", queryString)) != NULL) {
      sscanf(tmpValue, "%d", &index);
      //cout << "just read in: " << index << endl;

      free(tmpValue);
      if(index < 1) {
	printf("<font color='red'>Invalid peptide index: '%d'!</font><br/>\n", index);
	printf("</body></html>\n");
	fflush(stdout);
	free(queryString);
	free(xmlFile);
	free(baseName);
	free(timeStamp);
	return 1;
      }
    }
    else {
      printf("<font color='red'>No peptide index passed!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      free(baseName);
      free(timeStamp);
      return 1;
    }

    // peptide
    ASAPRatioPeptideCGIDisplayParser* parser = new ASAPRatioPeptideCGIDisplayParser(xmlFile, baseName, timeStamp, index, -2, zeroBG, mzBound);
        
    if(parser == NULL || ! parser->found()) {
      cout << "Error: could not find entry for " << index << " index with basename " << baseName << " in xmlfile: " << xmlFile << endl;
      exit(1);
    }

    peptide = new pepDataStrct();
    *peptide = parser->getPepDataStruct();

    
    // accepted ratio
    accRatio[0] = peptide->pepRatio[0];
    accRatio[1] = peptide->pepRatio[1];
    
    free(timeStamp);
    delete parser;
  } //   if(strcmp(getenv("REQUEST_METHOD"), "GET") == 0) 

  // collect information from CGI
  if(strcmp(getenv("REQUEST_METHOD"), "GET") != 0) {
    // xmlFile
    if((xmlFile = getHtmlFieldValue("Xmlfile", queryString)) == NULL){
      printf("<font color='red'>No input for pepXML file!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      return 1;
    }

    // spectrumName
    if((spectrumName = getHtmlFieldValue("Spectrum", queryString)) == NULL){
      printf("<font color='red'>No input for spectrum!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      return 1;
    }

    // baseName
    if((baseName = getHtmlFieldValue("Basename", queryString)) == NULL){
      printf("<font color='red'>No input for .mzML file!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      return 1;
    }
    /*
    // pepIndx
    if((tmpValue = getHtmlFieldValue("AsapIndex", queryString)) != NULL) {
      sscanf(tmpValue, "%d", &pepIndx);
      free(tmpValue);
      if(pepIndx < 1) {
	printf("<font color='red'>Invalid peptide index: '%d'!</font><br/>\n", pepIndx);
	printf("</body></html>\n");
	fflush(stdout);
	free(queryString);
	free(xmlFile);
	free(baseName);
	return 1;
      }
    }
    */
    // pepIndx
    if((tmpValue = getHtmlFieldValue("Indx", queryString)) != NULL) {
      sscanf(tmpValue, "%d", &index);
      free(tmpValue);
      if(index < 1) {
	printf("<font color='red'>Invalid peptide index: '%d'!</font><br/>\n", index);
	printf("</body></html>\n");
	fflush(stdout);
	free(queryString);
	free(xmlFile);
	free(baseName);
	return 1;
      }
    }
    else {
      printf("<font color='red'>No peptide index passed!</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      free(baseName);
      return 1;
    }

    // peptide
    if((peptide = getPepDataStrctFromQueryString(queryString)) == NULL){
      printf("<font color='red'>Cannot construct pepDataStrct from CGI.</font><br/>\n");
      printf("</body></html>\n");
      fflush(stdout);
      free(queryString);
      free(xmlFile);
      free(baseName);
      return 1;
    }
    
    // accRatio
    for (i = 0; i < 2; ++i) {
      sprintf(tmpString, "accRatio_%d", i);
      if ((tmpValue = getHtmlFieldValue(tmpString, queryString)) != NULL) {
	sscanf(tmpValue, "%lf", &(accRatio[i]));
	free(tmpValue);
      }
      else
	accRatio[i] = peptide->pepRatio[i];
    }
  } //   if(strcmp(getenv("REQUEST_METHOD"), "GET") != 0) 

  // ratioType
  if ((tmpValue = getHtmlFieldValue("ratioType", queryString)) != NULL) {
    if(sscanf(tmpValue, "%d", &ratioType) != 1
       || ratioType < 0 
       || ratioType > 2) {
      ratioType = 0;
    }
    free(tmpValue);
  }
  else
    ratioType = 0;

  // quantHighBG
  if ((tmpValue = getHtmlFieldValue("quantHighBG", queryString)) != NULL) {
    if(sscanf(tmpValue, "%d", &quantHighBG) != 1
       || quantHighBG < 0 
       || quantHighBG > 1) {
      quantHighBG = 0;
    }
    free(tmpValue);
  }
  else
    quantHighBG = 0;

  // zeroBG
  if ((tmpValue = getHtmlFieldValue("zeroBG", queryString)) != NULL) {
    if(sscanf(tmpValue, "%d", &zeroBG) != 1
       || zeroBG < 0 
       || zeroBG > 1) {
      zeroBG = 0;
    }
    free(tmpValue);
  }
  else
    zeroBG = 0;


  // mzBound
  if ((tmpValue = getHtmlFieldValue("mzBound", queryString)) != NULL) {
    if(sscanf(tmpValue, "%lf", &mzBound) != 1
       || mzBound <= 0 
       || mzBound > 1) {
      mzBound = 0.5;
    }
    free(tmpValue);
  }
  else
    mzBound = 0.5;
  
  strcpy(pngFilePath, xmlFile);
  char* tmpStr = findRightmostPathSeperator(pngFilePath); 
  if (tmpStr++)
     *tmpStr = 0;
  strcat(pngFilePath, "ASAPRatio");
  replace_path_with_webserver_tmp(pngFilePath,sizeof(pngFilePath)); // write this in tmpdir if we have one
  strcpy(pngFileLink, pngFilePath);
  translate_absolute_filesystem_path_to_relative_webserver_root_path(pngFileLink); // so /inetpub/wwwroot/foo becomes /foo

  // submit
  if((tmpValue = getHtmlFieldValue("submit", queryString)) != NULL) {
     int len;
    mzXMLFile = (char *) calloc(len=strlen(baseName)+strlen(xmlFile)+10, sizeof(char));
    rampConstructInputPath(mzXMLFile, len, xmlFile, baseName); // .mzML, mzXML or .mzData
    wrtIndx = 1;
    if(strcmp(tmpValue, "Interim_Ratio") == 0) { // accept ratio
      peptide->indx = 2;
      getPepDataStrct(peptide, mzXMLFile, NULL, 0,  quantHighBG, zeroBG, mzBound, wavelet);
    }
    else if(strcmp(tmpValue, "0:1") == 0) { // set to 0:1
      peptide->indx = 2;
      if(ratioType == 1)
	peptide->pepRatio[0] = -1.;
      else
	peptide->pepRatio[0] = 0.;
      peptide->pepRatio[1] = 0.;
    }
    else if(strcmp(tmpValue, "1:0") == 0) { // set to 1:0
      peptide->indx = 2; 
      if(ratioType == 1)
	peptide->pepRatio[0] = 0.;
      else 
	peptide->pepRatio[0] = -1.;
      peptide->pepRatio[1] = 0.;
    }
    else if(strcmp(tmpValue, "Unknown") == 0) { // set to 0:0
      peptide->indx = 2; 
      peptide->pepRatio[0] = -2.;
      peptide->pepRatio[1] = 0.;
    }
    else {
      wrtIndx = 0;
    }
    free(tmpValue);
    free(mzXMLFile);

    // update XML file
    if(wrtIndx == 1) {
      // accepted ratio
      accRatio[0] = peptide->pepRatio[0];
      accRatio[1] = peptide->pepRatio[1];

      ASAPRatioPeptideUpdateParser* updateParser = new ASAPRatioPeptideUpdateParser(xmlFile, baseName, index, generateXML(*peptide, 0, True, True));

      if(updateParser != NULL && updateParser->update()) {
	//cout << "<h5>Changes made to " << xmlFile << ", refresh browser to view</h5>" << endl;
	cout << "<script language='JavaScript'>"
	     << "messages.push('Changes made to  " << xmlFile << "');"
	     << "messages.push('Refresh browser to view updated analysis.');"
	     << "</script>" << endl;
      }
      else {
	//cout << "<h5>Error: no changes written to " << xmlFile << "</h5>" << endl;
	cout << "<script language='JavaScript'>"
	     << "messages.push('Error: No changes written to " << xmlFile << "');"
	     << "</script>" << endl;
      }

    } // if(wrtIndx == 1) {
  } //  if((tmpValue = getHtmlFieldValue("submit", queryString)) != NULL) {

  free(queryString);
  printf("</div>\n"); // filler

  printf("<script language='JavaScript'>\n");
  printf(" peptide.scan = %d;\n",peptide->scan);
  printf("</script>\n");

  // cgi display
  displayPepDataStrctInCgi(*peptide, cgiAction, xmlFile, baseName, index, ratioType, accRatio,
			   pngFileLink, pngFilePath, spectrumName, quantHighBG, zeroBG, mzBound, wavelet);


  printf("<br style='clear:both'>\n<br/><br/><br/><br/>");
  printf("<footer id='tppPageFooter'>ASAPRatio/Peptide Ratio :: <b>%s</b><br/>\n", spectrumName);
  printf("<b>%s</b><br>\n", xmlFile);
  printf("%s<br/></footer>\n", szTPPVersionInfo);
  printf("</div></body>\n</html>");
  fflush(stdout);

  // free memory
  free(xmlFile);
  free(baseName);
  delete peptide;

  exit(0);
}
