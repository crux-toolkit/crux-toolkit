#ifndef PROSPECTOR_RESULT_H
#define PROSPECTOR_RESULT_H

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "SearchResult.h"

/*
Program       : ProspectorResult for discr_calc of PeptideProphet
Author        : Peter Baker <pbaker@cgl.ucsf.edu>
Date          : 11.07.17
SVN Info      : $Id$

Copyright (C) 2017 Peter Baker

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

Peter Baker
Mass Spectrometry Facility
University of California San Francisco
600 16th Street,
Genentech Hall, suite N472A
San Francisco, CA 94143-2240
USA
pbaker@cgl.ucsf.edu

*/

class ProspectorResult : public SearchResult {

public:
	ProspectorResult(char* szBuf, Boolean preexisting_probs);
	ProspectorResult(Array<Tag*>* tags);
	void process(char* szBuf, Boolean preex_probs, char* val);
	const char* getName();

	double expect_;
	double pvalue_;
	double disc_score_;
	double ion_score_;
	double ion_score_diff_;

	protected:
};

#endif
