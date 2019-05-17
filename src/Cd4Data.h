#ifndef _CPPTEST_CD4DATA_H_
#define _CPPTEST_CD4DATA_H_

#include "read_array.h"
#include "consts.h"

struct Cd4Data {
	array_3d progression;
	array_3d initialDist;
	array_3d mortality;
	array_3d paedSurvDist;

	Cd4Data(array_3d cd4Prog, array_3d cd4InitDist, array_3d cd4Mort,
	        array_3d paedSurvCd4Distrib)
		: progression(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES - 1]),
		  initialDist(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  mortality(boost::extents[SEXES][AGE_GROUPS][CD4_STAGES]),
		  paedSurvDist(boost::extents[PROJECTION_YEARS][SEXES][CD4_STAGES]) {

		progression = cd4Prog;
		initialDist = cd4InitDist;
		mortality = cd4Mort;
		paedSurvDist = paedSurvCd4Distrib;
	}

};

#endif